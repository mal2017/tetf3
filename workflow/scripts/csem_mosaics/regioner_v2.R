Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(nullranges)
library(plyranges)
library(tidyverse)
library(regioneR)
library(excluderanges)
library(AnnotationHub)
library(patchwork)
library(RcppHMM)

# latest version of this script uses the fixed reference insertions file made upstream,
# so we only have to do it once

# chrom sizes for regioneR, nullranges
g <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6
seqlevelsStyle(g) <- "NCBI"
g_r <- unlist(tileGenome(seqlengths = seqlengths(g),ntile = 1))

# get list of tes coex w/ each gene, union M+F
coex_json_fl <- "results/resources/coexpressed_tes.json"
coex_json_fl <- snakemake@input$coex_json
coex_te_names <- jsonlite::read_json(coex_json_fl,simplifyVector = T)$all

# all reference insertions - just used to control for te density for some
# of the regioneR analysis variants
ref_ins_fl <- "upstream/reference_insertions.bed"
ref_ins_fl <- snakemake@input$ref_ins
ref_ins <- read_bed(ref_ins_fl) |>
  filter(!str_detect(name,"\\(.+\\)n"))
ref_ins <- sort(ref_ins)

# populate sequence lenghts field, required for some functions in regioner
csizes <- seqinfo(g) |> as.data.frame() |> DataFrame() |> 
  as_tibble(rownames="seqnames") |>
  dplyr::select(seqnames, seqlengths) |>
  deframe()

seqlengths(ref_ins) <- csizes[seqlevels(ref_ins)]

# get putatively fixed insertions
fixed_ref_ins_fl <- "results/resources/putative_fixed_insertions.rds"
fixed_ref_ins_fl <- snakemake@input$fixed_ref_ins
fixed_ref_ins <- read_rds(fixed_ref_ins_fl)

# -----------------------------------------------------------------------------------------------------------
# get CSEM/mosaics peak sets for each factor of interest
# -----------------------------------------------------------------------------------------------------------
get_replicated_pks <- function(glb, extraction_pattern = "(?<=mosaics/mosaics\\/).+(?=\\/)") {
  # takes a glob pattern matching peak bed files
  # and regex pattern for extracting sample names
  pks <- Sys.glob(glb)
  
  names(pks) <- str_extract(pks, extraction_pattern)  
  
  pks <- imap(pks, ~mutate(read_bed(.x), sample_name=.y)) |>
    GRangesList()
  
  pks <- pks |> 
    unlist() |>
    reduce_ranges() %>%
    plyranges::mutate(n_overlaps = count_overlaps(.,pks)) |>
    filter(n_overlaps > 1)
  
  pks
}

replicated_pks <- list(pan="upstream/csem_mosaics/mosaics/mosaics/pan_*rep*/pan*.mosaics.bed",
     gro="upstream/csem_mosaics/mosaics/mosaics/gro*rep*/gro*rep*.mosaics.bed",
     H3K9Me3="upstream/csem_mosaics/mosaics/mosaics/*_H3K9Me3_ChIPSeq_*/*H3K9Me3_ChIPSeq_*.mosaics.bed",
     Nnk="upstream/csem_mosaics/mosaics/mosaics/CG17802*/*mosaics.bed",
     vvl="upstream/csem_mosaics/mosaics/mosaics/vvl*/*mosaics.bed",
     Odj="upstream/csem_mosaics/mosaics/mosaics/CG7357_*_rep*/*rep*.mosaics.bed",
     NfI="upstream/csem_mosaics/mosaics/mosaics/NfI*rep*/*rep*.mosaics.bed",
     CG16779="upstream/csem_mosaics/mosaics/mosaics/CG16779*rep*/*rep*.mosaics.bed") |>
  map(get_replicated_pks)


# ------------------------------------------------------------------------------
# fixed ref ins for tes coex with each factor
# ------------------------------------------------------------------------------
factor_specific_fixed_ref_ins <- names(replicated_pks) |>
  set_names() |>
  map(~{fixed_ref_ins[fixed_ref_ins$name %in% coex_te_names[[.x]], ]})

# ----------------------------------------------------------------------------------------------------------
# generate segmentations used for controlling randomization further
# ---------------------------------------------------------------------------------------------------------
set.seed(5)
L_s <- 1e6

te_density_seg <- segmentDensity(ref_ins, n = 2, L_s = L_s, type = "hmm")

g_segments <- plotSegment(te_density_seg[seqnames(te_density_seg) %in% c("2L","2R","3L","3R")])  /
  plotSegment(te_density_seg[seqnames(te_density_seg) %in% c("X","Y","4")]) + plot_layout(guides = "collect")

segmentation_export <- list(g_segments = g_segments,gr_te_density_seg = te_density_seg)

# find which state name corresponds to high TE density, ie hetchrom
high_density_state <- split(te_density_seg,te_density_seg$state) |>
  map(~{.x$counts}) |>
  map(~{mean(.x,na.rm=T)}) %>%
  {names(.)[which.max(.)]}


# ------------------------------------------------------------------------------
# construct results table
# ------------------------------------------------------------------------------
# every combo of ChIP, te set (specific or all fixed TEs), and masking
# skeleton of results
df <- expand_grid(ChIP = names(replicated_pks),
            masking = c("none","mask_heterochromatin"),
            te_set = c("factor_specific","all")) |>
  filter(!((ChIP=="H3K9Me3") & te_set == "factor_specific")) |>
  mutate(mask_gr = map(masking,~{
    if (.x == 'none') {
      return(NA)
    } else if (.x == "mask_heterochromatin") {
      return(te_density_seg[te_density_seg$state==high_density_state])
    }
    })) |>
  mutate(te_gr = pmap(list(ts=te_set,f=ChIP,m=mask_gr),
                      function(ts,f,m){ # function gets either all tes, or factor-specific tes, then subsets them depending on masking
                        if (ts=="factor_specific") {
                            gr <- factor_specific_fixed_ref_ins[[f]]
                          } else {
                            gr <-fixed_ref_ins
                          }
                        if (is(m,"GRanges")) {
                          return(subsetByOverlaps(gr, m, invert = T))
                        } else {
                          return(gr)
                        }
                        })) |>
  mutate(pk_set = map2(ChIP,mask_gr, ~{
    gr <- replicated_pks[[.x]]
    if (is(.y,"GRanges")) {
      return(subsetByOverlaps(gr, .y, invert = T))
    } else {
      return(gr)
    }
  }))



df <- df |>
  mutate(regioneR_results = pmap(list(pk=pk_set,te=te_gr,m=mask_gr), 
                                 .f=function(pk,te,m) {
                                   set.seed(20)
                                   permTest(A=pk, 
                                            B=te, 
                                            randomize.function=randomizeRegions,
                                            evaluate.function=numOverlaps,
                                            alternative = "greater", 
                                            mask=m, 
                                            genome = g_r, # this comes from the top of the script, always the same for each test
                                            ntimes = 1000,
                                            mc.set.seed=FALSE)
                                   }))

# function for nicely plotting randomization results
plot_regioner <- function(x) {
  permuted_overlaps <- tibble(overlaps = x$numOverlaps$permuted)
  observed_overlaps <- tibble(overlaps = x$numOverlaps$observed)
  pval <- x$numOverlaps$pval
  z <- x$numOverlaps$zscore
  alt <- x$numOverlaps$alternative
  
  ggplot(permuted_overlaps, aes(overlaps)) +
    geom_histogram() +
    geom_vline(data = observed_overlaps, aes(xintercept = overlaps), color="darkgreen") +
    annotate("text", x=-Inf, y=Inf, label=sprintf("z=%s, p=%s, alternative=%s", format(z, digits=3), format.pval(pval,digits = 2), alt), hjust=0, vjust=1)
}

extract_regioner_results <- function(x){
  tibble(pval=x$numOverlaps$pval,
         z = x$numOverlaps$zscore,
         hypothesis = x$numOverlaps$alternative)
}

# make this easy to keep together with the dataset and use downstream
# save the plotting function with each row
df <- df |> mutate(plot_results_function = map(ChIP,~return(plot_regioner)))

df <- df |> mutate(regioner_tbl = map(regioneR_results, extract_regioner_results)) |>
  unnest(regioner_tbl)

write_rds(df,snakemake@output$rds)
write_rds(segmentation_export, snakemake@output$segmentation)