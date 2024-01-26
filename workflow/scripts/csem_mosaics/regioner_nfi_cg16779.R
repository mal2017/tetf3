library(nullranges)
library(plyranges)
library(tidyverse)
library(regioneR)
library(excluderanges)
library(AnnotationHub)
library(patchwork)

# this is a hefty script, but because this comprises a complicated set of analyses
# with shared dependencies
# I think that's better than trying to split up into components - that could get complicated
# may want to go back and split into smaller tasks for easy use with snakemake

g <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6
seqlevelsStyle(g) <- "NCBI"

g_r <- unlist(tileGenome(seqlengths = seqlengths(g),ntile = 1))

coex_te_names <- read_tsv("upstream/final-models.collected-info.tsv.gz") |> 
  filter(significant_x) |> 
  dplyr::select(gene_id = feature.y, gene_symbol, te=feature.y) |>
  distinct()

# get reference insertions
ref_ins <- read_bed("upstream/reference_insertions.bed") |>
  filter(!str_detect(name,"\\(.+\\)n"))
ref_ins <- sort(ref_ins)

# populate sequence lenghts field, required for some functions in regioner
csizes <- seqinfo(g) |> as.data.frame() |> DataFrame() |> 
  as_tibble(rownames="seqnames") |>
  dplyr::select(seqnames, seqlengths) |>
  deframe()

seqlengths(ref_ins) <- csizes[seqlevels(ref_ins)]

# get penetrance from tidal data
tidal_ins <- rtracklayer::import("~/amarel-matt/tetf/subworkflows/tetf_tidal/results/beds/dgrp_tidal_insertions.unique.bb")

# from the tidal insetions, subset to only retain plausibly fixed (in this case means)
# present at the same penetrance as the most common insertions
fixed_ref_ins <- tidal_ins |> filter(score == max(tidal_ins$score)) |>
  mutate(fixed =T) |>
  join_overlap_left(ref_ins,y=_) |>
  filter(name.x==name.y & fixed) |>
  group_by(name.x) |>
  reduce_ranges(name_list=name.x) |>
  ungroup() |>
  sort() |>
  mutate(name=name.x)

# -----------------------------------------------------------------------------------------------------------
# get pan as the primary peak set, gro as a sort of control
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

nnk_pks <- get_replicated_pks("~/amarel-matt/tetf/subworkflows/tetf_csem_mosaics/results/csem_mosaics/mosaics/CG17802*/*mosaics.bed")

vvl_pks <- get_replicated_pks("~/amarel-matt/tetf/subworkflows/tetf_csem_mosaics/results/csem_mosaics/mosaics/vvl*/*mosaics.bed")

odj_pks <- get_replicated_pks("~/amarel-matt/tetf/subworkflows/tetf_csem_mosaics/results/csem_mosaics/mosaics/CG7357_*_rep*/*rep*.mosaics.bed")

nfi_pks <- get_replicated_pks("~/amarel-matt/tetf/subworkflows/tetf_csem_mosaics/results/csem_mosaics/mosaics/NfI*rep*/*rep*.mosaics.bed")

cg16779_pks <- get_replicated_pks("~/amarel-matt/tetf/subworkflows/tetf_csem_mosaics/results/csem_mosaics/mosaics/CG16779*rep*/*rep*.mosaics.bed")

# ----------------------------------------------------------------------------------------------------------
# generate segmentations used for controlling randomization further on
# ---------------------------------------------------------------------------------------------------------
set.seed(5)
L_s <- 1e6
te_density_seg <- segmentDensity(ref_ins, n = 2, L_s = L_s, type = "hmm")

g_segments <- plotSegment(te_density_seg[seqnames(te_density_seg) %in% c("2L","2R","3L","3R")])  /
  plotSegment(te_density_seg[seqnames(te_density_seg) %in% c("X","Y","4")]) + plot_layout(guides = "collect")


# ------------------------------------------------------------------------------------------------------------
# run regioner and plot
# ----------------------------------------------------------------------------------------------------------

# does factor X overlap tes more than expected by chance
nnk_all_te_overlap.pt <- permTest(A=nnk_pks, B=fixed_ref_ins, randomize.function=randomizeRegions,
                                  evaluate.function=numOverlaps,alternative = "greater", mask=NA, genome = g_r, ntimes = 1000)

vvl_all_te_overlap.pt <- permTest(A=vvl_pks, B=fixed_ref_ins, randomize.function=randomizeRegions,
                                  evaluate.function=numOverlaps,alternative = "greater", mask=NA, genome = g_r, ntimes = 1000)

odj_all_te_overlap.pt <- permTest(A=odj_pks, B=fixed_ref_ins, randomize.function=randomizeRegions,
                                  evaluate.function=numOverlaps,alternative = "greater", mask=NA, genome = g_r, ntimes = 1000)

nfi_all_te_overlap.pt <- permTest(A=nfi_pks, B=fixed_ref_ins, randomize.function=randomizeRegions,
                                  evaluate.function=numOverlaps,alternative = "greater", mask=NA, genome = g_r, ntimes = 1000)

cg16779_all_te_overlap.pt <- permTest(A=cg16779_pks, B=fixed_ref_ins, randomize.function=randomizeRegions,
                                  evaluate.function=numOverlaps,alternative = "greater", mask=NA, genome = g_r, ntimes = 1000)


res <- list(Odj_all_te = odj_all_te_overlap.pt,
            NfI_all_te = nfi_all_te_overlap.pt,
            vvl_all_te = vvl_all_te_overlap.pt,
            Nnk_all_te = nnk_all_te_overlap.pt,
            CG16779_all_te = cg16779_all_te_overlap.pt,
            segmentation = list(gg = g_segments, gr = te_density_seg))

write_rds(res,snakemake@output$rds)