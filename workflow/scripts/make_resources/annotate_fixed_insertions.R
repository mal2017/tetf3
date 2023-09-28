library(tidyverse)
library(plyranges)
library(rtracklayer)

lms <- ifelse(exists("snakemake"),snakemake@input[["lms"]],
              "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv() %>%
  filter(significant_x) %>%
  dplyr::rename(sex=model)

# -------- get TFs -------------------------------------------------------------

remap_simple <- ifelse(exists("snakemake"),snakemake@input[["remap"]],
                     "resources/remap2022_nr_macs2_dm6_v1_0.bed.gz") %>%
  import()

remap_simple$ChIP <- remap_simple$name %>% str_extract(".+(?=:)")

# join overlapping ranges within the same TF
remap_simple <- remap_simple %>% 
  split(.,.$ChIP) %>%
  GenomicRanges::reduce()

remap_simple <- remap_simple[names(remap_simple) %in% lms$gene_symbol]

seqlevelsStyle(remap_simple) <- "NCBI"


# -------- get TEs -------------------------------------------------------------
# Which set of TE insertions to use is important to consider.. 

#insertions_path <- "upstream/dgrp_tidal_insertions.unique.bb"
insertions_path <- snakemake@input[["insertions"]]

#all_ins_path <- "upstream/dgrp_tidal_insertions.bb"
all_ins_path <- snakemake@input[["all_ins"]]

ref_ins <- import(all_ins_path) %>%
  filter(str_detect(name,"reference"))

ref_ins <- ref_ins %>% split(.,str_extract(.$name,".+(?=\\.DGRP)"))

tes_p1 <- import(insertions_path) %>%
  mutate(nOL = score)

strand(tes_p1) <- "*"

n_strains <- max(tes_p1$nOL)

# For INE-1, not all are fixed - most likely just small fragments aren't
fixed_in_dgrp <- tes_p1 %>% 
  filter(nOL == n_strains) %>% 
  split(.,.$name) %>% as.list()

# this only retains insertions that are 'fixed_in_dgrp' (ie nOL > x)
fixed_in_ref <- map(names(fixed_in_dgrp),.f=function(x){ 
  print(x)
  subsetByOverlaps(ref_ins[[x]],fixed_in_dgrp[[x]]) %>%
    mutate(name = x)
})

tes_p2 <-  plyranges::select(unlist(GRangesList(fixed_in_ref)), repeat_element = name) %>%
  split(.,.$repeat_element) %>%
  GenomicRanges::reduce() %>% 
  unlist() %>%
  mutate(.,repeat_element = names(.))

tes_p2 <- tes_p2[seqnames(tes_p2) %in% seqlevels(remap_simple)]

strand(tes_p2) <- "*" 

tes <- subsetByOverlaps(tes_p1,tes_p2) %>% filter(nOL==n_strains)


# annotate by feature ----------------------------------------------------------
library(GenomicFeatures)
txdb <- loadDb(ifelse(exists("snakemake"),snakemake@input[["txdb"]],
                                                "results/resources/txdb"))

ex <-exons(txdb) %>% GenomicRanges::reduce() %>% unstrand()

intr <- intronsByTranscript(txdb) %>% unlist() %>% GenomicRanges::reduce() %>% unstrand()

fiveUTR <- fiveUTRsByTranscript(txdb) %>% unlist() %>% GenomicRanges::reduce() %>% unstrand()

threeUTR <- threeUTRsByTranscript(txdb) %>% unlist() %>% GenomicRanges::reduce() %>% unstrand()

prom2kbup <- GenomicFeatures::promoters(txdb,upstream = 2000,downstream = 0) %>% GenomicRanges::reduce() %>% unstrand()

genes <- genes(txdb)

down5kb <- flank_downstream(genes,width = 5000)

tes <- tes %>% 
  mutate(.,
         exonic = unstrand(.) %over% ex,
         intronic = unstrand(.) %over% intr,
         utr5 = unstrand(.) %over% fiveUTR,
         utr3 = unstrand(.) %over% threeUTR,
         prom = unstrand(.) %over% prom2kbup,
         dn5kb = unstrand(.) %over% down5kb)

tes <- tes %>% 
  plyranges::add_nearest_distance(name = "nearest.tss" ,promoters(txdb,upstream = 0,downstream = 0))

# -------- get GC -------------------------------------------------------------
library(BSgenome.Dmelanogaster.UCSC.dm6)

seqlevelsStyle(Dmelanogaster) <- "NCBI" 

tes$GC <- letterFrequency(getSeq(Dmelanogaster, tes), "GC", as.prob=T)[,1]


# ------------------------------------------------------------------------------
# annotate insertion
tes2 <- tes %>% plyranges::mutate(size = width)

# good spot to annotate with all TF overlaps
names(tes2) <- NULL

chip_overlap_mat <- lapply(remap_simple, FUN=function(x){countOverlaps(tes2,x)>0}) %>% 
  do.call(cbind,.,)

# wtf -  Gal here is mislabeled by Remap22 -  it is actually H3K27me3... 
# I traced it to this study: https://www.biorxiv.org/content/10.1101/127951v1.full.pdf
chip_overlap_mat <- chip_overlap_mat[,colnames(chip_overlap_mat)!="Gal"]

mcols(tes2) <- cbind(mcols(tes2),chip_overlap_mat)

write_rds(tes2,snakemake@output[["rds"]])
write_rds(remap_simple,snakemake@output[["remap"]])

