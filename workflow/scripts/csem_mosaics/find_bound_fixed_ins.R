Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(plyranges)
library(nullranges)
library(plyranges)
library(tidyverse)
library(AnnotationHub)
library(patchwork)
library(Rsamtools)
library(GenomicAlignments)

# ------------------------------------------------------------------------------
# get pan pks
# ------------------------------------------------------------------------------

pks <- "upstream/pan_ENCSR058DSI_rep*.mosaics.bed" |>
  Sys.glob() |>
  map(read_bed) |>
  map(reduce_ranges) |>
  GRangesList() |>
  unlist() |>
  reduce_ranges(n=plyranges::n()) |>
  filter(n>=2)


# ------------------------------------------------------------------------------
# get fixed, ref ins of TEs coex'd with pan
# ------------------------------------------------------------------------------
# consider splitting
g <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6
seqlevelsStyle(g) <- "NCBI"

coex_te_names <- jsonlite::read_json("results/resources/coexpressed_tes.json") |>
  pluck('all') |>
  pluck('pan') |>
  unlist()

# get reference insertions
ref_ins <- read_bed("upstream/reference_insertions.bed") |>
  filter(!str_detect(name,"\\(.+\\)n")) |>
  sort()

# populate sequence lenghts field, required for some functions in regioner
csizes <- seqinfo(g) |> as.data.frame() |> DataFrame() |> 
  as_tibble(rownames="seqnames") |>
  dplyr::select(seqnames, seqlengths) |>
  deframe()

seqlengths(ref_ins) <- csizes[seqlevels(ref_ins)]

# get penetrance from tidal data
tidal_ins <- read_rds("results/resources/annotated_fixed_insertions.gr.rds")[,1:12]

# from the tidal insetions, subset to only retain plausibly fixed. in this case this means
# present in at least N DGRP lines and the reference
fixed_ref_ins <- tidal_ins |> 
  filter(score > 100) |>
  filter_by_overlaps(ref_ins)
  
# only the reference insertions that are also families coexpressed with pan
pan_ref_ins <- fixed_ref_ins[fixed_ref_ins$name %in% coex_te_names, ]

# for tracking left and right sides of each insertion later
pan_ref_ins <- mutate(pan_ref_ins, id= paste(as.character(pan_ref_ins), pan_ref_ins$name, sep="__"))

# filter to exclude pan ref insertions that abut another TE too closely to  convince us that it is actually bound
# turns out this doesn't much matter, but kept in to be safe
non_fixed_ref_ins <- filter_by_non_overlaps(ref_ins, fixed_ref_ins)
pan_ref_ins <- add_nearest_distance(pan_ref_ins, non_fixed_ref_ins)

pan_ref_ins <- pan_ref_ins |> subsetByOverlaps(pks)

# ------------------------------------------------------------------------------
# get read counts overlapping each edge of the fixed pan coex tes
# ------------------------------------------------------------------------------

# get each edge of the insertions
lefts <- plyranges::flank_left(pan_ref_ins, width = 1)
rights <- plyranges::flank_right(pan_ref_ins, width = 1)

# get couts
lefts_rep1 <- readGAlignments("upstream/pan_ENCSR058DSI_rep1.sorted.bam",param = ScanBamParam(which=lefts,mapqFilter = 3),with.which_label = T)
rights_rep1 <- readGAlignments("upstream/pan_ENCSR058DSI_rep1.sorted.bam",param = ScanBamParam(which=rights, mapqFilter = 3),with.which_label = T)
lefts_rep2 <- readGAlignments("upstream/pan_ENCSR058DSI_rep2.sorted.bam",param = ScanBamParam(which=lefts,mapqFilter = 3),with.which_label = T)
rights_rep2 <- readGAlignments("upstream/pan_ENCSR058DSI_rep2.sorted.bam",param = ScanBamParam(which=rights, mapqFilter = 3),with.which_label = T)

# gets list of high penetrance insertions annotated by the number of reads overlapping the left and right flanks
# then filters for at least 5
candidates <- list(lefts_rep1=lefts_rep1, rights_rep1=rights_rep1, lefts_rep2=lefts_rep2, rights_rep2=rights_rep2) |>
  map(.f = \(x) GRanges(x) |> join_overlap_inner(pan_ref_ins)) |>
  map_df(.f=\(x) x |> as_tibble() |> dplyr::count(insertion = id),.id="side") |>
  pivot_wider(names_from = "side", values_from = "n",values_fill = 0) |>
  filter(if_all(contains("_rep"),~{.x > 3})) |>
  separate(insertion, into=c("chr","start","end","te"), sep="[:-]|__") |>
  mutate(across(c(start, end),as.integer))

# ------------------------------------------------------------------------------
# annotate these by closest gene
# ------------------------------------------------------------------------------

txdb <- AnnotationDbi::loadDb("results/resources/txdb")
seqlevelsStyle(txdb) <- "NCBI"

txs <- transcriptsBy(txdb) |>
  unlist() |>
  anchor_5p() |>
  mutate(width=1)
  
txs$gene_id <- names(txs)

cd2 <- GRanges(candidates) |>
  join_nearest(txs,distance = "dist")

cd2 <- subsetByOverlaps(cd2,pks)

#Biostrings::getSeq(g,cd2) %>% 
#  {names(.) <- as.character(cd2); .} |>
#  writeXStringSet("~/work/231107_csem_semi-fixed-ref-bound-tes/bound.fasta")

#x <- read_rds("results/deg/ourKD.de.grs.rds")
#x2 <- x$adjusted$knockdown2_pan_female_head_Mef2.R_control_female_head_Mef2.R
#x2$gene_id <- names(x2)
#x2[unique(cd2$gene_id[cd2$gene_id %in% names(x2)]),] |>
#  filter(padj < 0.1)


write_rds(cd2, "resources/putatively_bound_insertions.rds")


#write_bed(cd2,"~/Downloads/cd2.bed")
#write_bed(pks,"~/Downloads/pks.bed")
#write_bed(rights,"~/Downloads/rights.bed")
