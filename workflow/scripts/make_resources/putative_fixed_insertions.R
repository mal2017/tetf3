library(plyranges)
library(tidyverse)

# get reference insertions
ref_ins_fl <- "upstream/reference_insertions.bed"
ref_ins_fl <- snakemake@input$ref_ins
ref_ins <- read_bed(ref_ins_fl) |>
  filter(!str_detect(name,"\\(.+\\)n"))
ref_ins <- sort(ref_ins)

# get penetrance from tidal data
#tidal_ins_fl <- "upstream/dgrp_tidal_insertions.unique.bb"
tidal_ins_fl <- snakemake@input$uniq_tidal
tidal_ins <- rtracklayer::import(tidal_ins_fl)

# from the tidal insertions, subset to only retain plausibly fixed
# in this case means fully in all DGRP strains and the reference
fixed_ref_ins <- tidal_ins |> filter(score == max(tidal_ins$score)) |>
  mutate(fixed =T) |>
  join_overlap_left(ref_ins,y=_) |>
  filter(name.x==name.y & fixed) |>
  group_by(name.x) |>
  reduce_ranges(name_list=name.x) |>
  ungroup() |>
  sort() |>
  mutate(name=name.x)


write_rds(fixed_ref_ins, snakemake@output$rds)