Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)
library(Biostrings)

te_fa <- "upstream/Tidalbase_transposon_sequence.dmel.fixed.sensechecked.fa.gz"
te_fa <- snakemake@input$te_fa
tes <- Biostrings::readDNAStringSet(te_fa)

rip_fl <- "results/ripseq/unr_ripseq.tsv.gz"
rip_fl <- snakemake@input$rip
rip <- read_tsv(rip_fl)

rip_tes <- rip |>
  filter(type == "TE" & status == "bound") |>
  pull(feature) |>
  unique()

non_rip_tes <- rip |>
  filter(type == "TE" & status != "bound") |>
  pull(feature) |>
  unique()

writeXStringSet(tes[rip_tes],snakemake@output$fa)
writeXStringSet(tes[non_rip_tes],snakemake@output$non_bound_fa)
