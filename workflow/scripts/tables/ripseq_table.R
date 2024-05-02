library(tidyverse)
library(writexl)

attta <- read_tsv("results/ripseq/unr_ripseq_features_attta_sites_in_region.tsv.gz")

at <- read_tsv("results/ripseq/unr_ripseq_features_au_content_in_region.tsv.gz")

rip <- read_tsv("results/ripseq/unr_ripseq.tsv.gz") |> filter(comparison == "WT_vs_CTRL")

list(`RIP-seq enrichment results`=rip,
     `ATTTA sites in highlighted region`=attta,
     `AT content in highlighted region`=at) |>
  write_xlsx(snakemake@output$xlsx)
