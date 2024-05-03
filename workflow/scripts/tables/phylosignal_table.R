Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(writexl)

fully_filtered <- read_tsv("results/phylosignal/phylosignal_filtered_hits.tsv.gz") |>
  filter(score_type == "score" & n_tests_sig >=3) |>
  dplyr::select(TF) |>
  distinct()

x <- list(`coexpression score phyloSignal results`=read_rds("results/phylosignal/phylosignal_df.rds"),
     `Unr RIP-seq phyloSignal`=read_rds("results/ripseq/unr_ripseq_phylosignal.tbl.rds"),
     `>=3 significant tests M or F`=fully_filtered)


write_xlsx(x,snakemake@output$xlsx)