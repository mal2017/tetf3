Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(writexl)

fully_filtered <- read_tsv("results/phylosignal/phylosignal_filtered_hits.tsv.gz") |>
  filter(score_type == "score" & n_tests_sig >=1) |>
  dplyr::select(TF) |>
  distinct()

te_regs <- read_tsv("results/resources/pirna_pathway.tsv")

te_regs_intersection <- fully_filtered |>
  filter(TF %in%te_regs$gene_symbol)

x <- list(`coexpression score phyloSignal results`=read_rds("results/phylosignal/phylosignal_df.rds") |> filter(metric!="Lambda"),
     `Unr RIP-seq phyloSignal`=read_rds("results/ripseq/unr_ripseq_phylosignal.tbl.rds"),
     `>=1 significant tests M or F`=fully_filtered,
     `intersection with piRNA regulators`=te_regs_intersection)


write_xlsx(x,snakemake@output$xlsx)