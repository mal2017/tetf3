library(tidyverse)
library(clusterProfiler)
library(writexl)

f_path <- "results/enrichment/sig_main_female_max_abs_estimate_qnorm.gg_gsea.rds"
m_path <- "results/enrichment/sig_main_male_max_abs_estimate_qnorm.gg_gsea.rds"
f_path <- snakemake@input[["gsea_tbl_f"]]
m_path <- snakemake@input[["gsea_tbl_m"]]

list(male=m_path,
     female = f_path) |>
  map(read_rds) |>
  map(as_tibble) |>
  write_xlsx(snakemake@output$xlsx)
