library(tidyverse)
library(clusterProfiler)
library(writexl)

#gsea_tbl_path <- "results/signatures/ourKD_gsea.rds"
gsea_tbl_path <- snakemake@input[["gsea_tbl"]]

gsea_tbl <- read_rds(gsea_tbl_path) 

res <- gsea_tbl %>%
  dplyr::select(-c(signature,gsea,ID,Description)) %>%
  mutate(signature_name = if_else(kd == signature_name,"factor-specific coexpressed TEs",signature_name))

split(res,res$signature_name) |>
  write_xlsx(snakemake@output$xlsx)
