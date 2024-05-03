Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(org.Dm.eg.db)
library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "dmelanogaster_gene_ensembl")

ids <- getBM(filters = "interpro",
    values = "IPR012934",
    attributes = c("ensembl_gene_id"),
    mart = ensembl)

ids <- ids %>% 
  as_tibble %>%
  distinct() %>%
  mutate(gene_symbol = mapIds(org.Dm.eg.db, 
                              keys=ensembl_gene_id, 
                              column="SYMBOL", 
                              keytype="ENSEMBL",
                              multiVals="first")) %>%
  distinct()

write_tsv(ids,snakemake@output$tsv)
