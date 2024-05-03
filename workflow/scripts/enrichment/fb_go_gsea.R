Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)
library(clusterProfiler)
library(org.Dm.eg.db)
library(enrichplot)

ranking <- ifelse(exists("snakemake"),  snakemake@wildcards$ranking, "abs")

rnks <- ifelse(exists("snakemake"),
               snakemake@input$tsv,
               "results/rankings/sig_main_female_max_abs_estimate_qnorm.tsv.gz") %>% 
  read_tsv()

rnks <- rnks %>% dplyr::select(-gene_id) %>% arrange(-value) %>% deframe()

possibly_gsea <- possibly(function(.x,.y) {
  # note that per the function help page, the pv cutoff is for adjusted values
  # also - as written here it expects already unnested inputs
  set.seed(1)
  
  st <- case_when(str_detect(.y,"abs|max|sum")~"pos",
                  str_detect(.y,"min")~"neg",
                  T~"std")
  gseGO(.x, ont = "ALL", OrgDb = org.Dm.eg.db, keyType = "SYMBOL", seed=2022,pvalueCutoff = 0.1, minGSSize = 15,maxGSSize = 300, pAdjustMethod = "BH", scoreType=st, eps=0)
},otherwise = NULL)

res <- possibly_gsea(rnks, ranking)

write_rds(res, snakemake@output$rds)

