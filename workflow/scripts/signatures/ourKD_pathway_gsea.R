library(clusterProfiler)
library(tidyverse)
library(GenomicRanges)
library(tidyverse)
library(clusterProfiler)
library(org.Dm.eg.db)


gene_universe <- read_tsv("resources/fbgn_fbtr_fbpp_expanded_fb_2021_04.tsv.gz", skip = 4)

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()

# ------------------------------------------------------------------------------
# pathways to test
# ------------------------------------------------------------------------------

pirna_tbl <- ifelse(exists("snakemake"),snakemake@input$pirna,"results/resources/pirna_pathway.tsv") %>% read_tsv()

pirna_tbl <- pirna_tbl %>% pivot_longer(cols = contains("in."),names_to = "gs_name", values_to = "is") %>%
  mutate(gs_name = 'piRNA') %>%
  dplyr::select(gs_name,ensembl_gene = gene_ID)  |>
  distinct()

go_tbl <- "https://ftp.flybase.net/releases/FB2023_06/precomputed_files/go/gene_association.fb.gz"

go_tbl <- readr::read_tsv(go_tbl,skip=5,col_names = NA)
go_tbl2 <- go_tbl |> dplyr::select(gs_name=X5,ensembl_gene=X2) |>
  distinct() |>
  arrange(gs_name)

# ------------------------------------------------------------------------------
# get rankings
# ------------------------------------------------------------------------------

res_path <- "results/deg/ourKD.de.grs.rds"
res_path <- snakemake@input[["res"]]

res <- read_rds(res_path) %>%
  map_df(~map_df(.x,as_tibble,.id="comparison"),.id="adjustment")

# add gene symbol name for each feature
res <- res %>% 
  filter(adjustment %in% c('adjusted')) %>%
  left_join(lkup, by=c(feature="gene_ID")) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol),feature, gene_symbol)) %>%
  dplyr::select(comparison,adjustment,score=log2FoldChange, feature, gene_symbol)

res <- mutate(res, kd = str_extract(comparison,"(?<=knockdown2_).+?(?=_)"))

# make df with a row for every comparison we want
# per kd, one row for the knocked-down gene's TE signature
# plus one for all TEs
res <- res %>%
  dplyr::select(comparison,adjustment,kd,feature,score) %>%
  arrange(-score) %>% 
  nest(data=c(feature,score))

# ------------------------------------------------------------------------------
# run test
# ------------------------------------------------------------------------------

# run gsea
possibly_gsea <- possibly(function(.x){GSEA(deframe(.x), TERM2GENE = go_tbl2,seed=2022,pvalueCutoff = 1,minGSSize = 10, eps=0)},otherwise = NULL)
#possibly_gsea <- possibly(function(.x){gseGO(deframe(.x), seed=2022,OrgDb = org.Dm.eg.db,ont="BP",keyType = "FLYBASE",pvalueCutoff = 1,minGSSize = 10, eps=0)},otherwise = NULL)

set.seed(2)

res <- res |>
  mutate(gsea = map(data,possibly_gsea)) %>%
  mutate(gsea.tidy = map(gsea,as_tibble)) %>%
  dplyr::select(-data)
  #unnest(gsea.tidy)

saveRDS(res, snakemake@output[["rds"]])
