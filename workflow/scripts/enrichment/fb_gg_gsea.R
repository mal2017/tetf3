Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)
library(clusterProfiler)


ranking <- ifelse(exists("snakemake"),  snakemake@wildcards$ranking, "abs")

rnks <- ifelse(exists("snakemake"),
               snakemake@input$tsv,
               "results/rankings/nofilt_main_female_mean_abs_estimate_qnorm.tsv.gz") %>% 
  read_tsv()

rnks <- rnks %>% dplyr::select(-gene_id) %>% arrange(-value) %>% deframe()

zad_tbl <- ifelse(exists("snakemake"),snakemake@input$zad,"results/resources/zad_genes.tsv") %>% read_tsv()

zad_tbl <- zad_tbl %>% 
  mutate(gs_name="ZAD_ZNF") %>%
    dplyr::select(gs_name,ensembl_gene = gene_symbol)

pirna_tbl <- ifelse(exists("snakemake"),snakemake@input$pirna,"results/resources/pirna_pathway.tsv") %>% read_tsv()

pirna_tbl <- pirna_tbl %>% pivot_longer(cols = contains("in."),names_to = "gs_name", values_to = "is") %>%
  filter(is) %>%
  mutate(gs_name = str_remove(gs_name,"in.")) %>%
  dplyr::select(gs_name,ensembl_gene = gene_symbol)

animaltfdb <- read_tsv(ifelse(exists("snakemake"),snakemake@input$animaltfdb_tfs,"resources/Drosophila_melanogaster_TF.txt")) %>%
  mutate(gs_name = "AnimalTFDB") %>%
  dplyr::select(gs_name,ensembl_gene = Symbol)

grps <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/gene_group_data_fb_2022_04.tsv.gz",skip=7)

t2g <- grps %>%
  dplyr::select(gs_name=FB_group_name,ensembl_gene=Group_member_FB_gene_symbol) %>%
  distinct() %>%
  drop_na() %>%
  group_by(gs_name) %>%
  filter(n()>10) %>%
  ungroup() %>%
  bind_rows(pirna_tbl) %>%
  bind_rows(zad_tbl) %>%
  bind_rows(animaltfdb)

#t2g %>% filter(str_detect(gs_name,"HIGH MOBIL")) %>% pull(gs_name) %>% unique %>% walk(message)

possibly_gsea <- possibly(function(.x,.y) {
  # note that per the function help page, the pv cutoff is for adjusted values
  # also - as written here it expects already unnested inputs
  set.seed(1)
  st <- case_when(str_detect(.y,"abs|max|sum")~"pos",
                  str_detect(.y,"min")~"neg",
                  T~"std")
  GSEA(.x, TERM2GENE = t2g,seed=2022,pvalueCutoff = 0.1, minGSSize = 15,pAdjustMethod = "BH", scoreType=st, eps=0)
},otherwise = NULL)

res <- possibly_gsea(rnks, ranking)

write_rds(res, snakemake@output$rds)



