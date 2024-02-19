library(tidyverse)
library(clusterProfiler)

lkup_fl <- "~/work/tetf3/results/resources/gene_symbol_lookup.tsv.gz"
lkup_fl <- snakemake@input$lkup
lkup <- read_tsv(lkup_fl)

rip_fl <- "results/ripseq/unr_ripseq.tsv.gz"
rip_fl <- snakemake@input$rip
rip <- read_tsv(rip_fl)

gs <- rip |> 
  unite("gs",status,type) |>
  filter(gs!="not_bound_gene") |> # this would be an absurdly large gs
  left_join(lkup, by=c(feature="gene_symbol")) |>
  group_by(gs) |>
  #slice_max(rip_stat,n=100,with_ties = F) |>
  mutate(gene_ID = if_else(is.na(gene_ID),feature,gene_ID)) |>
  dplyr::select(gs,gene_ID) |>
  arrange(gs)

kd_fl <- "results/deg/ourKD.de.df.rds"
kd_fl <- snakemake@input$kd
kd <- read_rds(kd_fl) |>
  filter(group == "head" & comparison == "knockdown2_Unr_female_head_Mef2.R")

rnk <- kd |>
  dplyr::select(feature,stat) |>
  arrange(-stat) |>
  distinct() |>
  deframe() |>
  sort(decreasing = T)

gsea_res <- GSEA(rnk,TERM2GENE = gs,pvalueCutoff = 1,maxGSSize = Inf,eps=0)

#gsea_res |>
#  as_tibble()

#enrichplot::gseaplot2(gsea_res,geneSetID = "bound_gene")

write_rds(gsea_res,snakemake@output$rds)