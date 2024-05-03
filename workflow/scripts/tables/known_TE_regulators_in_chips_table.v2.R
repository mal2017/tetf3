Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(GenomicRanges)
library(gt)
library(plyranges)
library(writexl)

te_regulators <- read_tsv("results/resources/pirna_pathway.tsv")
lkup <-read_tsv("results/resources/gene_symbol_lookup.tsv.gz")

d <- read_rds("results/pirna/encode_peaks_dist_to_pirna.gr.rds") |>
  filter(ChIP !="gro")

oi <- c("arx","piwi","Panx","nxf2","mael",
        "fs(1)Yb",
        "shu","SoYb","vret","armi","mino","zuc","Gasz","daed","papi",
        "cuff","moon","del","Boot","rhi","Nxf3",
        "vas","qin","squ","tapas","BoYb","krimp","spn-E","tej","AGO3","aub",
        "egg","mael","nxf2","Nxt1","arx","Su(var)205","piwi") |> unique()

oi <- filter(lkup, gene_symbol %in% oi) |> dplyr::select(gene_id=gene_ID,gene_symbol)

dx <- d |>
  as_tibble() |>
  filter(embryo.expressed) |>
  inner_join(oi) |>
  dplyr::select(feature="gene_id",gene_symbol,ChIP,distance) |>
  distinct() |>
  filter(distance < 5000) |>
  pivot_wider(names_from = ChIP, values_from = distance) |>
  dplyr::rename(`gene ID`=feature, `symbol`=gene_symbol) |>
  dplyr::select(-`gene ID`)

write_xlsx(list(`Closest ChIP peak to piRNA pathway promoters`=dx),snakemake@output$xlsx)