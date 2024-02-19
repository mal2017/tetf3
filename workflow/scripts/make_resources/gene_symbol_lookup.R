library(tidyverse)

gene_universe <- read_tsv("resources/fbgn_fbtr_fbpp_expanded_fb_2021_04.tsv.gz", skip = 4)

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol, annotation_ID,gene_type) %>% distinct()

write_tsv(lkup,snakemake@output[["tsv"]])