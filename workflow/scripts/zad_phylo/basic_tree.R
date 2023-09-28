library(tidytree)
library(treeio)
library(ggtree)
library(tidyverse)


lkup <- read_tsv("http://ftp.flybase.net/releases/FB2021_04/precomputed_files/genes/fbgn_annotation_ID.tsv.gz",skip = 4) %>%
  dplyr::select(gene_symbol = `##gene_symbol`,gene_id = `primary_FBgn#`) %>%
  distinct()

tree <- read.iqtree("results/iqtree/domain-hits/PF07776/PF07776.contree")

tree <- tree %>%
  mutate(gene_id = str_extract(label,"FBgn.+(?=\\.)"))
  
tree <- tidytree::left_join(tree,lkup, by="gene_id")

tree %>%
  ggtree(ladderize = T) +
  geom_nodelab(aes(label=UFboot),size=3,color="blue",hjust=0) +
  geom_tiplab(aes(label=gene_symbol))


tree %>% as_tibble() %>% filter(gene_id == "FBgn0039740")
