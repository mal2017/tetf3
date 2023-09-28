# ----------------------------------------------------------------------------
#
#
# not used
#
# 
# ----------------------------------------------------------------------------

library(tidyverse)

mash_dist_tbl <- ifelse(exists("snakemake"),snakemake@input[["mash_dist"]],"results/te_sequence_similarity/te_mash.txt") %>%
  read_tsv(col_names = c("seqA","seqB","mash_dist","pval","shared_hashes"))


kmer.dist <- mash_dist_tbl %>%
  dplyr::select(seqA,seqB,mash_dist) %>%
  pivot_wider(names_from = seqB,values_from = mash_dist,values_fill = 1) %>%
  column_to_rownames("seqA") %>%
  as.dist()

saveRDS(kmer.dist, snakemake@output[["rds"]])