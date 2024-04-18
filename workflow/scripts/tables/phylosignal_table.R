library(tidyverse)
library(writexl)

x <- list(`coexpression score phyloSignal results`=read_rds("results/phylosignal/phylosignal_df.rds"),
     `Unr RIP-seq phyloSignal`=read_rds("results/ripseq/unr_ripseq_phylosignal.tbl.rds"))


write_xlsx(x,snakemake@output$xlsx)