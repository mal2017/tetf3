Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(Biostrings)
library(tidyverse)

fa_fl <- "upstream/Tidalbase_transposon_sequence.dmel.fixed.sensechecked.fa.gz"
fa_fl <- snakemake@input$fa

fa <- readDNAStringSet(fa_fl)

res_fl <- "results/deg/ourKD.de.df.rds"
res_fl <- snakemake@input$res
res <- read_rds(res_fl)

sig_dn <- res |>
  filter(comparison == "knockdown2_CG16779_male_gonad_aTub") |>
  filter(padj < 0.1 & !str_detect(feature,"FBgn")) |>
  filter(log2FoldChange < 0) |>
  pull(feature)


ns <- res |>
  filter(comparison == "knockdown2_CG16779_male_gonad_aTub") |>
  filter(padj >= 0.1 & !str_detect(feature,"FBgn")) |>
  pull(feature)

stopifnot(length(intersect(sig_dn,ns))==0)

writeXStringSet(fa[sig_dn],snakemake@output$fa)
writeXStringSet(fa[ns],snakemake@output$nfa)
