# ----------------------------------------------------------------------------
#
#
# not used
#
# 
# ----------------------------------------------------------------------------

library(Biostrings)
library(tidyverse)

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]
mods <- read_tsv(mods_path)

#te_fa <- "upstream/Tidalbase_transposon_sequence.dmel.fixed.sensechecked.fa"
te_fa <- snakemake@input$te_fasta
te_seqs <- readDNAStringSet(te_fa)

te_seqs <- te_seqs[names(te_seqs) %in% mods$feature.y]

pid_df <- names(te_seqs) |>
  tibble(x=_) |>
  crossing(y=x) |>
  filter(x!=y)


pid_df2 <- left_join(pid_df, tibble(te = names(te_seqs),seq = as.character(te_seqs)),by=c(x="te")) |>
  left_join(tibble(te = names(te_seqs),seq = as.character(te_seqs)),by=c(y="te")) |>
  mutate(aln = map2(seq.x,seq.y, .f=function(x,y){pairwiseAlignment(x,y,type="global")})) |>
  mutate(pid1 = map_dbl(aln, pid, type="PID1"),
         pid2 = map_dbl(aln, pid, type="PID2"),
         pid3 = map_dbl(aln, pid, type="PID3"),
         pid4 = map_dbl(aln, pid, type="PID4"))

write_rds(dplyr::select(pid_df2,-seq.x,-seq.y, -aln), snakemake@output$rds)
