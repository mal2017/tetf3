library(tidyverse)
library(memes)

options(meme_bin = "/home/mlawlor/local/meme/bin/")

# get the representative motifs
sim <- read_rds("results/motifs/comparison/pan_denovo_comparison.meme.rds")
to_use <- sim |> filter(padj < 0.1 & (name == "pan" | str_detect(known,"Archbold"))) |> pull(denovo) |>
  str_remove("denovo::") |>
  unique()

meme_res <- memes::importMeme("results/motifs/meme_per_tf/pan/meme.txt") |>
  filter(name %in% to_use)

meme_res <- universalmotif::to_list(meme_res)

cons <- Biostrings::readDNAStringSet("upstream/Tidalbase_transposon_sequence.dmel.fixed.sensechecked.fa.gz")

memes::runFimo(sequences = cons,
               motifs = meme_res, text=F,
               outdir = "results/motifs/fimo_on_tes/denovo/pan/", 
               thresh=0.1, 
               qv_thresh=T)

x <- importFimo("results/motifs/fimo_on_tes/denovo/pan/fimo.tsv")


x$motif_id
