library(biomaRt)
library(tidyverse)
library("org.Dm.eg.db")
library(Biostrings)
library(rtracklayer)

# "results/peptides/all.fasta.gz"
peps <- readAAStringSet(snakemake@input$fasta)

# split by gene
peps <- split(peps,str_extract(names(peps),"FBgn[0-9]+")) %>%
  as.list()

# select longest per gene
peps <- peps %>%
  map(~{.x[which.max(width(.x))]}) %>%
  AAStringSetList() %>%
  unlist()

names(peps) <- str_extract(names(peps),"^FBgn.+(?= type)")

Biostrings::writeXStringSet(peps, snakemake@output$fasta)