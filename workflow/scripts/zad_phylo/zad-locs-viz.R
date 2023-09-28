library(ggbio)
library(rtracklayer)
library(Biostrings)
library(tidyverse)
library(GenomicFeatures)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(BSgenome.Dmelanogaster.UCSC.dm6)

txdb <- loadDb("~/work/tetf_downstream/results/resources/txdb")
genome <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6
seqlevelsStyle(genome) <- "NCBI"

zad_peps <- readAAStringSet("results/hmmer/domain-hits/PF07776.filtered.fasta")

names(zad_peps) <- names(zad_peps) %>% str_extract(.,"FBgn\\d+(?=\\.)")

zad_coords <- genes(txdb)[names(zad_peps)]
seqlevels(zad_coords) <- seqlevelsInUse(zad_coords)
seqlengths(zad_coords) <- seqlengths(genome)[seqlevels(zad_coords)]

ggbio::autoplot(zad_coords, layout="karyogram",ideogram=T)
