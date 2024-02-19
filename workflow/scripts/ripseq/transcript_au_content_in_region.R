library(rtracklayer)
library(jsonlite)
library(plyranges)
library(tidyverse)
library(Biostrings)


relpos <- 50
relpos <- snakemake@params$relpos
relpos <- relpos/100

# see my tetf_refs pipeline
te_fa <- "upstream/Tidalbase_transposon_sequence.dmel.fixed.sensechecked.fa.gz"
te_fa <- snakemake@input$te_fa
te_seqs <- Biostrings::readDNAStringSet(te_fa)

# get the longest tx from each gene as representative
tx_fa <- "resources/dmel-all-transcript-r6.41.fasta.gz"
tx_fa <- snakemake@input$tx_fa
tx_seqs <- rtracklayer::import(tx_fa)

representative_mRNA_seqs <- split(tx_seqs,str_extract(names(tx_seqs),"(?<=name=).+(?=-R)")) |>
  as.list() |>
  map(function(x){
    return(x[[which.max(width(x))]])
  }) |>
  DNAStringSet()

# the missing ones are nc
# filter(rip,!feature %in% names(representative_mRNA_seqs))
combined_seqs <- c(representative_mRNA_seqs, te_seqs)

# get the region of the tx we found enriched for AT content
combined_seqs <- names(combined_seqs) |>
  set_names() |>
  map(~{combined_seqs[[.x]]}) |>
  map(~{
    return(.x[round(length(.x)*relpos):length(.x)])
  }) |>
  DNAStringSet()

# nt freq for each sequence
x <- combined_seqs |> letterFrequency(letters = "AT",as.prob = T) |>
  as_tibble() |>
  mutate(feature=names(combined_seqs)) |>
  dplyr::select(feature,nt.content=`A|T`)

write_tsv(x,snakemake@output$tsv)
