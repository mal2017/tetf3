Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(rtracklayer)
library(jsonlite)
library(plyranges)
library(tidyverse)
library(Biostrings)

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

# get average nt fequency within 100nt sliding window over length of each sequence
x <- combined_seqs |>
  as.list() |>
  map_df(function(.x) {
    letterFrequencyInSlidingView(.x,view.width = 100,c("AT"),as.prob = T) |>
      as_tibble() |>
      dplyr::rename(nt.content='A|T') |>
      mutate(position = row_number())
  },.id="feature")

# rescale position to percent out of 100
x2 <- x |>
  group_by(feature) |>
  mutate(position.scaled=scales::rescale_max(position,to=c(0,100))) |>
  ungroup() |>
  mutate(position.scaled=ceiling(position.scaled)) |>
  group_by(feature,position.scaled) |>
  summarize(nt.content=mean(nt.content),approx.real.position=mean(position),.groups = "drop")

write_tsv(x2,snakemake@output$tsv)
