library(Biostrings)
library(tidyverse)

relpos <- 50
relpos <- snakemake@params$relpos
relpos <- relpos/100

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


te_fa <- "upstream/Tidalbase_transposon_sequence.dmel.fixed.sensechecked.fa.gz"
te_fa <- snakemake@input$te_fa
te_seqs <- Biostrings::readDNAStringSet(te_fa)

combined_seqs <- c(te_seqs,representative_mRNA_seqs)


# get the region of the tx we found enriched for AT content
combined_seqs <- names(combined_seqs) |>
  set_names() |>
  map(~{combined_seqs[[.x]]}) |>
  map(~{
    return(.x[round(length(.x)*relpos):length(.x)])
  }) |>
  DNAStringSet()


z <- combined_seqs |>
  as.list() |>
  map(function(.x) {
    countPattern("ATTTA",.x,fixed = F)
  })

combined_lens <- width(combined_seqs) |> set_names(names(combined_seqs)) |> enframe("feature","seqlen")

z2 <- z |> 
  as_vector() |> 
  enframe(name="feature", value="sites") |>
  left_join(combined_lens) |>
  mutate(site_density = sites/(seqlen/100))

write_tsv(z2,snakemake@output$tsv)