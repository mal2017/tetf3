library(universalmotif)
library(Biostrings)

#seed = as.integer(runif(n=1, min=1, max=1000))
seed = as.numeric(snakemake@params$seed)

e <- Biostrings::readDNAStringSet(snakemake@params$fasta)

shuf <- shuffle_sequences(e, k = 2, method = "markov", window = T, rng.seed = seed)

Biostrings::writeXStringSet(shuf, snakemake@output$p)
