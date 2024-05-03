Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)
library(rtracklayer)


fa <- rtracklayer::import(ifelse(exists("snakemake"),snakemake@params[["genome_fa"]],"resources/dmel-all-chromosome-r6.41.fasta.gz"))

names(fa) <- names(fa)%>% str_extract(".+(?=\\s+type)")

export(fa,snakemake@output[["genome_fa"]],format="fasta")