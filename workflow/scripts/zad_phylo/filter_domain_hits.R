library(Biostrings)
library(tidyverse)

# this assumes the sequence header is in the format
# returned when running hmmsearch > esl-reformat

n <- ifelse(exists("snakemake"),snakemake@params$n,1) %>% as.integer()

message(n)

fa <- ifelse(exists("snakemake"),snakemake@input$fasta,"results/hmmer/domain-hits/PF07776.fasta")
fa <- readAAStringSet(fa)

to_keep <- tibble(name = names(fa)) %>%
  mutate(starts = str_extract(name, "(?<=\\/)\\d+(?=-)") %>% as.integer()) %>%
  mutate(ix = row_number()) %>%
  mutate(gene_id = str_extract(name,"^.+?(?=\\.FBpp)")) %>%
  group_by(gene_id) %>%
  slice_min(ix,n=n, with_ties = F)

fa <- fa[to_keep$ix]  

#names(fa) <- str_extract(names(fa),"^.+?(?=\\.FBpp)") %>% make.unique()

writeXStringSet(fa, snakemake@output$fasta)