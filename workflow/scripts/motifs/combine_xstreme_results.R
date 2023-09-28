library(magrittr)
library(tidyverse)

#tsvs <- Sys.glob("results/motifs/xstreme_per_tf/pan/xstreme.tsv")
tsvs <- snakemake@input[["memes"]] %>% paste0("/xstreme.tsv")
tsvs <- tsvs %>% set_names(.,str_extract(.,"(?<=tf\\/).+(?=\\/xstreme.tsv)"))
res <- tsvs %>% 
  map_df(read_tsv, comment="#",.id = "te_group")


write_tsv(res, snakemake@output[["tsv"]])
