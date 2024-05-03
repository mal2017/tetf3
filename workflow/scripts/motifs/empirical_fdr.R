Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(magrittr)
library(tidyverse)

#tsvs <- Sys.glob("results/motifs/streme_per_tf/pan/sequences.tsv")
tsvs <- snakemake@input[["streme"]] %>% paste0("/sequences.tsv")
tsvs <- tsvs %>% set_names(.,str_extract(.,"(?<=tf\\/).+(?=\\/sequences.tsv)"))

res <- tsvs %>% 
  map_df(read_tsv, comment="#",.id = "te_group") |>
  dplyr::select(te_group, motif_ID, motif_Score) |>
  distinct()


#shuf <- Sys.glob("results/motifs/streme_per_tf_shuffled/pan/*/sequences.tsv")
shuf <- snakemake@input[["shuf"]] %>% paste0("/sequences.tsv")
shuf <- shuf %>% set_names(.,str_extract(.,"(?<=\\/)\\d+(?=\\/sequences.tsv)"))

res_shuf <- shuf %>% 
  map_df(read_tsv, comment="#",.id = "te_group") |>
  dplyr::select(te_group, motif_ID, motif_Score) |>
  distinct()

# add empirical fdr
res <- res |> mutate(fdr = map_dbl(motif_Score, ~{mean(res_shuf$motif_Score < .x)}))

dat <- bind_rows(res,res_shuf)

write_tsv(dat, snakemake@output[["tsv"]])
