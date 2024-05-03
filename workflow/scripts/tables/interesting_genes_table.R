Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(writexl)

mods <- read_tsv("upstream/final-models.collected-info.tsv.gz")


# genes with more than 1 coex te
n_coex <- mods |> 
  filter(significant_x) |>
  dplyr::select(gene_symbol,feature.y) |>
  distinct() |>
  group_by(gene_symbol) |>
  summarise(n=n(),.groups = "drop")

more_than_1_coex <- n_coex |> filter(n > 1)

# TEs with number of coexpressed TEs > pan
more_coex_tes_than_pan <- n_coex |> filter(n>pull(filter(n_coex,gene_symbol == "pan"),n))

# genes with median abs score > pan & >=3 coex tes
med_abs_scores <-mods |> 
  filter(significant_x) |>
  group_by(model,feature.x,gene_symbol) |>
  summarise(across(c(estimate.qnorm),.fns=list(median_abs=~median(abs(.x))),.names = "{.fn}_{col}"),n=n()) |>
  ungroup() |>
  filter(n >= 3)

pan_male_med_abs_coex <- pull(filter(med_abs_scores,gene_symbol == "pan" & model=="male"),median_abs_estimate.qnorm)
pan_female_med_abs_coex <- pull(filter(med_abs_scores,gene_symbol == "pan" & model=="female"),median_abs_estimate.qnorm)

male_med_abs_greater_than_pan <- med_abs_scores |> filter(model == "male" & median_abs_estimate.qnorm > pan_male_med_abs_coex)
female_med_abs_greater_than_pan <- med_abs_scores |> filter(model == "female" & median_abs_estimate.qnorm > pan_male_med_abs_coex)


res <- list(`genes with more than 1 coexpressed TE`=more_than_1_coex,
     `genes with more coexpressed TEs than pan`=more_coex_tes_than_pan,
     `genes with higher male median abs. coex score than pan`=male_med_abs_greater_than_pan,
     `genes with higher female median abs. coex score than pan`=female_med_abs_greater_than_pan)

write_xlsx(res,snakemake@output$xlsx)