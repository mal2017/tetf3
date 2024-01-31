library(jsonlite)
library(tidyverse)

#mod_fl <- "upstream/final-models.collected-info.tsv.gz"
mod_fl <- snakemake@input$mods

# get signif hits for each sex, with gene symbol and id included
coex_te_names <- read_tsv(mod_fl) |> 
  filter(significant_x) |> 
  dplyr::select(gene_id = feature.x, gene_symbol, te=feature.y, model) |>
  distinct()

# add the deduplicated set of tes from both sexes for each gene
coex_te_names <- coex_te_names |>
  mutate(model = "all") |>
  distinct() |>
  bind_rows(coex_te_names)

# allows lookup by sex (or both sexes) and by either gene id or gene symbol
coex_te_names <- coex_te_names |>
  pivot_longer(c(gene_symbol, gene_id),names_to = "feature_type",values_to = "feature") |>
  dplyr::select(-feature_type) |>
  relocate(feature, .before="te")

nested_list <- coex_te_names |> 
  group_by(feature,model) |>
  summarize(te=list(te)) |>
  ungroup() |>
  nest(data=-model) |>
  deframe() |>
  map(deframe)

jsonlite::write_json(nested_list,snakemake@output$json,pretty=T,simplifyVector=F)
