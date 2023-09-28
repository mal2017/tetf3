library(tidyverse)

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]

mods <- vroom::vroom(mods_path)

res <- mods |>
  mutate(var.explained.by.gex = sumsq_anova_x/total_variance) |>
  filter(significant_model) |>
  dplyr::select(model, gene_symbol, feature.y, estimate, estimate.qnorm, var.explained.by.gex, p.value_anova_x)

write_rds(res, snakemake@output[["rds"]])

#ggplot(aes(p.value_anova_x,fill=model)) +
#  geom_histogram() +
#  facet_wrap(~model)