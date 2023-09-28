library(tidyverse)
library(ggdensity)

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]

mods <- read_tsv(mods_path)

mods <- mods %>% filter(significant_x)

probs <- c(0.3,0.6,0.9)

genes_df <- mods %>%
  count(model,feature=feature.x,relationship=if_else(sign(estimate.qnorm)>0,"pos","neg")) %>%
  pivot_wider(names_from = relationship, values_from = n,values_fill = 0)

te_df <- mods %>%
  count(model,feature=feature.y, relationship=if_else(sign(estimate.qnorm)>0,"pos","neg")) %>%
  pivot_wider(names_from = relationship, values_from = n,values_fill = 0)

res <- bind_rows(gene = genes_df, TE = te_df, .id = "feature.type")

                 
write_rds(res,snakemake@output[["rds"]])
