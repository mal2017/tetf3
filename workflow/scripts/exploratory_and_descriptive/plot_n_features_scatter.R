Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)
library(ggdensity)

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]

mods <- read_tsv(mods_path)

mods <- mods %>% filter(significant_x)

genes_df <- mods %>%
  count(model,feature=feature.x,relationship=if_else(sign(estimate.qnorm)>0,"pos","neg")) %>%
  pivot_wider(names_from = relationship, values_from = n,values_fill = 0)

te_df <- mods %>%
  count(model,feature=feature.y, relationship=if_else(sign(estimate.qnorm)>0,"pos","neg")) %>%
  pivot_wider(names_from = relationship, values_from = n,values_fill = 0)

res <- bind_rows(gene = genes_df, TE = te_df, .id = "feature.type")

g_n_tes_scatter <- res |>
  dplyr::rename(sex=model) |>
  filter(feature.type == "gene") |> 
  ggplot(aes(neg,pos)) +
  geom_point(size=rel(0.5)) +
  xlab("N negatively coexpressed TEs") +
  ylab("N positively coexpressed TEs") +
  facet_wrap(~sex, nrow=1) +
  labs(title="genes coexpressed with N TEs")

g_n_genes_scatter <- res |>
  dplyr::rename(sex=model) |>
  filter(feature.type == "TE") |>
  ggplot(aes(neg,pos)) +
  geom_point(size=rel(0.5)) +
  xlab("N negatively coexpressed genes") +
  ylab("N positively coexpressed genes") +
  facet_wrap(~sex, nrow = 1) +
  labs(title="TEs coexpressed with N genes")
                 
write_rds(g_n_genes_scatter,snakemake@output[["n_genes_scatter"]])
write_rds(g_n_tes_scatter,snakemake@output[["n_tes_scatter"]])
