library(patchwork)
library(plotgardener)
library(tidyverse)
library(clusterProfiler)
library(ggpubr)
library(phylosignal)
library(phylobase)
library(ggplotify)
library(ggtree)
library(tidytree)
library(ggtreeExtra)
library(ggnewscale)
library(paletteer)
library(ggdensity)

# normalization
mods <- read_tsv("upstream/final-models.collected-info.tsv.gz")

library(ggridges)
# mods %>%
#   mutate(across(contains("mean_qtiled.gene"),as.factor)) %>%
#   ggplot(aes(x=estimate.qnorm,y=as.factor(mean_qtile.gene),fill=mean_qtile.gene)) +
#   geom_density_ridges2() +
#   #scale_fill_cyclical(values = c("darkgray", "lightgray")) +
#   coord_cartesian(xlim=c(-1,1))
# 
# mods %>%
#   mutate(across(contains("mean_qtiled.gene"),as.factor)) %>%
#   ggplot(aes(x=estimate,y=as.factor(mean_qtile.te),fill=mean_qtile.gene)) +
#   geom_density_ridges2() +
#   #scale_fill_cyclical(values = c("darkgray", "lightgray")) +
#   coord_cartesian(xlim=c(-1,1))

g_qnorm <- mods |>
  ggplot(aes(estimate.qnorm, estimate)) +
  ggdensity::geom_hdr() +
  xlab("quantile-normalized coexpression score") +
  ylab("raw coexpression score")

rm(mods); gc()

# replicate correlation
repcor <- read_rds("results/coexpression_replication/intermediate/replicate_dataset_correlation.rds")
g_repcor <- filter(repcor, result_set=="main_data") |>
  dplyr::select(model,gg) |>
  deframe() |>
  imap(~{.x + labs(title=.y)}) |>
  Reduce("+",x=_)


# male vs female
mf <- read_rds("results/coexpression_replication/intermediate/mf_dataset_correlation.rds")
g_mf <- mf |>
  filter(result_set %in% c("unfiltered","replicated")) |>
  dplyr::select(result_set,gg) |>
  deframe() |>
  imap(~{.x + labs(title=.y)}) |>
  Reduce("+",x=_)


# variance explained
g_var_exp <- read_rds("results/exploratory_and_descriptive/g_variance_explained.rds") + 
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle=25, hjust=1))


# plotting --------------------------------------------------------------------
theme_set(theme_classic() + theme(text=element_text(size=unit(5,"pt"))))

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 8.5, height = 11)

pageCreate(width =8.5, height = 11, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = g_qnorm, x = 0.25, y=0.05, width = 3.75, height=3)

pd <- plotGG(plot = g_var_exp, x = 4.1, y=0.05, width = 4, height=3.5)

pb <- plotGG(plot = g_repcor, x = 0.25,  y=3.75, width = 8, height=3)

pc <- plotGG(plot = g_mf, x = 0.25,  y=6.75, width = 8, height=3)


plotText(label = "A", fontsize = 7,
         x = 0.25, y = 0.25, just = "center", default.units = "inches")

plotText(label = "D", fontsize = 7,
         x = 4.25, y = 0.25, just = "center", default.units = "inches")

plotText(label = "B", fontsize = 7,
         x = 0.25, y = 3.75, just = "center", default.units = "inches")

plotText(label = "C", fontsize = 7,
         x = 0.25, y = 7, just = "center", default.units = "inches")


dev.off()
