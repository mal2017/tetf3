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
library(ggridges)

# ------------------------------------------------------------------------------
# read in data
# ------------------------------------------------------------------------------

mods <- read_tsv("upstream/final-models.collected-info.tsv.gz")
coex_dist = read_rds("results/exploratory_and_descriptive/coex_metric_distributions.rds")
n_models = read_rds("results/exploratory_and_descriptive/n_models_per_filtering_step.rds")
repcor <- read_rds("results/coexpression_replication/replicate_dataset_correlation.rds")
mf <- read_rds("results/coexpression_replication/mf_dataset_correlation.rds")
var_exp <- read_rds("results/exploratory_and_descriptive/g_variance_explained.rds") 

# ------------------------------------------------------------------------------
# plotting
# ------------------------------------------------------------------------------

# filtering
g_filt <- n_models |>
  dplyr::rename(sex=model) |>
  mutate(subset = str_wrap(subset,width=20)) |>
  mutate(subset = fct_reorder(subset, n)) |> #pull(subset) %>% .[12]
  ggplot(aes(n, subset,fill=sex)) +
  geom_col(position = "dodge") +
  geom_text(data = \(x) filter(x, subset == subset[which.min(n)]), 
            aes(label=paste0("n=", n), x=300000),
            size=rel(2),
            position = position_dodge(width = 0.75)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) +
  ylab("filtering step") +
  xlab("N TE/gene pairs") +
  scale_fill_grey() +
  theme(legend.position = c(1,0), legend.justification = c("right","bottom"))



# replicate correlation
g_repcor <- filter(repcor, result_set=="main_data") |>
  dplyr::select(model,gg) |>
  deframe() |>
  imap(~{.x + labs(title=.y)})

# male vs female
g_mf <- mf |>
  filter(result_set %in% c("unfiltered","replicated")) |>
  dplyr::select(result_set,gg) |>
  deframe() |>
  imap(~{.x + labs(title=.y)})

# variance explained
g_var_exp <- var_exp + 
  #theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle=25, hjust=1)) +
  scale_fill_grey(start = 0.4, end=0.8)


# plotting --------------------------------------------------------------------
theme_set(theme_classic() + theme(text=element_text(size=unit(7,"pt"))))

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 8.5, height = 11)

pageCreate(width =8.5, height = 11, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = g_filt, x = 0.5, y=0.5, width = 3, height=2.75)
plotText(label = "A", x = 0.5, y = 0.5)

pb <- plotGG(plot = g_var_exp, x = 4, y=0.5, width = 4.5, height=2.75)
plotText(label = "B", x = 4, y = 0.5)



pc <- plotGG(plot = g_repcor$male + scale_fill_distiller(palette = 6) + guides(fill="none"), x = 0.4,  y=3.75, width = 1.75, height=1.75)
pd <- plotGG(plot = g_repcor$female+ scale_fill_distiller(palette = 6) + guides(fill="none"), x = 2.3,  y=3.75, width = 1.75, height=1.75)
pe <- plotGG(plot = g_mf$unfiltered+ scale_fill_distiller(palette = 6) + guides(fill="none"), x = 4.3,  y=3.75, width = 1.75, height=1.75)
pf <- plotGG(plot = g_mf$replicated+ scale_fill_distiller(palette = 6) + guides(fill="none"), x = 6.3,  y=3.75, width = 1.75, height=1.75)


plotText(label = "C", x = 0.5, y = 3.75)
plotText(label = "D", x = 2.3, y = 3.75)
plotText(label = "E", x = 4.3, y = 3.75)
plotText(label = "F", x = 6.3, y = 3.75)


dev.off()
