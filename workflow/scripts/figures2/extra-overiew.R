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

# previously reported TE regulators
teregs <- read_tsv("results/resources/pirna_pathway.tsv") |> pull(gene_ID)


# nofilt max and mean absolute coexpression of piRNA genes
nofilt_max_score <- read_tsv('results/rankings/nofilt_main_male_max_abs_estimate_qnorm.tsv.gz')
nofilt_mean_score <- read_tsv('results/rankings/nofilt_main_male_mean_abs_estimate_qnorm.tsv.gz')

# male only
prev_reported_scores <- bind_rows(`max abs. score`=nofilt_max_score, `mean abs. score`=nofilt_mean_score,.id = "score.type") |>
  mutate(prev.reported = gene_id %in% teregs)


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

g_nhits_prev_rep_teregs <- mods |>
  filter(model == "male") |>
  group_by(feature.x,model) |>
  summarise(n_hits = sum(significant_x),.groups = "drop") |>
  mutate(prev.reported = feature.x %in% teregs) |>
  ggplot(aes(prev.reported,n_hits)) +
  geom_boxplot(outlier.shape = NA) +
  ggpubr::stat_compare_means(label.y = 9, size=2) +
  xlab("known TE regulator") +
  ylab("coexpressed TEs") +
  coord_cartesian(ylim=c(0,10))

g_score_prev_rep_teregs <- ggplot(prev_reported_scores, aes(prev.reported,value)) +
  geom_boxplot(outlier.shape = NA) +
  ggpubr::stat_compare_means(label.y = 0.45,size=2) +
  facet_wrap(~score.type, scales = "free_y") +
  ylab("coexpression score") +
  xlab("known TE regulator") +
  coord_cartesian(ylim=c(0,0.5))


g_sharedness <- mods |> 
  dplyr::select(feature.x,feature.y,model,significant_x) |>
  pivot_wider(names_from = model, values_from = significant_x) |>
  filter(male | female) |>
  mutate(status = map2_chr(male,female, ~case_when(
    .x & .y ~ "shared",
    .x & !.y ~ "male-private",
    .y & !.x ~ "female-private",
    !(.y | .x) ~ "n.s.",
    T ~ "wtf"
  ))) |>
  ggplot(aes(status)) +
  geom_bar()



# plotting --------------------------------------------------------------------
theme_set(theme_classic() + theme(text=element_text(size=unit(7,"pt"))))

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 8.5, height = 11)

pageCreate(width =8.5, height = 11, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = g_var_exp, x = 0.5, y=0.5, width = 7.5, height=2.75)
plotText(label = "A", x = 0.5, y = 0.5)

pb <- plotGG(plot = g_mf$unfiltered+ scale_fill_distiller(palette = 6) + guides(fill="none"), x = 0.4,  y=3.75, width = 1.75, height=1.75)
pc <- plotGG(plot = g_mf$replicated+ scale_fill_distiller(palette = 6) + guides(fill="none"), x = 2.3,  y=3.75, width = 1.75, height=1.75)
pd <- plotGG(plot = g_repcor$male + scale_fill_distiller(palette = 6) + guides(fill="none"), x = 4.3,  y=3.75, width = 1.75, height=1.75)
pe <- plotGG(plot = g_repcor$female+ scale_fill_distiller(palette = 6) + guides(fill="none"), x = 6.3,  y=3.75, width = 1.75, height=1.75)


plotText(label = "B", x = 0.5, y = 3.75)
plotText(label = "C", x = 2.3, y = 3.75)
plotText(label = "D", x = 4.3, y = 3.75)
plotText(label = "E", x = 6.3, y = 3.75)

plotGG(g_sharedness, x = 0.5, y=5.75, width = 2.5,height = 2)
plotText("F",  x = 0.5, y=5.75)

plotGG(g_score_prev_rep_teregs, x = 2.65, y=5.75, width = 2.75,height = 2)
plotText("G",  x = 2.65, y=5.75)

plotGG(g_nhits_prev_rep_teregs, x = 5.5, y=5.75, width = 2.5,height = 2)
plotText("H",  x = 5.5, y=5.75)

dev.off()
