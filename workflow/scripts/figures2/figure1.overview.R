library(tidyverse)
library(patchwork)
library(plotgardener)
library(DiagrammeRsvg)
library(rsvg)
library(vtree)
library(DiagrammeR)

# ------------------------------------------------------------------------------
# data req'd for plotting
# ------------------------------------------------------------------------------

lms <- read_tsv("upstream/final-models.collected-info.tsv.gz")

ncoex_scatter_and_hist = read_rds("results/exploratory_and_descriptive/ncoex_scatter_and_hist.rds")

piRNA_prop <- read_rds("results/pirna/pirna_coex_w_te_prop.rds")

# previously reported TE regulators
teregs <- read_tsv("results/resources/pirna_pathway.tsv") |> pull(gene_ID)


n_models = read_rds("results/exploratory_and_descriptive/n_models_per_filtering_step.rds")

# nofilt max and mean absolute coexpression of piRNA genes
# nofilt max and mean absolute coexpression of piRNA genes
nofilt_max_score <- read_tsv('results/rankings/nofilt_main_female_max_abs_estimate_qnorm.tsv.gz')
nofilt_mean_score <- read_tsv('results/rankings/nofilt_main_female_mean_abs_estimate_qnorm.tsv.gz')

prev_reported_scores <- bind_rows(`max abs. score`=nofilt_max_score, `mean abs. score`=nofilt_mean_score,.id = "score.type") |>
  mutate(prev.reported = gene_id %in% teregs)

# ------------------------------------------------------------------------------
# plot diagram
# ------------------------------------------------------------------------------

g_a <- grViz("
digraph boxes_and_circles {

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]

  # several 'node' statements
  node [shape = box,
        fontname = Arial]
  model [label='TEX ~ GEX + TE copies + gene/insertion proximity + wolbachia'];
  qnorm [label='quantile normalize coefficients\nby feature expression level'];
  filter [label='extract correlated\ngene/TE pairs'];
  quant [label='TE/gene expression\n(TEX/GEX)'];
  infect [label='wolbachia infection status'];
  cnv [label='WGS-based\nTE copy number estimate'];
  novel [label='novel TE regulating genes'];

  quant -> model;
  infect -> model;
  cnv -> model;
  model -> qnorm;
  qnorm -> filter;
  filter -> novel;
}
")

cartoon_temp <- tempfile()
#grVizToPNG(g_a, filename = "project_overview.png")
export_svg(g_a) |>
  charToRaw() |>
  rsvg_svg(cartoon_temp)

g_a_cartoon <- magick::image_read_svg(cartoon_temp) |> magick::image_ggplot()


# ------------------------------------------------------------------------------
# other figures
# ------------------------------------------------------------------------------

g_n_tes_scatter <- ncoex_scatter_and_hist |>
  dplyr::rename(sex=model) |>
  filter(feature.type == "gene") |> 
  ggplot(aes(neg,pos)) +
  geom_point(size=rel(0.5)) +
  xlab("N negatively coexpressed TEs") +
  ylab("N positively coexpressed TEs") +
  facet_wrap(~sex, nrow=1) +
  labs(title="genes coexpressed with N TEs")

g_n_genes_scatter <- ncoex_scatter_and_hist |>
  dplyr::rename(sex=model) |>
  filter(feature.type == "TE") |>
  ggplot(aes(neg,pos)) +
  geom_point(size=rel(0.5)) +
  xlab("N negatively coexpressed genes") +
  ylab("N positively coexpressed genes") +
  facet_wrap(~sex, nrow = 1) +
  labs(title="TEs coexpressed with N genes")

g_pirna <- piRNA_prop |> 
  mutate(data = map(data, ~rownames_to_column(.x, "class"))) |>
  unnest(data) |>
  mutate(proportion = n_coex/(n_coex + n_not_coex)) |>
  ggplot(aes(str_wrap(class,10),proportion)) +
  geom_col() +
  geom_text(data=piRNA_prop,aes(x=-Inf, y=1.1, label=label), hjust=0, size=rel(1.5)) +
  facet_wrap(~model) +
  coord_cartesian(ylim=c(0,1.3)) +
  xlab("gene type") +
  ylab("proportion >=1 significant TE correlation")

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



g_nhits_prev_rep_teregs <- lms |>
  filter(model == "female") |>
  group_by(feature.x,model) |>
  summarise(n_hits = sum(significant_x),.groups = "drop") |>
  mutate(prev.reported = feature.x %in% teregs) |>
  ggplot(aes(prev.reported,n_hits)) +
  geom_boxplot(outlier.shape = NA) +
  ggpubr::stat_compare_means(label.y = 9) +
  xlab("known TE regulator") +
  ylab("coexpressed TEs") +
  coord_cartesian(ylim=c(0,10))

g_score_prev_rep_teregs <- ggplot(prev_reported_scores, aes(prev.reported,value)) +
  geom_boxplot(outlier.shape = NA) +
  ggpubr::stat_compare_means(label.y = 0.4) +
  facet_wrap(~score.type, scales = "free_y") +
  ylab("coexpression score") +
  xlab("known TE regulator") +
  coord_cartesian(ylim=c(0,0.5))


# ------------------------------------------------------------------------------
# create page

theme_set(theme_classic() + 
            theme(text = element_text(size=7)) +
            theme(plot.title = element_text(hjust = 0.5))
          )

dir.create("results/figures2/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotText("A", x = 0.5, y=0.5)

plotGG(g_a_cartoon, x = 0.75, y=0.5, width = 3,height = 2.75)

plotGG(g_filt, x = 4.5, y=0.75, width = 2.75,height = 2.5)
plotText("B", x = 4.25, y=0.5)

plotGG(g_score_prev_rep_teregs, x = 0.5, y=3.5, width = 4,height = 2)
plotText("C",  x = 0.5, y=3.5)

plotGG(g_nhits_prev_rep_teregs, x = 5, y=3.5, width = 3,height = 2)
plotText("D",  x = 5, y=3.5)

plotGG(g_n_tes_scatter, x = 0.5, y=6, width = 3.5,height = 2)
plotText("E",  x = 0.5, y=6)

plotGG(g_n_genes_scatter, x = 4.25, y=6, width = 3.5,height = 2)
plotText("F",  x = 4.25, y=6)


dev.off()

