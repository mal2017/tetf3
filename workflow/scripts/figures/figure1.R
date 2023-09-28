library(tidyverse)
library(patchwork)
library(plotgardener)
library(DiagrammeRsvg)
library(rsvg)
library(vtree)
library(DiagrammeR)

lms <- read_tsv("upstream/final-models.collected-info.tsv.gz")

coex_dist = read_rds("results/exploratory_and_descriptive/coex_metric_distributions.rds")
n_models = read_rds("results/exploratory_and_descriptive/n_models_per_filtering_step.rds")
zads <- read_tsv("results/resources/zad_genes.tsv")$ensembl_gene_id
ncoex_scatter_and_hist = read_rds("results/exploratory_and_descriptive/ncoex_scatter_and_hist.rds")
piRNA_prop <- read_rds("results/pirna/pirna_coex_w_te_prop.rds")

g_b <- n_models |>
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

g_c <- lms |>
  group_by(model, feature.x) |>
  summarise(n=sum(significant_x)) |>
  ggplot(aes(n)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~model) +
  coord_flip() +
  xlab("N coexpressed TEs") +
  ylab("N genes")

g_c_label <- g_c$data |>
  filter(n>0) |>
  count(model) |>
  mutate(label = sprintf("\t\tN genes coexpressed\n\t\twith >=1 TE: %s",n))

g_c <- g_c +
  geom_text(data = g_c_label, aes(label=label, x=30, y=1000),size=rel(1.5))



g_d <- ncoex_scatter_and_hist |>
  dplyr::rename(sex=model) |>
  filter(feature.type == "gene") |> 
  ggplot(aes(neg,pos)) +
  geom_point(size=rel(0.5)) +
  xlab("N negatively coexpressed TEs") +
  ylab("N positively coexpressed TEs") +
  facet_wrap(~sex, nrow=1) +
  labs(title="genes coexpressed with N TEs")

g_e <- ncoex_scatter_and_hist |>
  dplyr::rename(sex=model) |>
  filter(feature.type == "TE") |>
  ggplot(aes(neg,pos)) +
  geom_point(size=rel(0.5)) +
  xlab("N negatively coexpressed genes") +
  ylab("N positively coexpressed genes") +
  facet_wrap(~sex, nrow = 1) +
  labs(title="TEs coexpressed with N genes")

g_f <- piRNA_prop |> 
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

g_a <- grViz("
digraph boxes_and_circles {

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]

  # several 'node' statements
  node [shape = box,
        fontname = Arial]
  model [label='TEX ~ GEX + TE copies + gene/insertion proximity + wolbachia'];
  qnorm [label='quantile normalize GEX coeficients\nby feature expression level'];
  filter [label='extract correlated\ngene/TE pairs'];
  quant [label='TE/gene expression\n(TEX/GEX)'];
  infect [label='wolbachia infection status'];
  cnv [label='WGS-based\nTE copy number estimate'];
  novel [label='novel TE regulating genes'];
  enrich [label='enriched gene groups'];
  rewire [label='TE-mediated \ngene regulation changes'];

  quant -> model;
  infect -> model;
  cnv -> model;
  model -> qnorm;
  qnorm -> filter;
  filter -> novel;
  filter -> enrich;
  filter -> rewire;
}
")

cartoon_temp <- tempfile()
#grVizToPNG(g_a, filename = "project_overview.png")
export_svg(g_a) |>
  charToRaw() |>
  rsvg_svg(cartoon_temp)

g_a_cartoon <- magick::image_read_svg(cartoon_temp) |> magick::image_ggplot()

# ------------------------------------------------------------------------------
# create page

theme_set(theme_classic() + 
            theme(text = element_text(size=7)) +
            theme(plot.title = element_text(hjust = 0.5))
          )

dir.create("results/figures/")

pdf("results/figures/figure1.pdf",width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotText("A", x = 0.5, y=0.5)

plotGG(g_a_cartoon, x = 0.5, y=0.5, width = 2.75,height = 2)

plotGG(g_b, x = 3.3, y=0.5, width = 2.5,height = 2.1)
plotText("B", x = 3.5, y=0.5)

plotText("C",  x = 5.95, y=0.5)
plotGG(g_c, x = 6, y=0.5, width = 2,height = 2.1)


plotText("D",  x = 0.5, y=2.85)
plotGG(g_d, x = 0.5, y=3, width = 2.5,height = 1.75)

plotText("E" , x = 3.25, y=2.85)
plotGG(g_e, x = 3, y=3, width = 2.5,height = 1.75)

plotText("F", x = 5.625, y=2.85)

plotGG(g_f, x = 5.5, y=3.125, width = 2.5,height = 2)

dev.off()

