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


# ------------------------------------------------------------------------------
# create page

theme_set(theme_classic() + 
            theme(text = element_text(size=7)) +
            theme(plot.title = element_text(hjust = 0.5))
          )

dir.create("results/figures2/")

pdf("results/figures2/figure1.pdf",width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotText("A", x = 0.5, y=0.5)

plotGG(g_a_cartoon, x = 0.75, y=0.5, width = 3,height = 2.75)

plotGG(g_pirna, x = 4.75, y=0.75, width = 2.5,height = 2.1)
plotText("B", x = 4.25, y=0.5)

plotGG(g_n_tes_scatter, x = 0.5, y=3.5, width = 3.5,height = 2.25)
plotText("C",  x = 0.5, y=3.5)


plotGG(g_n_genes_scatter, x = 4.75, y=3.5, width = 3.5,height = 2.25)
plotText("D",  x = 4.75, y=3.5)

dev.off()

