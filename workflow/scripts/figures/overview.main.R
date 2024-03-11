library(tidyverse)
library(patchwork)
library(plotgardener)
library(DiagrammeRsvg)
library(rsvg)
library(vtree)
library(DiagrammeR)

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
  rsvg_svg(cartoon_temp,width = 2000,height = 2000)

g_a_cartoon <- magick::image_read_svg(cartoon_temp) |> magick::image_ggplot(interpolate = T)


# ------------------------------------------------------------------------------
# other figures
# ------------------------------------------------------------------------------

g_n_tes_scatter <- read_rds('results/exploratory_and_descriptive/n_tes_scatter.gg.rds')

g_n_genes_scatter <- read_rds('results/exploratory_and_descriptive/n_genes_scatter.gg.rds')

g_pirna <- read_rds("results/pirna/te_silencer_overrepresentation.gg.rds")

g_filt <- read_rds("results/exploratory_and_descriptive/n_models_per_filtering_step.gg.rds")

g_nhits_prev_rep_teregs <- read_rds("results/pirna/te_silencer_n_hits_boxplot.females.gg.rds")

g_score_prev_rep_teregs <- read_rds("results/pirna/te_silencer_scores_boxplot.females.gg.rds")


# ------------------------------------------------------------------------------
# create page

theme_set(theme_classic() + 
            theme(text = element_text(size=7)) +
            theme(plot.title = element_text(hjust = 0.5))
          )

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotText("A", x = 0.5, y=0.5)

plotGG(g_a_cartoon, x = 0.75, y=0.4, width = 3.5,height = 3)

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

