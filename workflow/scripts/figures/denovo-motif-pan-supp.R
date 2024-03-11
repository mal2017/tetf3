library(tidyverse)
library(memes)
library(universalmotif)
library(gridExtra)
library(plotgardener)
library(patchwork)


# ------------------------------------------------------------------------------
g_denovo_motifs <- read_rds("results/motifs/comparison/pan_within_denovo.gg.rds")

g_meme_boxplot <- read_rds("results/motifs/n_denovo_vs_sig_coef.pan.rds") |>
  filter(str_detect(motif,"MARY")) |>
  pull(g_boxplot) |>
  pluck(1)

g_upstream_sea <- read_rds("results/motifs/upstream_csem_known_pan_sea.gg.rds")

# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=7))
)

dir.create("results/figures/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_meme_boxplot, x=0.75, y=0.75, width = 2, height=2.25)
plotText("A", x = 0.5, y=0.5)

plotGG(g_denovo_motifs, x = 3, y=0.5, width = 5,height = 2.8)
plotText("B", x = 3, y=0.5)

plotGG(g_upstream_sea, x = 0.5, y=3.5, width = 7.5,height = 2.75)
plotText("C", x = 0.5, y=3.5)

dev.off()