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
library(patchwork)

# size of the supercells found
g_supercell_size <- read_rds("results/calderon22/g_supercell_size.rds")

g_supercell_size <- g_supercell_size + xlab("barcodes per supercell")

gg_spqn <- read_rds("results/calderon22/gg_spqn.rds")

g_spqn <- (gg_spqn$gg_coefs_raw[[1]] + guides(fill="none")) +
  gg_spqn$gg_coefs_spqn[[1]] & 
  plot_layout(guides="collect") &
  theme(text=element_text(size=5), axis.text = element_text(size=5))

g_all_cell_overlapping_features_coexpressed <- read_rds("results/calderon22/g_all_cell_overlapping_features_coexpressed.rds")
g_pan_highly_corr_with_tes <- read_rds("results/calderon22/g_pan_highly_corr_with_tes.rds")

g_pan_highly_corr_with_tes <- g_pan_highly_corr_with_tes + theme(legend.position = "bottom")

g_poscon_and_hits_panel <- read_rds("results/calderon22/g_poscon_and_hits_panel.rds")
g_negcon_panel <- read_rds("results/calderon22/g_negcon_panel.rds")


# plotting page 1 --------------------------------------------------------------------
theme_set(theme_classic() + theme(text=element_text(size=unit(5,"pt"))))

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 8.5, height = 11)

pageCreate(width =8.5, height = 11, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = g_supercell_size, x = 0.5, y=0.5, width = 2.2, height=2)

pb1 <- plotGG(plot = g_spqn, x = 3,  y=0.5, width = 5, height=2)

pc <- plotGG(plot = g_all_cell_overlapping_features_coexpressed, x = 0.5, y=3, width = 4, height=2)

pd <- plotGG(plot = g_pan_highly_corr_with_tes, x = 4.75,  y=3, width = 3, height=2)

plotText(label = "A", fontsize = 7, x = 0.5, y = 0.5)

plotText(label = "B", fontsize = 7, x = 3, y = 0.5)

plotText(label = "C", fontsize = 7, x = 0.5, y = 3)

plotText(label = "D", fontsize = 7, x = 4.75, y = 3)

dev.off()

# plotting page 2 --------------------------------------------------------------------

if (!interactive()) pdf(snakemake@output[["pdf2"]],width = 8.5, height = 11)

pageCreate(width =8.5, height = 11, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = g_poscon_and_hits_panel, x = 0.5,  y=0.5, width = 7.5, height=4.25)

pb <- plotGG(plot = g_negcon_panel, x = 0.5,  y=5.25, width = 7.5, height=4.25)

plotText(label = "A", fontsize = 7, x = 0.5, y = 0.5)

plotText(label = "B", fontsize = 7, x = 0.5, y = 5)

dev.off()