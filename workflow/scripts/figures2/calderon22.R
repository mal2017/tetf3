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
g_supercell_size <- read_rds("results/calderon22/g_supercell_size.rds") +
  xlab("barcodes per supercell")

g_pan_highly_corr_with_tes <- read_rds("results/calderon22/g_pan_highly_corr_with_tes.rds") +
  theme(legend.position = "bottom") +
  scale_fill_manual(values=c("pan"="red","other"="darkgray"),name="gene")

g_pan_vs <- read_rds("results/calderon22/g_pan_vs.rds")

# plotting page 1 --------------------------------------------------------------------
theme_set(theme_classic() + theme(text=element_text(size=unit(5,"pt"))))

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 8.5, height = 11)

pageCreate(width =8.5, height = 11, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = g_supercell_size, x = 0.5, y=0.5, width = 3, height=2.25)
plotText(label = "A", x = 0.5, y = 0.5)

pb <- plotGG(plot = g_pan_highly_corr_with_tes, x = 4,  y=0.5, width = 3.5, height=2.25)
plotText(label = "B", x = 4, y = 0.5)

plotGG(g_pan_vs$nkd + 
         g_pan_vs$`1360` + 
         g_pan_vs$invader2 + 
         plot_layout(guides="collect") & 
         theme(title=element_text(size=7),legend.position = "bottom"), x=0.5, y=3.5, width=7.5, height = 3.5)


dev.off()