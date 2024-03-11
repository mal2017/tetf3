library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(ggtext)
library(patchwork)
library(plotgardener)

# ------------------------------------------------------------------------------
# get de results
# ------------------------------------------------------------------------------

g_a <- read_rds("results/deg/de_volcanos.gg.rds")

# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=7))
)

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_a, x = 0.5, y=0.5, width = 7.5,height = 10)

dev.off()

