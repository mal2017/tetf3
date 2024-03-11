library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(ggtext)
library(patchwork)
library(plotgardener)
library(tidyverse)

# ------------------------------------------------------------------------------
# target-specific all-te hit barchart
# ------------------------------------------------------------------------------

g_a <- read_rds("results/signatures/ourKD_gsea_barplots.gg_list.rds")$all_tes

# ------------------------------------------------------------------------------
# exemplary all-te random walk
# ------------------------------------------------------------------------------
gg <- read_rds("results/signatures/ourKD_gsea_randomwalks.gg_df.rds")
g_bcd <- gg |>
  filter(signature_name == "all_tes") |>
  dplyr::select(comparison, gg) |>
  deframe() |>
  map( ~{ .x & theme(axis.title = element_text(size=5), axis.text = element_text(size=5), plot.title = element_text(size=7, hjust=0.5))})

# ------------------------------------------------------------------------------
# create page
  # ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=7))
)

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_a, x = 0.5, y=0.5, width = 3.75,height = 2.5)
plotText("A", x = 0.5, y=0.5)

plotGG(g_bcd$knockdown2_NfI_female_head_Mef2.R_control_female_head_Mef2.R, 
       x = 4.5, y=0.5, width = 3.25,height = 2.5)
plotText("B", x = 4.5, y=0.5)

plotGG(g_bcd$knockdown2_pan_female_head_Mef2.R_control_female_head_Mef2.R, x = 0.75, y=3.5, width = 3.25,height = 2.5)
plotText("C",  x = 0.5, y=3.5)

plotGG(g_bcd$knockdown2_Unr_female_head_Mef2.R_control_female_head_Mef2.R, x = 4.5, y=3.5, width = 3.25,height = 2.5)
plotText("D",  x = 4.5, y=3.5)

dev.off()

