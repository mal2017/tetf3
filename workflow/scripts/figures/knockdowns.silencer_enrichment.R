library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(plotgardener)

# ------------------------------------------------------------------------------
# rws for te regs and sirna
# ------------------------------------------------------------------------------
rws <- read_rds("results/signatures/ourKD_gsea_randomwalks.gg_df.rds")

teregs <- filter(rws,signature_name == "TE.regulators")
sirna <- filter(rws,signature_name == "siRNA")

g_bcdef <- teregs |>
  dplyr::select(comparison, gg) |>
  deframe() |>
  map( ~{ .x & theme(axis.title = element_text(size=5), axis.text = element_text(size=5), plot.title = element_text(size=7, hjust=0.5))})

g_hijkl <- teregs |>
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
figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

# ------------------------------------------------------------------------------
# te regs
plotGG(g_bcdef$knockdown2_NfI_female_head_Mef2.R_control_female_head_Mef2.R,
       x = 0.5, y=0.5, width = 2.25,height = 2)
plotText("A", x = 0.5, y=0.5)

plotGG(g_bcdef$knockdown2_vvl_female_head_Mef2.R_control_female_head_Mef2.R, 
       x = 3, y=0.5, width = 2.25,height = 2)
plotText("B", x = 3, y=0.5)

plotGG(g_bcdef$knockdown2_Unr_female_head_Mef2.R_control_female_head_Mef2.R, 
       x = 5.5, y=0.5, width = 2.25,height = 2)
plotText("C", x = 5.5, y=0.5)

plotGG(g_bcdef$knockdown2_pan_female_head_Mef2.R_control_female_head_Mef2.R,
       x = 0.5, y=2.75, width = 2.25,height = 2)
plotText("D", x = 0.5, y=2.75)

plotGG(g_bcdef$knockdown2_CG16779_female_head_Mef2.R_control_female_head_Mef2.R, 
       x = 3, y=2.75, width = 2.25,height = 2)
plotText("E", x = 3, y=2.75)

plotGG(g_bcdef$knockdown2_CG16779_male_gonad_aTub_control_male_gonad_aTub, 
       x = 5.5, y=2.75, width = 2.25,height = 2)
plotText("F", x = 5.5, y=2.75)


# ------------------------------------------------------------------------------
# sirna

plotGG(g_hijkl$knockdown2_NfI_female_head_Mef2.R_control_female_head_Mef2.R,
       x = 0.5, y=5, width = 2.25,height = 2)
plotText("G", x = 0.5, y=5)

plotGG(g_hijkl$knockdown2_vvl_female_head_Mef2.R_control_female_head_Mef2.R, 
       x = 3, y=5, width = 2.25,height = 2)
plotText("H", x = 3, y=5)

plotGG(g_hijkl$knockdown2_Unr_female_head_Mef2.R_control_female_head_Mef2.R, 
       x = 5.5, y=5, width = 2.25,height = 2)
plotText("I", x = 5.5, y=5)

plotGG(g_hijkl$knockdown2_pan_female_head_Mef2.R_control_female_head_Mef2.R,
       x = 0.5, 7, width = 2.25,height = 2)
plotText("J", x = 0.5, y=7)

plotGG(g_hijkl$knockdown2_CG16779_female_head_Mef2.R_control_female_head_Mef2.R, 
       x = 3, y=7, width = 2.25,height = 2)
plotText("K", x = 3, y=7)

plotGG(g_hijkl$knockdown2_CG16779_male_gonad_aTub_control_male_gonad_aTub, 
       x = 5.5, y=7, width = 2.25,height = 2)
plotText("L", x = 5.5, y=7)

dev.off()

