Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


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

g_b <- g_bcd$knockdown2_NfI_female_head_Mef2.R_control_female_head_Mef2.R
g_c <- g_bcd$knockdown2_pan_female_head_Mef2.R_control_female_head_Mef2.R
g_d <- g_bcd$knockdown2_Unr_female_head_Mef2.R_control_female_head_Mef2.R
g_e <- g_bcd$knockdown2_vvl_female_head_Mef2.R_control_female_head_Mef2.R
g_f <- g_bcd$knockdown2_CG16779_male_gonad_aTub_control_male_gonad_aTub

theme_set(theme_classic() + 
            theme(text = element_text(size=7))
)

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())
figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

plotGG(g_a, x = 0.5, y=0.5, width = 3.75,height = 2.5)
plotText("A", x = 0.5, y=0.5)

plotGG(g_b, 
       x = 4.5, y=0.5, width = 3.25,height = 2.5)
plotText("B", x = 4.5, y=0.5)

plotGG(g_c, x = 0.75, y=3.5, width = 3.25,height = 2.5)
plotText("C",  x = 0.5, y=3.5)

plotGG(g_d, x = 4.5, y=3.5, width = 3.25,height = 2.5)
plotText("D",  x = 4.5, y=3.5)


plotGG(g_e, x = 0.75, y=6.5, width = 3.25,height = 2.5)
plotText("E",  x = 0.5, y=6.5)

plotGG(g_f, x = 4.5, y=6.5, width = 3.25,height = 2.5)
plotText("F",  x = 4.5, y=6.5)


dev.off()


writexl::write_xlsx(list(A=g_a$data,
                         B=g_b$data,
                         C=g_c$data,
                         D=g_d$data,
                         E=g_e$data,
                         `F`=g_f$data),
                    path = ifelse(exists("snakemake"),
                                  snakemake@output$xlsx,
                                  "~/Downloads/test.xlsx"))