Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


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

g_a <- g_bcdef$knockdown2_NfI_female_head_Mef2.R_control_female_head_Mef2.R
plotGG(g_a,
       x = 0.5, y=0.5, width = 2.25,height = 2)
plotText("A", x = 0.5, y=0.5)

g_b <- g_bcdef$knockdown2_vvl_female_head_Mef2.R_control_female_head_Mef2.R
plotGG(g_b, 
       x = 3, y=0.5, width = 2.25,height = 2)
plotText("B", x = 3, y=0.5)

g_c <- g_bcdef$knockdown2_Unr_female_head_Mef2.R_control_female_head_Mef2.R
plotGG(g_c, 
       x = 5.5, y=0.5, width = 2.25,height = 2)
plotText("C", x = 5.5, y=0.5)

g_d <- g_bcdef$knockdown2_pan_female_head_Mef2.R_control_female_head_Mef2.R
plotGG(g_d,
       x = 0.5, y=2.75, width = 2.25,height = 2)
plotText("D", x = 0.5, y=2.75)

g_e <- g_bcdef$knockdown2_CG16779_female_head_Mef2.R_control_female_head_Mef2.R
plotGG(g_e, 
       x = 3, y=2.75, width = 2.25,height = 2)
plotText("E", x = 3, y=2.75)

g_f <- g_bcdef$knockdown2_CG16779_male_gonad_aTub_control_male_gonad_aTub
plotGG(g_f, 
       x = 5.5, y=2.75, width = 2.25,height = 2)
plotText("F", x = 5.5, y=2.75)


# ------------------------------------------------------------------------------
# sirna

g_g <- g_hijkl$knockdown2_NfI_female_head_Mef2.R_control_female_head_Mef2.R
plotGG(g_g,
       x = 0.5, y=5, width = 2.25,height = 2)
plotText("G", x = 0.5, y=5)

g_h <- g_hijkl$knockdown2_vvl_female_head_Mef2.R_control_female_head_Mef2.R
plotGG(g_h, 
       x = 3, y=5, width = 2.25,height = 2)
plotText("H", x = 3, y=5)

g_i <- g_hijkl$knockdown2_Unr_female_head_Mef2.R_control_female_head_Mef2.R
plotGG(g_i, 
       x = 5.5, y=5, width = 2.25,height = 2)
plotText("I", x = 5.5, y=5)

g_j <- g_hijkl$knockdown2_pan_female_head_Mef2.R_control_female_head_Mef2.R
plotGG(g_j,
       x = 0.5, 7, width = 2.25,height = 2)
plotText("J", x = 0.5, y=7)

g_k <- g_hijkl$knockdown2_CG16779_female_head_Mef2.R_control_female_head_Mef2.R
plotGG(g_k, 
       x = 3, y=7, width = 2.25,height = 2)
plotText("K", x = 3, y=7)

g_l <- g_hijkl$knockdown2_CG16779_male_gonad_aTub_control_male_gonad_aTub
plotGG(g_l, 
       x = 5.5, y=7, width = 2.25,height = 2)
plotText("L", x = 5.5, y=7)

dev.off()

writexl::write_xlsx(list(A=g_a$data,
                         B=g_b$data,
                         C=g_c$data,
                         D=g_d$data,
                         E=g_e$data,
                         `F`=g_f$data,
                         G=g_g$data,
                         H=g_h$data,
                         I=g_i$data,
                         J=g_j$data,
                         K=g_k$data,
                         L=g_l$data),
                    path = ifelse(exists("snakemake"),
                                  snakemake@output$xlsx,
                                  "~/Downloads/test.xlsx"))
