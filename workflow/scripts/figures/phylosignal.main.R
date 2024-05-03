Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)
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
library(ggpp)


# ------------------------------------------------------------------------------
# tf correlograms
# ------------------------------------------------------------------------------

crlg_gs <- read_rds("results/phylosignal/main_fig_correlograms.gg.rds")
ripseq_crlg_gs <- read_rds("results/ripseq/unr_ripseq_phylosignal.crlg.gg.rds")

g_b <- crlg_gs$random

g_c <- crlg_gs$bm

g_d <- crlg_gs$female_score_pan

g_e <- ripseq_crlg_gs

# ------------------------------------------------------------------------------
# get tree
# ------------------------------------------------------------------------------

tr <- read_rds("results/te_sequence_similarity/te_sketch_tidytree.rds")

tr <- tr |> mutate(repClass= if_else(label %in% c("Stalker3T","TLD2_LTR"),"LTR",repClass))

g_a <- ggtree(tr, layout = "rectangular",right = T, ladderize = T, size=0.25) +
  #geom_tiplab2(size=unit(1.1,"pt"),hjust = -0.5) + 
  geom_tippoint(aes(color=repClass)) +
  theme(#plot.margin = unit(c(-140, -140, -170, -140), "pt"),
        line = element_line(linewidth=0.001),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.position = c(0.65,0.45),legend.justification = c("left","top"), plot.background = element_blank())

# ------------------------------------------------------------------------------
# make page
# ------------------------------------------------------------------------------

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())
figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

plotGG(g_a, x = 0.5, y=0.5, width = 5,height = 2.5)
plotText("A", x = 0.5, y=0.5)

plotGG(g_b, x = 5.5, y=0.75, width = 2.5,height = 2)
plotText("B", x = 5.5, y=0.5)

plotGG(g_c, x = .5, y=3.25, width = 2.25,height = 2)
plotText("C", x = .5, y=3.25)

plotGG(g_d, x = 3.25, y=3.25, width = 2.25,height = 2)
plotText("D", x = 3.25, y=3.25)

plotGG(g_e, x = 5.75, y=3.25, width = 2.25,height = 2)
plotText("E", x = 5.75, y=3.25)

dev.off()

writexl::write_xlsx(list(A=g_a$data,
                         B=g_b$data,
                         C=g_c$data,
                         D=g_d$data,
                         E=g_e$data),
                    path = ifelse(exists("snakemake"),
                                  snakemake@output$xlsx,
                                  "~/Downloads/test.xlsx"))