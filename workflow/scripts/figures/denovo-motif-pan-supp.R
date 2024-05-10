Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)
library(memes)
library(universalmotif)
library(gridExtra)
library(plotgardener)
library(patchwork)


# ------------------------------------------------------------------------------
g_denovo_motifs <- read_rds("results/motifs/comparison/pan_within_denovo.gg.rds")
g_denovo_motifs <- g_denovo_motifs$CRCAKKSMCRARRAS

g_meme_boxplot <- read_rds("results/motifs/n_denovo_vs_sig_coef.pan.rds") |>
  filter(str_detect(motif,"CRCAK")) |>
  pull(g_boxplot) |>
  pluck(1)

g_upstream_sea <- read_rds("results/motifs/upstream_csem_known_pan_sea.gg.rds")

# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

g_a <- g_upstream_sea
g_b <- g_meme_boxplot
g_c <- g_denovo_motifs


theme_set(theme_classic() + 
            theme(text = element_text(size=7))
)

dir.create("results/figures/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())
figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

plotGG(g_a, x = 0.5, y=0.5, width = 7.5,height = 2.75)
plotText("A", x = 0.5, y=.5)

plotGG(g_b, x=0.75, y=3.75, width = 2, height=2.25)
plotText("B", x = 0.5, y=3.75)

plotGG(g_c, x = 3, y=3.75, width = 5,height = 2.8)
plotText("C", x = 3, y=3.75)

dev.off()

writexl::write_xlsx(list(A=g_a$data,
                         B=g_b$data,
                         C=g_c$data),
                    path = ifelse(exists("snakemake"),
                                  snakemake@output$xlsx,
                                  "~/Downloads/test.xlsx"))