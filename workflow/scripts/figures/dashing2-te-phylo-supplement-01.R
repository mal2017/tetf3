Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


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

# ------------------------------------------------------------------------------
# get
# ------------------------------------------------------------------------------

tr <- read_rds("results/te_sequence_similarity/te_sketch_tidytree.rds")


g <- ggtree(tr, layout = "rectangular", ) +
  geom_tiplab(size=rel(2))

g <- g + 
  geom_fruit(geom=geom_tile, aes(y=label, fill=repClass, position=1), pwidth = 0.01) + 
  scale_fill_brewer(palette = 1, type="qual")

g <- g + 
  new_scale_fill() + geom_fruit(geom=geom_tile, aes(y=label, fill=repFamily, position=2), pwidth = 0.01) +
  scale_fill_paletteer_d("pals::polychrome")


dir.create("results/figures/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())
figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

plotGG(g, x = 0.5, y=0.5, width = 7.5,height = 10)
plotText("A", x = 0.5, y=0.5)

dev.off()

writexl::write_xlsx(list(A=g$data),
                    path = ifelse(exists("snakemake"),
                                  snakemake@output$xlsx,
                                  "~/Downloads/test.xlsx"))