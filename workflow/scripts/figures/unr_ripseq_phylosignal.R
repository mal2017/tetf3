library(tidyverse)
library(plotgardener)

gtr <- read_rds("results/ripseq/unr_ripseq_phylosignal.tree.gg.rds")

g_crlg <- read_rds("results/ripseq/unr_ripseq_phylosignal.crlg.gg.rds")

# ------------------------------------------------------------------------------
# create page 2
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=7)) +
            theme(plot.title = element_text(hjust = 0.5)))

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

plotGG(gtr+guides(fill="none"), x = .5, y=0.5, width = 7,height = 5.5)
plotText("A", x = .5, y=0.5)

plotGG(g_crlg, x = 2, y=6.5, width = 4,height = 2)
plotText("B", x = 2, y=6.5)

dev.off()