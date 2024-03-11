library(plotgardener)
library(tidyverse)

g_kd_check <- read_rds("results/deg/check_kds_by_chip_prox.gg.rds") +
  facet_wrap(~lab,ncol=5,scales="free")

#------------------------------------------------------------------------------
# create page 1
# ------------------------------------------------------------------------------
theme_set(theme_classic() + 
            theme(text = element_text(size=7)) +
            theme(plot.title = element_text(hjust = 0.5))
)

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_kd_check, x = .5, y=0.5, width = 7.5,height = 3)
plotText("A", x = .5, y=0.5)

dev.off()