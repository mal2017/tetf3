Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(plotgardener)
library(tidyverse)

g_kd_check <- read_rds("results/deg/check_kds_by_chip_prox.gg.rds") +
  facet_wrap(~lab,ncol=5,scales="free") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

#------------------------------------------------------------------------------
# create page 1
# ------------------------------------------------------------------------------
theme_set(theme_classic() + 
            theme(text = element_text(size=7)) +
            theme(plot.title = element_text(hjust = 0.5))
)

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())
figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

plotGG(g_kd_check, x = .5, y=0.5, width = 7.5,height = 1.5)
plotText("A", x = .5, y=0.5)

dev.off()