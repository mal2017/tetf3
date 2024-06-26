Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(ggtext)
library(patchwork)
library(plotgardener)

# ------------------------------------------------------------------------------
# get de results
# ------------------------------------------------------------------------------

g_a <- read_rds("results/deg/de_volcanos.gg.rds")

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

plotGG(g_a, x = 0.5, y=0.5, width = 7.5,height = 5)

dev.off()

writexl::write_xlsx(list(A=g_a$data),
                    path = ifelse(exists("snakemake"),
                                  snakemake@output$xlsx,
                                  "~/Downloads/test.xlsx"))