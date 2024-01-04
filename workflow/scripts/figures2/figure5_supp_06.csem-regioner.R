library(tidyverse)
library(plotgardener)
library(patchwork)
library(GenomicFeatures)
library(plyranges)

# function for nicely plotting randomization results
plot_regioner <- function(x) {
  permuted_overlaps <- tibble(overlaps = x$numOverlaps$permuted)
  observed_overlaps <- tibble(overlaps = x$numOverlaps$observed)
  pval <- x$numOverlaps$pval
  z <- x$numOverlaps$zscore
  alt <- x$numOverlaps$alternative
  
  ggplot(permuted_overlaps, aes(overlaps)) +
    geom_histogram() +
    geom_vline(data = observed_overlaps, aes(xintercept = overlaps), color="darkgreen") +
    annotate("text", x=-Inf, y=Inf, label=sprintf("z=%s, p=%s, alternative=%s", format(z, digits=3), format.pval(pval,digits = 2), alt), hjust=0, vjust=1)
}

dat <- read_rds("results/csem_mosaics/regioner.rds")

pan_all_te <- plot_regioner(dat$pan_all_te)
pan_pan_te <- plot_regioner(dat$pan_pan_te)
gro_all_te <- plot_regioner(dat$gro_all_te)

# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=7))
)

dir.create("results/figures2/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(pan_all_te, x = 0.5, y=0.5, width = 3.5,height = 2.5)
plotText("A", x = 0.5, y=0.5)

plotGG(pan_pan_te, x = 4.5, y=0.5, width = 3.5,height = 2.5)
plotText("B", x = 4.5, y=.5)

plotGG(gro_all_te, x = 0.5, y=3.25, width = 3.5,height = 2.5)
plotText("C", x = 0.5, y=3.25)

dev.off()