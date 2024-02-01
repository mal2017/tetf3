library(tidyverse)
library(plotgardener)
library(patchwork)
library(GenomicFeatures)
library(plyranges)
library(gridExtra)
library(grid)

regioner_df <- read_rds("results/csem_mosaics/regioner.rds")

# function for nicely plotting randomization results
plot_regioner <- function(x,g,m,ts) {
  x <- filter(x,ChIP==g & masking==m & te_set == ts)$regioneR_results
  stopifnot(nrow(x)==1)
  x <- x[[1]]
  
  restricted_to <- if_else(m=="none",NA,"considered only euchromatic peaks/TEs and euchromatic regions")
  lb <- sprintf("%s CSEM/MOSAICS peaks; considering %s TEs", g, str_replace(ts,"_","-"))
  
  permuted_overlaps <- tibble(overlaps = x$numOverlaps$permuted)
  observed_overlaps <- tibble(overlaps = x$numOverlaps$observed)
  pval <- x$numOverlaps$pval
  z <- x$numOverlaps$zscore
  alt <- x$numOverlaps$alternative
  
  g2 <- ggplot(permuted_overlaps, aes(overlaps)) +
    geom_histogram() +
    geom_vline(data = observed_overlaps, aes(xintercept = overlaps), color="darkgreen") +
    annotate("text", x=-Inf, y=Inf, label=sprintf("z=%s, p<%s, alternative=%s", format(z, digits=3), format.pval(pval,digits = 2), alt), hjust=0, vjust=1) +
    labs(title=lb) +
    ylab("n")
  
  if (!is.na(restricted_to)) {
    g2 <- g2 + labs(subtitle = restricted_to)
  }
  
  return(g2)
}


h3k9me3_all_te <- plot_regioner(regioner_df, "H3K9Me3","none","all")
gro_all_te <- plot_regioner(regioner_df, "gro","none","all")
pan_all_te <- plot_regioner(regioner_df, "pan","none","all")
pan_pan_te <- plot_regioner(regioner_df, "pan","none","factor_specific")
pan_euch_all_te <- plot_regioner(regioner_df, "pan","mask_heterochromatin","all")
nfi_all_te <- plot_regioner(regioner_df, "NfI","none","all")
vvl_all_te <- plot_regioner(regioner_df, "vvl","none","all")
cg16779_all_te <- plot_regioner(regioner_df, "CG16779","none","all")



  

# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=7))
)

dir.create("results/figures2/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(h3k9me3_all_te, x = 0.5, y=0.5, width = 3.5,height = 2.5)
plotText("A", x = 0.5, y=0.5)

plotGG(gro_all_te, x = 4.5, y=0.5, width = 3.5,height = 2.5)
plotText("B", x = 4.5, y=.5)

plotGG(pan_pan_te, x = 0.5, y=3.25, width = 3.5,height = 2.5)
plotText("C", x = 0.5, y=3.25)

plotGG(pan_all_te, x = 4.5, y=3.25, width = 3.5,height = 2.5)
plotText("D", x = 4.5, y=3.25)


plotGG(pan_euch_all_te, x = 0.5, y=6, width = 3.5,height = 2.5)
plotText("E", x = 0.5, y=6)

plotGG(vvl_all_te, x = 4.5, y=6, width = 3.5,height = 2.5)
plotText("F", x = 4.5, y=6)

plotGG(nfi_all_te, x = 0.5, y=8.75, width = 3.5,height = 2.5)
plotText("G", x = 0.5, y=8.75)

plotGG(cg16779_all_te, x = 4.5, y=8.75, width = 3.5,height = 2.5)
plotText("H", x = 4.5, y=8.75)

dev.off()