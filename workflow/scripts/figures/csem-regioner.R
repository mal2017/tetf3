Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


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

dir.create("results/figures/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())
figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

g_a <- h3k9me3_all_te
plotGG(g_a, x = 0.5, y=0.5, width = 3.5,height = 2.25)
plotText("A", x = 0.5, y=0.5)

g_b <- gro_all_te
plotGG(g_b, x = 4.5, y=0.5, width = 3.5,height = 2.25)
plotText("B", x = 4.5, y=.5)

g_c <- pan_pan_te
plotGG(g_c, x = 0.5, y=3, width = 3.5,height = 2.25)
plotText("C", x = 0.5, y=3)

g_d <- pan_all_te
plotGG(g_d, x = 4.5, y=3, width = 3.5,height = 2.25)
plotText("D", x = 4.5, y=3)

g_e <- pan_euch_all_te
plotGG(g_e, x = 0.5, y=5.75, width = 3.5,height = 2.25)
plotText("E", x = 0.5, y=5.5)

g_f <- vvl_all_te
plotGG(g_f, x = 4.5, y=5.75, width = 3.5,height = 2.25)
plotText("F", x = 4.5, y=5.5)

g_g <- nfi_all_te
plotGG(g_g, x = 0.5, y=8, width = 3.5,height = 2.25)
plotText("G", x = 0.5, y=8)

g_h <- cg16779_all_te
plotGG(g_h, x = 4.5, y=8, width = 3.5,height = 2.25)
plotText("H", x = 4.5, y=8)

dev.off()

writexl::write_xlsx(list(A=g_a$data,
                         B=g_b$data,
                         C=g_c$data,
                         D=g_d$data,
                         E=g_e$data,
                         `F`=g_f$data,
                         G=g_g$data,
                         H=g_h$data),
                    path = ifelse(exists("snakemake"),
                                  snakemake@output$xlsx,
                                  "~/Downloads/test.xlsx"))