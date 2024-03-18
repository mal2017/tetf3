library(tidyverse)
library(plotgardener)
library(patchwork)

motif_fig_df <- read_rds("results/motifs/comparison/pan_denovo_comparison.meme.gg_df.rds") |>
  mutate(g_rnk = map2(denovo,g_rnk,~{.y + labs(title=.x)}))

maplot <- read_rds("results/ripseq/unr_ripseq_maplot.gg.rds")

au <- read_rds("results/ripseq/unr_ripseq_au_richness.boxplot.gg_list.rds")$TE

attta <- read_rds("results/ripseq/unr_ripseq_attta_richness.boxplot.gg_list.rds")$TE

# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + theme(text=element_text(size=unit(7,"pt"))))

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())
figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

# ------------------------------------------------------------------------------
# plot motif figs
# ------------------------------------------------------------------------------

g_a <- plotGG(motif_fig_df$g_aln[[1]], x = 0.5, y=0.5, width = 3.5,height = 1.5)
plotText("A", x = 0.5, y=0.5)

g_b <- plotGG(motif_fig_df$g_aln[[2]], x = 4.25, y=0.5, width = 3.5,height = 1.5)
plotText("B", x = 3.5, y=.5)

g_c <- plotGG(Reduce(`+`,motif_fig_df$g_rnk) + plot_layout(nrow=1,guides = "collect") & theme(legend.position = "right") & aes(color=class),
       x=0.75,y=2.25,width = 7,height=2)
plotText("C", x = 0.75, y=2.25)

g_d <- plotGG(maplot,x=0.5,y=4.5,width=2,height=1.5)
plotText("D", x = 0.5, y=4.5)
g_e <- plotGG(au,x=2.75,y=4.5,width=2.5,height=2.125)
plotText("E", x = 2.75, y=4.5)
g_f <- plotGG(attta,x=5.5,y=4.5,width=2.5,height=2.125)
plotText("F", x = 5.5, y=4.5)

dev.off()