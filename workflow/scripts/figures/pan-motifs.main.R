Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


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

g_a <- motif_fig_df$g_aln[[1]]
g_b <- motif_fig_df$g_aln[[5]]
g_cd <- motif_fig_df$g_rnk[[1]] + motif_fig_df$g_rnk[[5]] + plot_layout(nrow=1,guides = "collect") & theme(legend.position = "right")
g_e <- maplot
g_f <- au
g_g <- attta

theme_set(theme_classic() + theme(text=element_text(size=unit(7,"pt"))))

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())
figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

# ------------------------------------------------------------------------------
# plot motif figs
# ------------------------------------------------------------------------------

plotGG(g_a, x = 0.5, y=0.5, width = 3.5,height = 1.5)
plotText("A", x = 0.5, y=0.5)

plotGG(g_b, x = 4.25, y=0.5, width = 3.5,height = 1.5)
plotText("B", x = 4.25, y=.5)

plotGG(g_cd,
       x=0.75,y=2.25,width = 7,height=2)
plotText("C", x = 0.75, y=2.25)
plotText("D", x = 4.125, y=2.25)

plotGG(g_e,x=0.5,y=4.5,width=2,height=1.5)
plotText("E", x = 0.5, y=4.5)

plotGG(g_f,x=2.75,y=4.5,width=2.5,height=2.125)
plotText("F", x = 2.75, y=4.5)

plotGG(g_g,x=5.5,y=4.5,width=2.5,height=2.125)
plotText("G", x = 5.5, y=4.5)

dev.off()


writexl::write_xlsx(list(A=g_a$data,
                         B=g_b$data,
                         C=motif_fig_df$g_rnk[[1]]$data,
                         D=motif_fig_df$g_rnk[[5]]$data,
                         E=g_e$data,
                         `F`=g_f$data,
                         G=g_g$data),
                    path = ifelse(exists("snakemake"),
                                  snakemake@output$xlsx,
                                  "~/Downloads/test.xlsx"))