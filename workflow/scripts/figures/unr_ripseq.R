Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)
library(patchwork)
library(plotgardener)

relpos <- 50
relpos <- snakemake@params$relpos

# ------------------------------------------------------------------------------
# leading edge in knockdown
# ------------------------------------------------------------------------------

gsea <- read_rds("results/ripseq/unr_bound_tx_in_kd.gsea.rds")

g_bound_te_gsea <- gsea |> enrichplot::gseaplot2("bound_TE") & 
  theme(axis.title = element_text(size=5),axis.text = element_text(size=5))
g_bound_gene_gsea <- gsea |> enrichplot::gseaplot2("bound_gene") &
  theme(axis.title = element_text(size=5),axis.text = element_text(size=5))

# ------------------------------------------------------------------------------
# au richness
# ------------------------------------------------------------------------------
g_au_line_gene <- read_rds("results/ripseq/unr_ripseq_au_richness.lineplot.gg_list.rds")$gene
g_au_line_te <- read_rds("results/ripseq/unr_ripseq_au_richness.lineplot.gg_list.rds")$TE

g_au_box_gene <- read_rds("results/ripseq/unr_ripseq_au_richness.boxplot.gg_list.rds")$gene

# ------------------------------------------------------------------------------
# attta sites
# ------------------------------------------------------------------------------

g_attta_box_gene <- read_rds("results/ripseq/unr_ripseq_attta_richness.boxplot.gg_list.rds")$gene

# ------------------------------------------------------------------------------
# create page 1
# ------------------------------------------------------------------------------

g_a <- g_bound_gene_gsea + xlab("rank in Unr knockdown")
g_b <- g_bound_te_gsea + xlab("rank in Unr knockdown")
g_c <- g_au_line_gene
g_d <- g_au_line_te
g_e <- g_au_box_gene
g_f <- g_attta_box_gene

theme_set(theme_classic() + 
            theme(text = element_text(size=7)) +
            theme(plot.title = element_text(hjust = 0.5))
          )

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())
figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

plotGG(g_a, 
       x = 0.75, y=0.5, width = 3.5,height = 2.25)
plotText("Unr-bound genes",x=2.6,y=0.75,fontsize = 7)
plotText("A", x = 0.75, y=0.5)

plotGG(g_b, 
       x = 4.5, y=0.5, width = 3.5,height = 2.25)
plotText("Unr-bound TEs",x=6.5,y=0.75,fontsize = 7)
plotText("B", x = 4.5, y=0.5)

plotGG(g_c, x = 1, y=3, width = 6.5,height = 2)
plotText("C",  x = 1, y=3)

plotGG(g_d, x = 1, y=5.25, width = 6.5,height = 2)
plotText("D",  x = 1, y=5.25)


plotGG(g_e,
       x = 1.75, y=7.5, width = 2.25,height = 2.25)
plotText("E",  x = 1.65, y=7.5)

plotGG(g_f,
       x = 4.75, y=7.5, width = 2.25,height = 2.25)
plotText("F",  x = 4.75, y=7.5)


dev.off()

writexl::write_xlsx(list(A=g_a$data,
                         B=g_b$data,
                         C=g_c$data,
                         D=g_d$data,
                         E=g_e$data,
                         `F`=g_f$data),
                    path = ifelse(exists("snakemake"),
                                  snakemake@output$xlsx,
                                  "~/Downloads/test.xlsx"))
