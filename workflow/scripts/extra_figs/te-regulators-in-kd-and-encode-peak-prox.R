library(patchwork)
library(plotgardener)
library(tidyverse)
library(clusterProfiler)
library(ggpubr)
library(ggdensity)
library(plyranges)

# ------------------------------------------------------------------------------
# te regulators leading edge bar chart
# ------------------------------------------------------------------------------
x <- read_rds("results/signatures/ourKD_gsea.rds")

x <- x |> mutate(lab = str_replace_all(str_extract(comparison,"(?<=knockdown2_).+(?=_control)"), "_", " / "))

plot_barchart <- function(dat, cmp) {
  dat |>
    filter(signature_name == cmp) |>
    mutate(lab = fct_reorder(lab, pvalue)) |>
    ggplot(aes(-log10(pvalue), lab)) +
    geom_col() +
    geom_vline(xintercept = -log10(0.05),color="red",linetype="dashed") +
    ylab("RNAi / sex / sample / driver")
}

g_a <- plot_barchart(x, "TE.regulators")


# ------------------------------------------------------------------------------
# CHIP
# ------------------------------------------------------------------------------

# plot distance from each gene (incl piRNA pathway genes) to
# NfI, CG16779, and pan peaks from encode
gr <- read_rds("results/pirna/encode_peaks_dist_to_pirna.gr.rds")

ge<- gr |>
  as_tibble() |>
  filter(embryo.expressed & ChIP!="gro") |>
  mutate(type = if_else(is.piRNA.pathway,"TE regulators","other genes")) |>
  mutate(ChIP = paste(ChIP, "peaks")) |>
  ggplot(aes(type,log10(distance+1))) +
  geom_boxplot() +
  facet_wrap(~ChIP) +
  ylab("log10(dist. to nearest peak + 1)") +
  xlab("genes") +
  stat_compare_means(size=unit(2,"pt"),method = 'wilcox.test')

# ------------------------------------------------------------------------------
# te regulators random walk
# ------------------------------------------------------------------------------

# access data with example gbs$`CG16779 female gonad tj vs female gonad tj`[1]
gbcd_df <- read_rds("results/signatures/ourKD_gsea.rds") |>
  filter(signature_name == "TE.regulators") |>
  arrange(pvalue) |>
  mutate(padj = p.adjust(pvalue, method="BH")) |>
  filter(str_detect(comparison,"head")) |>
  filter(kd %in% c("NfI","CG16779","pan","Unr","vvl")) |>
  group_by(kd) |>
  slice_min(pvalue) |>
  ungroup() |>
  mutate(comparison = str_remove(comparison,"knockdown2_")) |>
  mutate(comparison = str_replace(comparison, "_control","_vs")) |>
  mutate(comparison = str_replace_all(comparison,"_"," "))

gbcd_p_lookup <- gbcd_df |>
  dplyr::select(comparison, padj) |>
  deframe() |>
  as.list() |>
  map(format.pval,3) |>
  map(~paste0("padj<",.x))

gbcd <- gbcd_df |>
  dplyr::select(comparison, gsea, padj) |>
  mutate(plot = pmap(list(x=comparison,y=gsea, z=padj),
       function(x, y, z) {
         enrichplot::gseaplot2(y, geneSetID="TE.regulators",
                               title=x, subplots=1:2, base_size = 7,
                               color=if_else(z < 0.1,"green","red"))
         })) |>
  dplyr::select(comparison,plot) |>
  deframe()


# ------------------------------------------------------------------------------
# make page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=unit(7,"pt"))) +
            theme(plot.title = element_text(hjust = 0.5))
)

dir.create("results/figures2/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())



plotGG(g_a, x = 0.5, y=0.5, width = 3.5,height = 2.25)
plotText("A", x = 0.5, y=0.5)


plotGG(gbcd$`NfI female head Mef2.R vs female head Mef2.R`, x = 4.5, y=0.5, width = 3.5,height = 2.25)
plotText(gbcd_p_lookup$`NfI female head Mef2.R vs female head Mef2.R`, x=6.75, y=1,fontsize = 7)
plotText("B", x = 4.5, y=0.5)


plotGG(gbcd$`CG16779 female head Mef2.R vs female head Mef2.R`, x = 0.5, y=3, width = 3.5,height = 2.25)
plotText(gbcd_p_lookup$`CG16779 female head Mef2.R vs female head Mef2.R`, x=2.75, y=3.5,fontsize = 7)
plotText("C", x = .5, y=3)


plotGG(gbcd$`pan female head Mef2.R vs female head Mef2.R`, x = 4.5, y=3, width = 3.5,height = 2.25)
plotText(gbcd_p_lookup$`pan female head Mef2.R vs female head Mef2.R`, x=6.75, y=3.5,fontsize = 7)
plotText("D", x = 4.5, y=3)

plotGG(ge, x = 0.5, y=5.5, width = 3.75,height = 2.5)
plotText("E", x = 0.5, y=5.5)

dev.off()
