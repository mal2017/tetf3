library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(ggtext)
library(patchwork)
library(plotgardener)
library(tidyverse)

# ------------------------------------------------------------------------------
# get signature enrichment results from our kd
# ------------------------------------------------------------------------------
x <- read_rds("results/signatures/ourKD_gsea.rds")

x <- x |> mutate(lab = str_replace_all(str_extract(comparison,"(?<=knockdown2_).+(?=_control)"), "_", " / "))
  

gg <- x |> 
  filter(padj < 0.1 & signature_name == "all_tes") |>
  group_by(signature_name, kd) |>
  slice_min(pvalue) |>
  ungroup() |>
  mutate(gg = pmap(list(lab, ID, gsea),
                   .f = function(lab,gs, obj) {
                     enrichplot::gseaplot2(obj, geneSetID=gs, title = lab)
                   }
  )) |>
  dplyr::select(comparison, signature_name, gg)

# ------------------------------------------------------------------------------
# target-specific all-te hit barchart
# ------------------------------------------------------------------------------

g_a <- x |>
  filter(signature_name == "all_tes") |>
  mutate(lab = fct_reorder(lab, pvalue)) |>
  ggplot(aes(-log10(padj), lab)) +
  geom_col() +
  geom_vline(xintercept = -log10(0.1),color="red",linetype="dashed") +
  ylab("RNAi / sex / sample / driver")


# ------------------------------------------------------------------------------
# exemplary all-te random walk
# ------------------------------------------------------------------------------

g_bcd <- gg |>
  dplyr::select(comparison, gg) |>
  deframe() |>
  map( ~{ .x & theme(axis.title = element_text(size=5), axis.text = element_text(size=5), plot.title = element_text(size=7, hjust=0.5))})

# ------------------------------------------------------------------------------
# create page
  # ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=7))
)

dir.create("results/figures2/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_a, x = 0.5, y=0.5, width = 3.75,height = 2.5)
plotText("A", x = 0.5, y=0.5)

plotGG(g_bcd$knockdown2_NFI_female_head_Mef2.R_control_female_head_Mef2.R, 
       x = 4.5, y=0.5, width = 3.25,height = 2.5)
plotText("B", x = 4.5, y=0.5)

plotGG(g_bcd$knockdown2_pan_female_gonad_tj_control_female_gonad_tj, x = 0.75, y=3.5, width = 3.25,height = 2.5)
plotText("C",  x = 0.5, y=3.5)

plotGG(g_bcd$knockdown2_CG16779_male_gonad_aTub_control_male_gonad_aTub, x = 4.5, y=3.5, width = 3.25,height = 2.5)
plotText("D",  x = 4.5, y=3.5)

dev.off()

