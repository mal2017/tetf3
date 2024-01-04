library(patchwork)
library(plotgardener)
library(tidyverse)
library(clusterProfiler)
library(ggpubr)
library(ggdensity)
library(plyranges)

# plot distance from each gene (incl piRNA pathway genes) to
# NfI, CG16779, and pan peaks from REMAP22
gr <- read_rds("results/pirna/remap_peaks_dist_to_pirna.gr.rds")

gd<- gr |>
  as_tibble() |>
  mutate(type = if_else(is.piRNA.pathway,"TE regulators","other genes")) |>
  mutate(ChIP = paste(ChIP, "peaks")) |>
  ggplot(aes(type,log10(distance+1))) +
  geom_violin() +
  facet_wrap(~ChIP) +
  ylab("log10(dist. to nearest peak + 1)") +
  xlab("genes") +
  stat_compare_means(size=unit(2,"pt"))


# access data with example gbs$`CG16779 female gonad tj vs female gonad tj`[1]
gabc_df <- read_rds("results/signatures/ourKD_gsea.rds") |>
  filter(signature_name == "TE.regulators") |>
  arrange(pvalue) |>
  mutate(padj = p.adjust(pvalue, method="BH")) |>
  filter(kd %in% c("NfI","CG16779","pan")) |>
  group_by(kd) |>
  slice_min(pvalue) |>
  ungroup() |>
  mutate(comparison = str_remove(comparison,"knockdown2_")) |>
  mutate(comparison = str_replace(comparison, "_control","_vs")) |>
  mutate(comparison = str_replace_all(comparison,"_"," "))

gabc_p_lookup <- gabc_df |>
  dplyr::select(comparison, padj) |>
  deframe() |>
  as.list() |>
  map(format.pval,3) |>
  map(~paste0("padj<",.x))

gabc <- gabc_df |>
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

dir.create("results/figures/")

pdf("results/figures/candidate_tfs_and_piRNA_pathway_supplement.pdf",width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(gabc$`NFI female head Mef2.R vs female head Mef2.R`, x = 0.5, y=0.5, width = 3.5,height = 2.25)
plotText(gabc_p_lookup$`NFI female head Mef2.R vs female head Mef2.R`, x=3, y=1,fontsize = 7)
plotText("A", x = 0.5, y=0.5)


plotGG(gabc$`CG16779 female gonad tj vs female gonad tj`, x = 4.5, y=0.5, width = 3.5,height = 2.25)
plotText(gabc_p_lookup$`CG16779 female gonad tj vs female gonad tj`, x=7, y=1,fontsize = 7)
plotText("B", x = 4.5, y=0.5)


plotGG(gabc$`pan female gonad nos vs female gonad nos`, x = 0.5, y=3, width = 3.5,height = 2.25)
plotText(gabc_p_lookup$`pan female gonad nos vs female gonad nos`, x=3, y=3.5,fontsize = 7)
plotText("C", x = 0.5, y=3)

plotGG(gd, x = 4.25, y=3, width = 3.75,height = 2.5)
plotText("D", x = 4.25, y=3)

dev.off()
