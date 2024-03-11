library(tidyverse)
library(tidytree)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(plotgardener)
library(patchwork)
library(GenomeInfoDb)
library(AnnotationDbi)
library(ggbio)
library(org.Dm.eg.db)
library(GenomicFeatures)
library(rtracklayer)
library(plyranges)
library(GenomicFeatures)
source("workflow/scripts/utils/plotting.R")
library(marge)
library(gridExtra)

# ------------------------------------------------------------------------------
## motif alignments
# ------------------------------------------------------------------------------
motif_comp <- read_rds("results/motifs/comparison/pan_denovo_comparison.homer.rds")

motif_comp <- motif_comp |>
  group_by(denovo) |>
  arrange(Pval) |>
  mutate(rank = row_number()) |>
  ungroup()

# highest p at bh adjusted p of 0.1
max_p <- filter(motif_comp, padj < 0.1) |> pull(Pval) |> max()

motif_comp2plot <- motif_comp |>
  dplyr::select(rank, denovo, known, name, Pval, padj) |>
  mutate(class = case_when(str_detect(known,"degenerate")~"Archbold 2014 degenerate",
                           str_detect(known,"::MA0") & name=="pan" ~"pan (jaspar)",
                           str_detect(known,"known::MA")~"jaspar (other)",
                           str_detect(known,"known::")~"Archbold 2014 HMG/helper",
                           T~"wut")) |>
  mutate(label = if_else(is.na(name),class,name))
#filter(class == "pan (jaspar)")

which.motifs <- motif_comp2plot |>
  group_by(denovo) |>
  slice_min(padj) |>
  filter(padj < 0.1 & (name == "pan")) |>
  pull(denovo)

g_a <- g_a <- motif_comp2plot |>
  filter(denovo %in% which.motifs) |>
  ggplot(aes(rank,-log10(Pval), color=class)) +
  geom_point(data=\(x){filter(x, padj >= 0.1)}, size=rel(0.75),color="gray") +
  geom_point(data=\(x){filter(x, padj < 0.1)}, size=rel(1)) +
  ggrepel::geom_text_repel(data= \(x) filter(x,rank <= 1), aes(label=label),color="black") +
  scale_color_brewer(type = "qual", palette = 2,name="") +
  geom_hline(yintercept = -log10(max_p), linetype="dashed", color="darkgray") +
  facet_wrap(~denovo, scales = "free", nrow=1) +
  theme(legend.position = "bottom")

motif_alns_to_plot <- motif_comp2plot |> filter(class!="jaspar (other)" & rank < 10)

g_b <- filter(motif_comp, padj <0.1 & (map2_lgl(name,known, {~str_detect(.y, paste0("^",.x,"$"))}) | name == "pan")) |>
  group_by(denovo) |>
  slice_min(Pval,with_ties = F) |> 
  filter(denovo %in% motif_alns_to_plot$denovo) |>
  pull(gg) |>
  Reduce(`+`,x=_ ) & theme_bw() & 
  guides(color="none", fill="none") & plot_layout(nrow=1) & 
  theme(text = element_text(size=5), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank())

# ------------------------------------------------------------------------------
# motif table
# ------------------------------------------------------------------------------

g_table <- group_by(motif_comp, denovo) |>
  slice_head(n=1) |>
  ungroup(denovo) |>
  dplyr::select(denovo, motifs) |>
  mutate(motifs = map(motifs,pluck,2)) |>
  mutate(`p-value`=map_dbl(motifs,`@`,pval)) |>
  arrange(`p-value`) |>
  mutate(motif_name = paste0("HOMER-",row_number())) |>
  mutate(sequence = str_extract(denovo,"(?<=::).+")) |>
  dplyr::select(motif_name, sequence, `p-value`)

g_table <- tableGrob(g_table) |> ggplotify::as.ggplot()


# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=5))
)

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_table, x=1.75, y=1.25, width = 5, height=2)
plotText("A", x = 2.25, y=0.5)

plotGG(g_a, x = 2.5, y=4.5, width = 3.5,height = 3)
plotText("B", x = 2.25, y=4.5)

plotGG(g_b, x = 2.5, y=7.5, width = 3.5,height = 1.5)
plotText("C", x = 2.25, y=7.5)

dev.off()
