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
library(gridExtra)
source("workflow/scripts/utils/plotting.R")


# ------------------------------------------------------------------------------
# meme result table
# ------------------------------------------------------------------------------

fdr_tbl <- read_tsv("results/motifs/streme_per_tf_empirical_fdr/pan_empirical_fdr.tsv") |>
  mutate(`motif sequence` = str_extract(motif_ID,"(?<=\\d{1,2}-).+"))

streme_res <- memes::importStremeXML("results/motifs/streme_per_tf/pan/streme.xml") |>
  as_tibble() |>
  dplyr::select(name=altname, motif,`motif sequence`=name, `p-value`=train_pval) |>
  mutate(`p-value`=as.numeric(`p-value`)) |>
  mutate(`motif sequence` = str_extract(`motif sequence`, "(?<=_).+$"))

g_table <- dplyr::select(streme_res,-motif) |>
  left_join(fdr_tbl, by = c("motif sequence",'p-value'="motif_Score")) |>
  dplyr::select(name, `motif sequence`, `p-value`, `empirical FDR`=fdr)


g_table <- tableGrob(g_table) |> ggplotify::as.ggplot()



# ------------------------------------------------------------------------------
## motif alignments
# ------------------------------------------------------------------------------
motif_comp <- read_rds("results/motifs/comparison/pan_denovo_comparison.streme.rds")

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
  slice_min(padj,n = 5) |>
  filter(padj < 0.1 & (name == "pan" | str_detect(known,"Archbold"))) |>
  pull(denovo)

g_a <- motif_comp2plot |>
  filter(denovo %in% which.motifs) |>
  ggplot(aes(rank,-log10(Pval), color=class)) +
  geom_point(size=rel(0.75)) +
  ggrepel::geom_text_repel(data= \(x) filter(x,class!="jaspar (other)" & rank < 5), aes(label=label)) +
  scale_color_brewer(type = "qual", palette = 2) +
  geom_hline(yintercept = -log10(max_p), linetype="dashed", color="darkgray") +
  facet_wrap(~denovo, scales = "free", nrow=1) +
  theme(legend.position = "bottom")


motif_alns_to_plot <- motif_comp2plot |> filter(class!="jaspar (other)" & padj < 0.1)

g_b <- filter(motif_comp, padj <0.1 & (map2_lgl(name,known, {~str_detect(.y, paste0("^",.x,"$"))}) | name == "pan" | str_detect(known,"degenerate"))) |>
  group_by(denovo) |>
  slice_min(Pval,with_ties = F) |> 
  filter(denovo %in% motif_alns_to_plot$denovo) |>
  pull(gg) |>
  Reduce(`+`,x=_ ) & theme_bw() & 
  guides(color="none", fill="none") & plot_layout(nrow=1) & 
  theme(text = element_text(size=5), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank())



# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=5))
)

dir.create("results/figures/")

pdf("results/figures/denovo-motifs-on-tes-streme-supplement-02.pdf",width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_table, x=1.75, y=0.5, width = 5, height=2)
plotText("A", x = 1.75, y=0.5)

plotGG(g_a, x = 0.5, y=2.75, width = 7,height = 3)
plotText("B", x = 0.5, y=2.5)


plotGG(g_b, x = 0.5, y=6, width = 7,height = 1.5)
plotText("C", x = 0.5, y=6)

dev.off()
