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

# ------------------------------------------------------------------------------
## motif alignments
# ------------------------------------------------------------------------------
motif_comp <- read_rds("results/motifs/comparison/pan_denovo_comparison.rds")

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

g_a <- motif_comp2plot |>
  ggplot(aes(rank,-log10(Pval), color=class)) +
  geom_point(size=rel(0.75)) +
  ggrepel::geom_text_repel(data= \(x) filter(x,class!="jaspar (other)" & rank < 5), aes(label=label)) +
  scale_color_brewer(type = "qual", palette = 2) +
  geom_hline(yintercept = -log10(max_p), linetype="dashed", color="darkgray") +
  facet_wrap(~denovo, scales = "free", nrow=1) +
  theme(legend.position = "bottom")


motif_alns_to_plot <- motif_comp2plot |> filter(class!="jaspar (other)" & rank < 5)

g_b <- filter(motif_comp, padj <0.1 & (map2_lgl(name,known, {~str_detect(.y, paste0("^",.x,"$"))}) | name == "pan" | str_detect(known,"degenerate"))) |>
  group_by(denovo) |>
  slice_min(Pval,with_ties = F) |> 
  filter(denovo %in% motif_alns_to_plot$denovo) |>
  pull(gg) |>
  Reduce(`+`,x=_ ) & theme_bw() & 
  guides(color="none", fill="none") & plot_layout(nrow=1) & 
  theme(text = element_text(size=5), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank())

# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

sl <- getChromInfoFromNCBI("GCF_000001215.4")

chroi <-  c("2L","2R","3L","3R","X","4")

# get full range of main chroms
wh <-  sl |>
  as_tibble() |>
  dplyr::filter(SequenceName %in%chroi) |>
  dplyr::select(chr=SequenceName,end=SequenceLength) |>
  mutate(start=1) |>
  GRanges()

# get names/paths of bigwigs

bws <- Sys.glob("~/amarel-matt/tetf/subworkflows/tetf_csem_mosaics/results/csem_mosaics/viz/pan_ENCSR058DSI_rep*.log2ratio.bw") |>
  c(Sys.glob("~/amarel-matt/tetf/subworkflows/tetf_csem_mosaics/results/csem_mosaics/viz/gro_ENCSR981DLO_rep1.log2ratio.bw")) |>
  c(Sys.glob("~/amarel-matt/tetf/subworkflows/tetf_csem_mosaics/results/csem_mosaics/viz/E0.4_H3K9Me3_ChIPSeq_1.log2ratio.bw"))

names(bws) <- str_extract(bws,"(?<=viz\\/).+(?=\\.log2)")

# get tiles and average within each tile - plotting takes forever with default 50bp windows
tiles <- tileGenome(deframe(sl[,c("SequenceName","SequenceLength")])[chroi],tilewidth = 50000, cut.last.tile.in.chrom=TRUE)

grs <- map(bws, import_and_summarize, tiles, wh) # this func from workflow/scripts/utils/plotting.R

gs <- grs |> GRangesList() |> plot_genome_signal()

g_tracks <- gs@ggplot 

# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=5))
)

dir.create("results/figures/")

pdf("results/figures/figure5.pdf",width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_a, x = 0.5, y=0.5, width = 7,height = 2.1)
plotText("A", x = 0.5, y=0.5)


plotGG(g_b, x = 0.5, y=2.5, width = 7,height = 1.5)
plotText("B", x = 0.5, y=2.5)

plotGG(g_tracks, x=0.5, y=4, width = 7.5, height = 2)
plotText("C", x = 0.5, y=4.1)


dev.off()