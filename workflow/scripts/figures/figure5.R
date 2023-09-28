library(tidyverse)
library(tidytree)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(plotgardener)
library(patchwork)

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
  
g_a <- motif_comp |>
  dplyr::select(rank, denovo, known, name, Pval) |>
  mutate(class = case_when(str_detect(known,"degenerate")~"Archbold 2014 degenerate",
                           str_detect(known,"::MA0") & name=="pan" ~"pan (jaspar)",
                           str_detect(known,"known::MA")~"jaspar (other)",
                           str_detect(known,"known::")~"Archbold 2014 HMG/helper",
                           T~"wut")) |>
  #filter(class == "pan (jaspar)")
  ggplot(aes(rank,-log10(Pval), color=class)) +
  geom_point() +
  geom_hline(yintercept = -log10(max_p), linetype="dashed", color="darkgray") +
  facet_wrap(~denovo, scales = "free", nrow=1) +
  theme(legend.position = "bottom")


g_b <- filter(motif_comp, padj <0.1 & (map2_lgl(name,known, {~str_detect(.y, paste0("^",.x,"$"))}) | name == "pan" | str_detect(known,"degenerate"))) |>
  group_by(denovo) |>
  slice_min(Pval,with_ties = F) |> pull(gg) |>
  Reduce(`+`,x=_ ) & theme_bw() & 
  guides(color="none", fill="none") & plot_layout(nrow=1) & 
  theme(text = element_text(size=5), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank())



# ------------------------------------------------------------------------------
## tree with known and denovo motif rings
# ------------------------------------------------------------------------------
tree_pan_motif_fig <- "results/integrative/motif_and_coex_on_tree.pan.plot.rds"
tree_pan_motif_fig <- read_rds(tree_pan_motif_fig)


g_c <- tree_pan_motif_fig +
  guides(color = guide_legend(override.aes = list(size = 0.75))) +
  guides(color=guide_legend(nrow=8, byrow=TRUE,keywidth = 0.1, keyheight = 0.1)) +
  guides(fill=guide_legend(nrow=8, byrow=TRUE,keywidth = 0.1, keyheight = 0.1))

# ------------------------------------------------------------------------------
## denovo motifs coex probability
# ------------------------------------------------------------------------------
#https://stats.oarc.ucla.edu/r/dae/logit-regression/

prob_df <- read_rds("results/integrative/n_denovo_vs_sig_coef.pan.rds")

g_d <- prob_df |>
  mutate(gg = map2(data, test_data, ~{
    ggplot(data=.x, aes(x=n)) +
      geom_jitter(width=0.1, aes(y=as.integer(is_coex)),height = 0.025) +
      geom_path(data=.y, aes(y=pred)) +
      xlab("N de novo motifs") +
      ylab("Probability of coexpression")
  })) |>
  pull(gg) |>
  pluck(1)

# ------------------------------------------------------------------------------
## repetitivess
# ------------------------------------------------------------------------------
repet <- "results/repetitiveness/chip_repetitiveness.rds"
repet <- read_rds(repet)

g_e_repetitiveness <- repet |>
  mutate(target = fct_reorder(target,estimate)) |>
  ggplot(aes(target,estimate)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ggpubr::stat_compare_means(method="anova",size=rel(2),label.x.npc = "center",label.y.npc = 0.9) +
  ylab("mapped read ratio:\n(IP TE/IP genomic) / (WCE TE/WCE genomic)")




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


plotGG(g_b, x = 0.25, y=2.5, width = 7.5,height = 1.5)
plotText("B", x = 0.5, y=2.5)

plotGG(g_c + theme(legend.text = element_text(size=5), legend.title = element_text(size=5)), x = 1.25, y=2.75, width = 6,height = 6)
plotText("C", x = 1.25, y=4.25)

plotGG(g_d, x = 0.75, y=7.1, width = 2.5,height = 1.7)
plotText("D", x = 0.5, y=7)

plotGG(g_e_repetitiveness, x = 4, y=7.1, width = 3,height = 1.93)
plotText("E", x = 4, y=7)

dev.off()