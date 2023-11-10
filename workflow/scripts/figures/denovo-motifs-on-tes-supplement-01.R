library(tidytree)
library(ggtree)
library(tidyverse)
library(memes)
library(universalmotif)
library(gridExtra)
library(plotgardener)
library(patchwork)


# ------------------------------------------------------------------------------
# meme result table
# ------------------------------------------------------------------------------

meme_res <- memes::importMeme("results/motifs/meme_per_tf/pan/meme.txt") |>
  as_tibble() |>
  filter(eval < 0.05) |>
  dplyr::select(name=altname, motif,`motif sequence`=name, `E-value`=eval)

g_table <- dplyr::select(meme_res,-motif)


g_table <- tableGrob(g_table) |> ggplotify::as.ggplot()



# ------------------------------------------------------------------------------
## tree with known and denovo motif rings
# ------------------------------------------------------------------------------
tree_pan_motif_fig <- "results/integrative/motif_and_coex_on_tree.pan.plot.rds"
tree_pan_motif_fig <- read_rds(tree_pan_motif_fig)


g_tree <- tree_pan_motif_fig + 
  geom_tiplab(offset = 17, size=unit(2,"pt")) +
  guides(color = guide_legend(override.aes = list(size = 0.3))) +
  theme(text=element_text(size=7))

g_motif <- universalmotif::view_motifs(meme_res$motif[[1]]) + labs(title = paste0("MEME E-value=",meme_res$`E-value`[[1]])) + theme(plot.title = element_text(hjust=1))

# ------------------------------------------------------------------------------
## denovo motifs coex probability
# ------------------------------------------------------------------------------
#https://stats.oarc.ucla.edu/r/dae/logit-regression/

prob_df <- read_rds("results/integrative/n_denovo_vs_sig_coef.pan.rds")

g_d <- prob_df |>
  filter(motif == "meme_6") |>
  mutate(gg = map2(data, smoothed_data, ~{
    ggplot(data=.x, aes(x=n)) +
      geom_jitter(width=0.1, aes(y=as.integer(is_coex)),height = 0.025) +
      geom_path(data=.y, aes(y=pred)) +
      xlab("N de novo motifs") +
      ylab("Probability of coexpression")
  })) |>
  pull(gg) |>
  pluck(1)
  
g_wilc <- prob_df |>
    filter(motif == "meme_6") |>
    pull(data) |>
  pluck(1) |>
  mutate(`pan coexpressed`=if_else(is_coex,"coexpressed","n.s.")) |>
  ggplot(aes(`pan coexpressed`, n)) +
  geom_boxplot() +
  ggpubr::stat_compare_means()


g_hist <- prob_df |>
  filter(motif == "meme_6") |>
  pull(data) |>
  pluck(1) |>
  mutate(`pan coexpressed`=if_else(is_coex,"coexpressed","n.s.")) |>
  ggplot(aes(n)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~`pan coexpressed`, ncol=1, scales="free_y")


# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=5))
)

dir.create("results/figures/")

pdf("results/figures/denovo-motifs-on-tes-supplement-01.pdf",width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_table, x=1.75, y=0, width = 5, height=2)
plotText("A", x = 2, y=0.5)

plotGG(g_tree, x = 0.75, y=2, width = 2*3.5,height = 2*2.5)
plotText("B", x = 1.5, y=2.5)

plotGG(g_wilc, x=1, y=7.25, width = 2, height = 2)
plotText("C", x = 1, y=7.125)

plotGG(g_d, x=3.25, y=7.25, width = 2, height = 2)
plotText("D", x = 3.25, y=7.125)

plotGG(g_hist, x=5.5, y=7.25, width = 2, height = 2)
plotText("E", x = 5.5, y=7.125)

dev.off()