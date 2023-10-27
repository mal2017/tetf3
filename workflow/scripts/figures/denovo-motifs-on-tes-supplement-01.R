library(tidytree)
library(ggtree)

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

