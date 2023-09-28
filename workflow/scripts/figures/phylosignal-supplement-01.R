library(patchwork)
library(plotgardener)
library(tidyverse)
library(clusterProfiler)
library(ggpubr)
library(phylosignal)
library(phylobase)
library(ggplotify)
library(ggtree)
library(tidytree)
library(ggtreeExtra)

# ------------------------------------------------------------------------------
# disttribution of distances
# ------------------------------------------------------------------------------
g_a <- read_rds("results/te_sequence_similarity/coex_vs_seq_similarity.rds") |>
  filter(model == "male") |>
  dplyr::select(te1, te2, seq_distance) |>
  ggplot(aes(seq_distance)) +
  geom_histogram() +
  ylab("TE/TE pairs") + xlab("dashing2 distance")

# ------------------------------------------------------------------------------
# phylosignal tf results by test
# ------------------------------------------------------------------------------
zads <- read_tsv("results/resources/zad_genes.tsv")

ps_df <- read_rds("results/phylosignal/phylosignal_df.rds")

n_tfs <- ps_df |>
  filter(score_type == "score") |>
  pull(TF) |>
  unique() |> length()

g_b <- ps_df |>
  filter(score_type =="control") |>
  ggplot(aes(-log10(pval),coef,fill=metric)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(type = "qual", palette = 2) +
  geom_vline(xintercept = -log10(0.05), linetype="dashed",color="darkgray") +
  ylab("control")

# ------------------------------------------------------------------------------
# control correlograms
# ------------------------------------------------------------------------------
crlgs <- read_rds("results/phylosignal/correlograms.rds")

get_col <- \(x) c("red", "orange", "darkblue")[match(x,c("bm","ou","random"))]

plot_crlg <-   function(dat,lab) {
    h0 <- -1/(dat$n - 1)
    as_tibble(dat$res) |>
      set_names(c("x","ci.upper","ci.lower","y")) |>
      ggplot(aes(x)) +
        geom_path(aes(y=y),color=get_col(lab)) +
        geom_path(aes(y=ci.upper), linetype="dotted") +
        geom_path(aes(y=ci.lower), linetype="dotted") +
        xlab("patristic distance") + ylab("correlation") +
        labs(title=lab) +
      geom_hline(yintercept = h0, color="gray")
}

g_crlgs <- crlgs |>
  imap(plot_crlg)

# convert to ggplot/grid
#gpfy <- function(p) {
#  as.ggplot(function() plot(p))
#}

g_c <- g_crlgs$random
g_d <- g_crlgs$bm
g_e <- g_crlgs$ou

# ------------------------------------------------------------------------------
# phylogenies + barplot of control values
# ------------------------------------------------------------------------------
ps <- read_rds("results/phylosignal/phylosignal.rds")

#g_c <- ggplotify::as.ggplot(~barplot.phylo4d(ps$p4d,trait = c("bm","ou","random"),tree.ladderize = T, tip.cex = 0.6))

tr <- read_rds("results/te_sequence_similarity/te_sketch_tidytree.rds")

tr_dat <- tibble(label = ps$p4d@label,
                 bm=ps$p4d@data$bm,
                 ou=ps$p4d@data$ou,
                 random=ps$p4d@data$random)

#tr <- tr |> left_join(tr_dat, by="label")

p <- ggtree(tr, layout="rectangular") +
  geom_tiplab(size=5, hjust = 0, linetype = NA,as_ylab = T)

# from https://rdrr.io/bioc/ggtreeExtra/man/geom_fruit.html:
# "...the column contained tree tip labels should be as y in mapping."
# geom_fruit can have an axis: see http://yulab-smu.top/treedata-book/chapter10.html#geom-fruit1
g_f <- p +
  #geom_facet(geom=geom_text, panel="TE", data=tr_dat, aes(x=1, label=label), orientation="y") +
  #geom_facet(geom=geom_col, panel="random", data=tr_dat, aes(x=random), orientation="y", fill="darkblue") +
  #geom_facet(geom=geom_col, panel="bm", data=tr_dat, aes(x=bm), orientation="y", fill="red") +
  #geom_facet(geom=geom_col, panel="ou", data=tr_dat, aes(x=ou), orientation="y", fill="orange") +
  geom_fruit(data=tr_dat, geom=geom_col, aes(y=label,x=random),offset = 0.15, pwidth = 0.5, fill="darkblue",axis.params= list(axis = "x", text.size=2, hjust = 1, vjust = 0.5, nbreak=1, limits=c(-0.2,0.2))) +
  geom_fruit(data=tr_dat,geom=geom_col, aes(y=label,x=bm),offset = 0.5, pwidth = 0.5, fill="red",axis.params = list(axis = "x",  hjust = 1, text.size=2, vjust = 0.5, nbreak = 1, limits=c(-0.2,0.2))) +
  geom_fruit(data=tr_dat,geom=geom_col, aes(y=label,x=ou),offset = 0.5,pwidth = 0.5, fill="orange",axis.params = list(axis = "x", hjust = 1, text.size=2, vjust = 0.5, nbreak = 1, limits=c(-0.2,0.2)))

# ------------------------------------------------------------------------------
# make page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=unit(7,"pt"))) +
            theme(plot.title = element_text(hjust = 0.5))
)

dir.create("results/figures/")

pdf("results/figures/phylosignal-supplement-01.pdf",width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_a, x = 0.5, y=0.5, width = 3,height = 2.75*5/9)
plotText("A", x = 0.5, y=0.5)


plotGG(g_b, x = 4, y=0.5, width = 3.5,height = 2.75*5/9)
plotText("B", x = 4, y=0.5)

# ---------
plotGG(g_c, x = 0.625, y=2, width = 2.25,height = 2)
plotText("C", x = 0.5, y=2.125)


plotGG(g_d, x=3, y = 2, width = 2.25, height = 2)
plotText("D", x = 3, y=2.125)


plotGG(g_e, x=5.5, y = 2, width = 2.25, height = 2)
plotText("E", x = 5.5, y=2.125)

# ---------

plotGG(g_f, x=0.5, y = 4.25-0.25, width = 7.5, height = 6.75+.25)
plotText("F", x = 0.5, y=4.25-0.125)

dev.off()

