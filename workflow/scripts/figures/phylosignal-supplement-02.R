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
# phylosignal tf results by test
# ------------------------------------------------------------------------------
zads <- read_tsv("results/resources/zad_genes.tsv")

ps_df <- read_rds("results/phylosignal/phylosignal_df.rds")

n_tfs <- ps_df |>
  filter(score_type == "score") |>
  pull(TF) |>
  unique() |> length()

g_a <- ps_df |>
  filter(score_type == "score" & padj < 0.1) |>
  ggplot(aes(metric, fill=sex)) +
  geom_bar(position="dodge") +
  scale_fill_grey() +
  ylab("sig. hits") +
  xlab("phylosignal test") +
  annotate("text", label=sprintf("N TFs tested=%s", n_tfs),x=1, y=200, size=rel(2), hjust=0)

cutoff <- 2
hits <- ps_df |>
  filter(score_type == "score") |>
  filter(metric!="Lambda") |>
  group_by(sex,TF) |>
  summarise(n_tests_passed = sum(padj < 0.1),  .groups="drop") |>
  mutate(gene.group = if_else(TF %in% zads$gene_symbol,"ZAD","other")) |>
  pivot_wider(names_from = "sex", values_from = "n_tests_passed") |>
  #filter(male>=2 & female>=2) |> 
  mutate(grp = case_when(female >= cutoff & male < cutoff ~ "female only",
                         male >=cutoff & female < cutoff ~ "male only",
                         female >=cutoff & male >= cutoff ~ "both",
                         T ~ "ns")) |>
  filter(grp!="ns")

g_b <- hits |>
  mutate(grp=fct_infreq(grp)) |>
  ggplot(aes(grp,fill=gene.group)) +
  geom_bar() +
  xlab("sex") + ylab(sprintf(">=%s phylosignal hits", cutoff)) +
  scale_fill_manual(values = c("ZAD"="steelblue", other="tomato"))


# ------------------------------------------------------------------------------
# phylogenies + barplot of control values
# ------------------------------------------------------------------------------
ps <- read_rds("results/phylosignal/phylosignal.rds")

#g_c <- ggplotify::as.ggplot(~barplot.phylo4d(ps$p4d,trait = c("bm","ou","random"),tree.ladderize = T, tip.cex = 0.6))


signal <- ps_df |>
  filter(score_type == "score") |>
  filter(metric!="Lambda") |>
  group_by(sex,TF) |>
  summarise(n_tests_passed = sum(padj < 0.1),  .groups="drop") |>
  filter(n_tests_passed >= cutoff) |>
  mutate(coef_name = sprintf("%s_score_%s", sex, TF)) |>
  mutate(res= map(coef_name, .f=function(x){
    tibble(label=ps$p4d@label, 
           score = ps$p4d@data[[x]])
  }
  )) |>
  unnest(res)

#signal <- hits |> 
#  dplyr::select(TF, grp) |>
#  mutate(grp = str_extract(grp,".+(?=\\s)|both")) |>
#  set_names() |>
#  map_df(.f=function(x){
#    tibble(label=ps$p4d@label, 
#           male_score = ps$p4d@data[[sprintf("male_score_%s",x)]],
#           female_score = ps$p4d@data[[sprintf("female_score_%s",x)]])
#  }, .id="TF")

signal2df <- \(sex2use) {
  signal |>
    filter(sex==sex2use) |>
    dplyr::select(label, TF, score) |>
    pivot_wider(names_from = TF, values_from = score) |>
    column_to_rownames("label")
}

signal_square_list <- c("male","female") |>
  set_names() |>
  map(signal2df)

tr <- read_rds("results/te_sequence_similarity/te_sketch_tidytree.rds")

p <- ggtree(tr, layout="rectangular", branch.length = "none") +
  geom_tiplab(size=2, hjust = 0, linetype = NA,as_ylab = T)


tf_hc_list <- map(signal_square_list, t) |> 
  map(dist) |> 
  map(hclust,method="complete")


g_c <- gheatmap(p, signal_square_list$female[,tf_hc_list$female$order], 
         colnames_position = "bottom", 
         font.size=1.9, width = 9,
         colnames_angle =90,
         hjust =1,color = NA) +
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 25)) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "red",midpoint = 0,name="coex. score") +
  #scale_y_continuous(expand=c(0, 0.3))
  theme_tree() +
  theme(legend.position="top",
        legend.justification="right",
        legend.box.spacing = unit(0,"pt"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,40,-40,-20),
        legend.key.height = unit(0.1,"in"))

g_d <- gheatmap(p, signal_square_list$male[,tf_hc_list$male$order], 
                colnames_position = "bottom", 
                font.size=2, width = 9,
                colnames_angle =90,
                hjust =1,color = NA) +
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 25)) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "red",midpoint = 0,name="coex. score") +
  #scale_y_continuous(expand=c(0, 0.3))
  theme_tree() +
  theme(legend.position="top",
        legend.justification="right",
        legend.box.spacing = unit(0,"pt"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,40,-40,-20),
        legend.key.height = unit(0.1,"in"))

#signal |> pull(TF) |> unique() |> walk(message)
# ------------------------------------------------------------------------------
# make page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=unit(7,"pt"))) +
            theme(plot.title = element_text(hjust = 0.5))
)

dir.create("results/figures/")

pdf("results/figures/phylosignal-supplement-02.pdf",width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())


plotGG(g_c, x = 0, y=2, width = 8.5,height = 4.25)
plotText("C", x = 0.5, y=2.5)

plotGG(g_d, x = 0, y=6.2, width = 8.5,height = 4.25)
plotText("D", x = 0.5, y=6.5)



plotGG(g_a, x = 0.5, y=0.5, width = 4.5,height = 2.75*5/9)
plotText("A", x = 0.5, y=0.5)


plotGG(g_b, x = 5, y=0.5, width = 3,height = 2.75*5/9)
plotText("B", x = 5, y=0.5)





dev.off()
