
res.tbl |>
  filter(coef == "zad_mean")


ggplot(res.tbl, aes(stat)) +
  geom_histogram() +
  facet_wrap(~metric)

hits <- res.tbl |> 
  filter(score_type == "score") |> 
  group_by(sex, metric) |>
  mutate(scaled.stat = scale(stat)[,1]) |>
  group_by(sex,TF, coef) |>
  summarise(n_sig_metrics= sum(padj < 0.1), mean_scaled.stat = mean(scaled.stat),.groups = "drop") |>
  filter(n_sig_metrics ==5) |>
  arrange(-abs(mean_scaled.stat)) |>
  mutate(lab = coef)

hits|>
  mutate(is.zad = TF %in% zad$gene_symbol) |>
  #filter(is.zad)
  ggplot(aes(sex, fill=is.zad)) +
  geom_bar(position = "stack") +
  labs(title="significant using 5/5\nphylosignal metrics")

#-------
#x <- x[,zads[zads %in% colnames(x)]]
#x <- x[,colnames(x) %in% head(hits$TE,50)]

hc_tfs <- dist(t(dat)) |> hclust(method="ward.D")


library(tidytree)
library(ggtree)

g <- ggtree(tree) + geom_tiplab(size=2)

lkup <- read_tsv("http://ftp.flybase.net/releases/current/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2023_04.tsv.gz", skip=4) |>
  dplyr::select(gene_ID,gene_symbol) |>
  deframe()

zad_order <- zad |>
  set_names() |>
  map_chr(.f=function(x){
    if (x %in% names(lkup)) {
      lkup[[x]]
    } else {
      x
    }
  })



gheatmap(g, (x[,zad_order[zad_order %in% colnames(x)]]), 
         colnames_position = "bottom", 
         font.size=2, width = 0.99,
         color = NA, colnames_angle = -45,
         hjust =0) +
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 25)) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkgreen",midpoint = 0) +
  #scale_y_continuous(expand=c(0, 0.3))
  theme_tree()



# -----
walk(hits$lab, .f = ~barplot.phylo4d(p4d,trait = .x, center=F, scale=F, main=.x, ladderize=T))

barplot.phylo4d(p4d,trait = c("female_score_CG11902","male_score_CG11902",
                              "female_score_nau","male_score_nau"), center=F, scale=F, main="", ladderize=T)

dat[,c("female_score_CG11902","male_score_CG11902",
       "female_score_nau","male_score_nau")] |>
  ggplot(aes(male_score_nau,male_score_CG11902)) +
  geom_point() +
  ggpubr::stat_cor() +
  geom_smooth(method="lm")


crlgBM <- phyloCorrelogram(p4d, trait = "bm")
crlgRand <- phyloCorrelogram(p4d, trait = "random")
crlgZad <- phyloCorrelogram(p4d, trait = "zad_mean")

par(mfrow=c(1,3))
plot(crlgBM)
plot(crlgRand)
plot(crlgZad)
par(mfrow = c(1, 1))


barplot.phylo4d(p4d,trait = c("random","bm"), center=T, scale=T, main="", ladderize=T)
barplot.phylo4d(p4d,trait = c("zad_mean"), center=T, scale=T, main="", ladderize=T)



