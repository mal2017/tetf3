
#https://cran.r-project.org/web/packages/phylosignal/vignettes/Demo_plots.htmllibr
library(phylosignal)
library(ape)
library(adephylo)
library(phylobase)
library(phylosignal)
library(tidyverse)

lkup <- read_tsv("http://ftp.flybase.net/releases/FB2021_04/precomputed_files/genes/fbgn_annotation_ID.tsv.gz",skip = 4) %>%
  dplyr::select(gene_id = `primary_FBgn#`,gene_symbol = `##gene_symbol`) %>%
  distinct() %>%
  deframe()

tree <- ape::read.tree("results/zad/iqtree/domain-hits/PF07776/PF07776.contree")

rename_tree_nodes <- function(x) {
  # get a weird warning later if itnernal nodes aren't unique
  # the example in phylosignal doesn't even have itnernal node names
  # so I remove them
  x$node.label <- NULL
  x$tip.label <- x$tip.label %>% str_extract(".+(?=\\.FBpp)") %>% lkup[.]
  return(x)
}

tree <- rename_tree_nodes(tree)

coefs <- read_tsv("~/work/tetf3/upstream/final-models.collected-info2.tsv.gz")

#considered adding neg.pscore=p.value_anova_scaled.copies.y as a negative control
# this is probably not a good negative control. We see many significant hits, but they're all due to blocks of expremely low/missing values
# put model first to avoid adding X to te names that start with numerics
coefs2  <- coefs %>%
  dplyr::select(model,
                id = gene_symbol,
                feature.y,
                pscore=p.value_anova_x,
                cscore=estimate.qnorm)  %>%
  mutate(pscore = -log10(pscore),
         cscore = abs(cscore)) %>%
  filter(id %in% tree$tip.label) %>%
  pivot_longer(c(pscore, cscore),names_to = "score_type",values_to = "score") %>%
  pivot_wider(names_from = c(model,score_type,feature.y), values_from = score)

dat <- tree$tip.label %>%
  tibble(id = .) %>%
  left_join(coefs2, by="id") %>%
  distinct() %>%
  mutate(across(where(is.numeric),replace_na,0)) %>%
  mutate(across(where(is.numeric),abs)) %>%
  column_to_rownames("id")

dat$random <- rnorm(dim(dat)[1],sd=10)
dat$bm <- rTraitCont(tree)

p4d <- phylo4d(tree, tip.data=dat,rownamesAsLabels=T)

res.ps <- phyloSignal(p4d)

# take the simple list of results and make into a tbl for easy querying
res.ps.tbl <- res.ps$pvalue %>%
  as_tibble(rownames="coef") %>%
  separate(coef,into=c("sex","score","TE"),sep="_",extra = "merge",remove = F) %>%
  mutate(sex = ifelse(is.na(sex),coef,sex),
         score = ifelse(is.na(score),coef,score),
         TE = ifelse(is.na(TE),coef,TE)) %>%
  pivot_longer(-c(TE,sex,coef,score),names_to = "metric",values_to = "pval") %>%
  group_by(sex,score,metric) %>%
  mutate(padj = p.adjust(pval,method="BH")) %>%
  ungroup()
  

hits <- res.ps.tbl %>%
  filter(!coef %in% c("bm","random")) %>%
  filter(score=="pscore") %>%
  arrange(pval) %>%
  filter(padj < 0.1) %>% 
  arrange(sex,TE)

f_hits2use <- filter(hits,score=="pscore" & sex == "female" & metric == "Lambda")$coef
m_hits2use <- filter(hits,score=="pscore" & sex == "male" & metric == "Lambda")$coef

m_carni.lipa <- lipaMoran(p4d, trait = c(m_hits2use, f_hits2use),alternative = "greater")


crlg0 <- phyloCorrelogram(p4d, trait = c("female_pscore_Doc"))
crlg1 <- phyloCorrelogram(p4d, trait = c("female_pscore_G4"))
crlg2 <- phyloCorrelogram(p4d, trait = c("female_pscore_Ivk"))
crlg3 <- phyloCorrelogram(p4d, trait = c("female_pscore_Stalker"))
crlg4 <- phyloCorrelogram(p4d, trait = c("female_pscore_Stalker2"))
crlg5 <- phyloCorrelogram(p4d, trait = c("female_pscore_gypsy9"))
crlg6 <- phyloCorrelogram(p4d, trait = c("female_pscore_invader2"))
crlg7 <- phyloCorrelogram(p4d, trait = c("female_pscore_transib1"))


plot(crlg0)
plot(crlg1)
plot(crlg2)
plot(crlg3)

plot(crlg7)

# =======================================================================


library(tidyverse)


zadtree <- ape::read.tree("results/zad/iqtree/domain-hits/PF07776/PF07776.contree")

zad <- zadtree$tip.label |> str_extract(".+(?=\\.)")

tfs <- read_tsv("resources/Drosophila_melanogaster_TF.txt")

zads <- read_tsv("resources/Drosophila_melanogaster_TF.txt") |>
  filter(Ensembl %in% zad) |>
  pull(Symbol) |>
  unique()

other_tfs <- read_tsv("resources/Drosophila_melanogaster_TF.txt") |>
  filter(!Ensembl %in% zad) |>
  pull(Ensembl) |>
  unique()

tree <- read_rds("results/te_sequence_similarity/te_sketch_dist.rds") |>
  filter(te1 %in% coefs$feature.y & te2 %in% coefs$feature.y) |>
  pivot_wider(names_from = "te1", values_from = "dissimilarity") |>
  column_to_rownames("te2") |>
  as.matrix() |>
  as.dist() |>
  ape::bionjs() |>
  ape::as.phylo()

coefs2  <- coefs %>%
  filter(feature.x %in% tfs$Ensembl | feature.x %in% zad) |>
  dplyr::select(model,
                id = feature.y,
                gene_symbol,
                cscore=estimate.qnorm)  %>%
  filter(id %in% tree$tip.label) %>%
  pivot_longer(c(cscore),names_to = "score_type",values_to = "score") %>%
  pivot_wider(names_from = c(model,score_type,gene_symbol), values_from = score)

dat <- tree$tip.label %>%
  tibble(id = .) %>%
  left_join(coefs2, by="id") %>%
  distinct() %>%
  mutate(across(where(is.numeric),replace_na,0)) %>%
  column_to_rownames("id")

x <- dat[,str_detect(colnames(dat),"female")]
colnames(x) <- str_remove(colnames(x),"female_cscore_")


dat$random <- rnorm(dim(dat)[1],sd=10)
dat$bm <- rTraitCont(tree)
dat$zad_mean <- rowMeans(x[,colnames(x) %in% zads])[rownames(dat)] 

p4d <- phylo4d(tree, tip.data=dat,rownamesAsLabels=T)

res.ps <- phyloSignal(p4d)

# take the simple list of results and make into a tbl for easy querying

res.stat.tbl <- res.ps$stat %>%
  as_tibble(rownames="coef") %>%
  separate(coef,into=c("sex","score","TE"),sep="_",extra = "merge",remove = F) %>%
  mutate(sex = ifelse(is.na(sex),coef,sex),
         score = ifelse(is.na(score),coef,score),
         TE = ifelse(is.na(TE),coef,TE)) %>%
  pivot_longer(-c(TE,sex,coef,score),names_to = "metric",values_to = "stat")

res.ps.tbl <- res.ps$pvalue %>%
  as_tibble(rownames="coef") %>%
  separate(coef,into=c("sex","score","TE"),sep="_",extra = "merge",remove = F) %>%
  mutate(sex = ifelse(is.na(sex),coef,sex),
         score = ifelse(is.na(score),coef,score),
         TE = ifelse(is.na(TE),coef,TE)) %>%
  pivot_longer(-c(TE,sex,coef,score),names_to = "metric",values_to = "pval") %>%
  group_by(sex,score, metric) %>%
  mutate(padj = p.adjust(pval,method="BH")) %>%
  ungroup()


res.tbl <- full_join(res.stat.tbl, res.ps.tbl, by = join_by(coef, sex, score, TE, metric)) |>
  filter(score == "cscore")

ggplot(res.tbl, aes(stat)) +
  geom_histogram() +
  facet_wrap(~metric)

hits <- res.tbl |> 
  filter(score == "cscore" & metric !="I") |> 
  group_by(sex, metric) |>
  mutate(scaled.stat = scale(stat)[,1]) |>
  group_by(sex,TE) |>
  summarise(n_sig_metrics= sum(padj < 0.1), mean_scaled.stat = mean(scaled.stat),.groups = "drop") |>
  filter(n_sig_metrics >=3) |>
  arrange(-abs(mean_scaled.stat)) |>
  mutate(lab = sprintf("%s_cscore_%s", sex, TE))

hits|>
  mutate(is.zad = TE %in% zads) |>
  #filter(is.zad)
  ggplot(aes(sex, fill=is.zad)) +
  geom_bar(position = "stack") +
  labs(title="significant using >2/5\nphylosignal metrics")

#-------
#x <- x[,zads[zads %in% colnames(x)]]
#x <- x[,colnames(x) %in% head(hits$TE,50)]

hc_tfs <- dist(t(x)) |> hclust(method="ward.D")


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

barplot.phylo4d(p4d,trait = c("female_cscore_CG11902",
                              "female_cscore_nau"), center=F, scale=F, main="", ladderize=T)


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





hc <- dat[, useful_hits] |>
  t() |>
  dist() |>
  hclust(method="ward.D2")


plot(hc, main= "a")

gridplot(p4d, 
         trait = useful_hits[hc$order],
         cell.col = viridis::rocket(1000),tree.ladderize=F)




gridplot(p4d, trait =  useful_hits,
         tree.type = "fan",
         tip.cex = 0.6, tree.open.angle = 160, trait.cex = 0.6,
         show.trait = T, show.tip = T,trait.font=1,
         cell.col = viridis::magma(10))

# =========================================
barplot.phylo4d(p4d,trait = useful_hits, bar.col = (m_carni.lipa$p.value < 0.05) + 1, center=T, scale=T)

barplot.phylo4d(p4d,trait = unique(hits$coef), center=T, scale=T)

gridplot(p4d, 
        trait = useful_hits,
        cell.col = viridis::magma(1000),tree.ladderize=F)

gridplot(p4d, 
         trait = f_hits2use,
         cell.col = viridis::magma(1000),tree.ladderize=F)



gridplot(carni.lipa.p4d, trait =  c(unique(hits$coef)),
         tree.type = "fan",
         tip.cex = 0.6, tree.open.angle = 160, trait.cex = 0.6,
         show.trait = T, show.tip = T,trait.font=1,
         cell.col = viridis::magma(10))
