library(phylosignal)
library(ape)
library(adephylo)
library(phylobase)
library(phylosignal)
library(tidyverse)
library(TreeTools)
library(R.devices)

#https://cran.r-project.org/web/packages/phylosignal/vignettes/Demo_plots.htmllibr

rnk <- ifelse(exists("snakemake"), snakemake@input$rnk, "results/rankings/nofilt_main_male_mean_abs_estimate_qnorm.tsv.gz") %>%
    read_tsv()

lkup <- read_tsv("http://ftp.flybase.net/releases/FB2021_04/precomputed_files/genes/fbgn_annotation_ID.tsv.gz",skip = 4) %>%
  dplyr::select(gene_id = `primary_FBgn#`,gene_symbol = `##gene_symbol`) %>%
  distinct() %>%
  deframe()

# import the representative or consensus tree
tree <- ifelse(exists("snakemake"),snakemake@input$contree,"results/zad/iqtree/domain-hits/PF07776/PF07776.contree") %>%
  ape::read.tree()

# import bootstrapped trees
boots <- ifelse(exists("snakemake"),snakemake@input$boots,"results/zad/iqtree/domain-hits/PF07776/PF07776.ufboot") %>%
  ape::read.tree()

rename_tree_nodes <- function(x) {
  # get a weird warning later if itnernal nodes aren't unique
  # the example in phylosignal doesn't even have itnernal node names
  # so I remove them
  x$node.label <- NULL
  x$tip.label <- x$tip.label %>% str_extract(".+(?=\\.FBpp)") %>% lkup[.]
  return(x)
}

tree <- rename_tree_nodes(tree)

boots <- boots %>% map(rename_tree_nodes) %>% as.multiPhylo()

dat <- tree$tip.label %>%
  tibble(id = .) %>%
  left_join(rnk, by=c(id="gene_symbol")) %>%
  distinct() %>%
  mutate(across(where(is.numeric),replace_na,0)) %>%
  mutate(across(where(is.numeric),abs)) %>%
  dplyr::select(-gene_id) %>%
  column_to_rownames("id")

# get internal controls
set.seed(2142)
dat$random <- rnorm(dim(dat)[1],sd=10)
dat$bm <- rTraitCont(tree)

# combine data and phylogeny into single object
p4d <- phylo4d(tree, tip.data=dat,rownamesAsLabels=T)

# run phylosignal and get it for boots
res.ps <- phyloSignal(p4d)
R.devices::suppressGraphics({
  res.ps.boots <- phyloSignalBS(p4d,boots)
})


make_pc_res_tbl <- function(x) {
  x$pvalue %>%
    as_tibble(rownames="score") %>%
    pivot_longer(-c(score),names_to = "metric",values_to = "pval") %>%
    left_join(pivot_longer(as_tibble(x$stat,rownames="score"),-c(score),names_to = "metric",values_to = "statistic"))
}

# take the simple list of results and make into a tbl for easy querying
res.ps.tbl <- res.ps %>% make_pc_res_tbl() %>% mutate(boot = "0")
res.ps.tbl <- res.ps.boots %>% 
  set_names(.,as.character(1:length(.))) %>%
  map_df(make_pc_res_tbl, .id="boot") %>%
  bind_rows(res.ps.tbl)

res.ps.tbl %>%
  mutate(score= fct_reorder(score,pval)) %>%
  ggplot(aes(score, -log10(pval))) +
  geom_boxplot(data = . %>% filter(boot!=0), outlier.shape = NA) +
  geom_point(data = . %>% filter(boot == 0), color="red", size=3) +
  facet_wrap(~metric,nrow = length(unique(res.ps.tbl$metric))) +
  geom_hline(yintercept = -log10(0.05), color="red",linetype="dashed")

# list of hits we'll call significant. expecting that bm will be sig and random won't
hits <- res.ps.tbl %>% filter(boot == 0 & pval < 0.1 & !score %in% c("bm","random"))

# simulate p values with pos con vector (brownian motion
#phylosim <- phyloSim(tree = tree, method = "all", nsim = 100, reps = 99)
#plot.phylosim(phylosim, what = "pval", stacked.methods = TRUE)

# correlogram
crlg <- phyloCorrelogram(p4d, trait = "value",ci.bs = 1000,n.points = 100,ci.conf = 0.95)
plot(crlg)

barplot.phylo4d(p4d,trait = c("value","bm","random"),center=F, scale=F, tree.ladderize = T)

carni.lipa <- lipaMoran(p4d)
carni.lipa.p4d <- lipaMoran(p4d, as.p4d = TRUE)

barplot.phylo4d(p4d, 
                trait = "value",
                bar.col = (carni.lipa$p.value < 0.05) + 1, center = T, scale = T, tree.ladderize = T)

gridplot(p4d, 
         trait =  "max_coex",
         tree.type = "fan",
         tip.cex = 0.6, tree.open.angle = 160, trait.cex = 0.6,
         show.trait = T, show.tip = T,trait.font=1,
         cell.col = viridis::magma(100))
