
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

group_ranking <- ifelse(exists("snakemake"),snakemake@input$group_ranking,"results/rankings/te_scores_by_groupings.tsv.gz") %>% read_tsv()

# filtered, max is too sparse - we get obviously spurious signal
coefs2 <- group_ranking %>% 
  filter(grouping_level == "superfamily" & filtering_approach == "nofilt") %>%
  dplyr::select(id=gene_symbol, model, grouping, score ="mean_abs_estimate.qnorm") %>%
  pivot_wider(names_from = c(model, grouping), values_from = "score") %>%
  mutate(across(where(is.numeric),function(x){replace_na(x,0)}))

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
  mutate(padj = p.adjust(pval,method="BH")) %>%
  ungroup()
  

hits <-  res.ps.tbl %>%
  filter(!coef %in% c("bm","random")) %>%
  arrange(pval) %>%
  filter(padj < 0.1) %>% 
  arrange(pval, sex,TE)

f_hits2use <- filter(hits,sex == "female")$coef %>% unique()
m_hits2use <- filter(hits,sex == "male")$coef %>% unique()

m_carni.lipa <- lipaMoran(p4d, trait = m_hits2use)
f_carni.lipa <- lipaMoran(p4d, trait =f_hits2use)


# correlogram
set.seed(1)
crlgM <- phyloCorrelogram(p4d, trait = m_hits2use[1])
plot(crlgM)

set.seed(1)
crlgBM <- phyloCorrelogram(p4d, trait = c("bm"))
plot(crlgBM)

set.seed(1)
crlgRND <- phyloCorrelogram(p4d, trait = c("random"))
plot(crlgRND)

barplot.phylo4d(p4d,trait = c(m_hits2use,f_hits2use,"bm","random"), center=F, scale=F)


