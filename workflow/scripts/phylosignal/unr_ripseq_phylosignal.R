library(tidyverse,quietly = T)
library(phylosignal,quietly = T)
library(phylobase,quietly = T)
library(ggtree)
library(tidytree)
library(ape)

# coexpression
coex_fl <- "~/work/tetf3/upstream/final-models.collected-info.tsv.gz"
coex_fl <- snakemake@input$coex
coex <- read_tsv(coex_fl) |>
  filter(gene_symbol == "Unr") |>
  filter(valid)

coex <- coex |>
  dplyr::select(model,estimate.qnorm,feature.y) |>
  pivot_wider(names_from = model,values_from = estimate.qnorm,values_fill = 0)

# rip signal
rip_fl <- "results/ripseq/unr_ripseq.tsv.gz"
rip_fl <- snakemake@input$rip
rip <- read_tsv(rip_fl) |> filter(type=="TE")

# get distance matrix
# filter to only include TEs with sufficient depth for inclusion in the upstream
# ripseq enrichment test
dist_fl <- "~/work/tetf3/results/te_sequence_similarity/te_female-coex_dist.rds"
dist_fl <- snakemake@input$dist
dist <- read_rds(dist_fl)
dist <- as.matrix(dist)
dist <- dist[rownames(dist) %in% rip$feature,colnames(dist) %in% rip$feature] |> as.dist()

# attta density - this is already within the region specified by relpos - that
# is handled by passing a param to the jop that produces this
attta_fl <- "results/ripseq/unr_ripseq_features_attta_sites.tsv.gz"
attta_fl <- snakemake@input$attta
density <- read_tsv(attta_fl) |>
  dplyr::select(feature,ARE_sites=sites)

# at content
relpos <- 50
relpos <- snakemake@params$relpos
at_content_fl <- "results/ripseq/unr_ripseq_features_au_content.tsv.gz"
at_content_fl <- snakemake@input$at_content
at_content <- read_tsv(at_content_fl)

lastXpct <- at_content |>
  filter(position.scaled>relpos) |>
  right_join(dplyr::select(rip,feature,group)) |>
  group_by(group, feature) |>
  summarise(nt.content=mean(nt.content,na.rm=T),.groups = "drop")

phylo <- ape::bionjs(dist)

tree <- as.treedata(phylo)

dat <- tree@phylo$tip.label %>%
  tibble(id = .) %>%
  left_join(dplyr::select(rip,feature,rip_stat,rip_log2FoldChange), by=c("id"="feature")) %>%
  distinct() %>%
  left_join(density,by=c("id"="feature")) |>
  left_join(lastXpct,by=c("id"="feature")) |>
  left_join(coex,by=c("id"="feature.y")) |>
  mutate(across(c(male,female),.fns=~replace_na(.x,0))) |>
  dplyr::select(-group) |>
  column_to_rownames("id")

tree <- tree |> as_tibble() |> left_join(as_tibble(dat,rownames="label")) |> as.treedata()

set.seed(1)
dat$random <- rnorm(dim(dat)[1],sd=sd(dat$rip_stat))

set.seed(1)
dat$bm <- rTraitCont(tree@phylo,model = "BM")

p4d <- phylo4d(tree@phylo, tip.data=dat[tree@phylo$tip.label,],rownamesAsLabels=T)

res.ps <- phyloSignal(p4d)

# take the simple list of results and make into a tbl for easy querying
res.tbl <- res.ps$pvalue %>%
  as_tibble(rownames="coef") %>%
  pivot_longer(-c(coef),names_to = "metric",values_to = "pval") %>%
  group_by(coef) %>%
  mutate(padj = p.adjust(pval,method="BH")) %>%
  ungroup()

print(res.tbl, n=Inf)

write_rds(res.tbl,snakemake@output$tbl)
write_rds(p4d,snakemake@output$p4d)
write_rds(tree,snakemake@output$tree)

