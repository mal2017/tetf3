library(tidyverse)
library(phylosignal)
library(phylobase)
library(ape)

zad_fl <- "results/resources/zad_genes.tsv"
zad_fl <- snakemake@input$zads
zad <- read_tsv(zad_fl)

tfs_fl <- "resources/Drosophila_melanogaster_TF.txt"
tfs_fl <-  snakemake@input$tfs
tfs <- read_tsv(tfs_fl)

tree_fl <- "results/te_sequence_similarity/te_sketch_tidytree.rds"
tree_fl <- snakemake@input$tree
tree <- read_rds(tree_fl)

coefs_fl <- "~/work/tetf3/upstream/final-models.collected-info.tsv.gz"
coefs_fl <- snakemake@input$mods
coefs <- read_tsv(coefs_fl)

coefs2  <- coefs %>%
  filter(feature.x %in% tfs$Ensembl | feature.x %in% zad$ensembl_gene_id) |>
  mutate(estimate.qnorm = valid * estimate.qnorm) |> # if not a valid model, no business interpreting coefficients anyway
  dplyr::select(model,
                id = feature.y,
                gene_symbol,
                score=estimate.qnorm)  %>%

  filter(id %in% tree@phylo$tip.label) %>%
  pivot_longer(c(score),names_to = "score_type",values_to = "score") %>%
  pivot_wider(names_from = c(model,score_type,gene_symbol), values_from = score)

dat <- tree@phylo$tip.label %>%
  tibble(id = .) %>%
  left_join(coefs2, by="id") %>%
  distinct() %>%
  mutate(across(where(is.numeric),replace_na,0)) %>%
  column_to_rownames("id")

dat$random <- rnorm(dim(dat)[1],sd=mean(matrixStats::colSds(as.matrix(dat))))
set.seed(1)
dat$bm <- rTraitCont(tree@phylo,model = "BM")
set.seed(1)
dat$ou <- rTraitCont(tree@phylo,model = "OU")
dat$zad_mean <- rowMeans(dat[,str_remove(colnames(dat),".+male_score_") %in% zad$gene_symbol])[rownames(dat)] 

p4d <- phylo4d(tree@phylo, tip.data=dat[tree@phylo$tip.label,],rownamesAsLabels=T)

res.ps <- phyloSignal(p4d)

# take the simple list of results and make into a tbl for easy querying

res.stat.tbl <- res.ps$stat %>%
  as_tibble(rownames="coef") %>%
  separate(coef,into=c("sex","score_type","TF"),sep="_",extra = "merge",remove = F) %>%
  mutate(sex = ifelse(is.na(sex),coef,sex),
         score_type = ifelse(is.na(score_type),"control",score_type),
         TF = ifelse(is.na(TF),coef,TF)) %>%
  pivot_longer(-c(TF,sex,coef,score_type),names_to = "metric",values_to = "stat")

res.ps.tbl <- res.ps$pvalue %>%
  as_tibble(rownames="coef") %>%
  separate(coef,into=c("sex","score_type","TF"),sep="_",extra = "merge",remove = F) %>%
  mutate(sex = ifelse(is.na(sex),coef,sex),
         score_type = ifelse(is.na(score_type),"control",score_type),
         TF = ifelse(is.na(TF),coef,TF)) %>%
  pivot_longer(-c(TF,sex,coef,score_type),names_to = "metric",values_to = "pval") %>%
  group_by(sex, metric) %>%
  mutate(padj = p.adjust(pval,method="BH")) %>%
  ungroup()

res.tbl <- full_join(res.stat.tbl, res.ps.tbl, by = join_by(coef, sex, score_type, TF, metric))

write_rds(res.tbl, snakemake@output$tab)
write_rds(list(p4d = p4d, phylosignal=res.ps), snakemake@output$phylosignal)

