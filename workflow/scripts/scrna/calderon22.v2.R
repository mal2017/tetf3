library(tidyverse)
library(SuperCell)
library(scater)
library(scran)
library(corrr)
library(vegan)
library(HiClimR)
library(ggtree)
library(modelr)
library(DescTools)
library(psych) # for weighted correlation
library(spqn)
library(patchwork)

# this is a big script, because I think it would make it too complicated to split up into a million chunks to make a smk workflow.
# it takes a saved SuperCell object from the upstream workflow, plus assorted other files.

tfs <- read_tsv("resources/Drosophila_melanogaster_TF.txt")
tfs <- unique(tfs$Symbol)

tes <- jsonlite::read_json("upstream/te_element_lookup.json") %>%
  names()

lms <- read_tsv("upstream/final-models.collected-info.tsv.gz") |>
  mutate(feature.y.scrna  = str_remove(feature.y,"-element"))

relationship_df <- lms |>
  filter(valid) |>
  mutate(class = case_when(significant_x ~ "coex",
                           str_detect(overlap_filter_type,"direct") ~"fixed.overlap",
                           str_detect(overlap_filter_type,"correlated") ~ "corr.gene.fixed.overlap",
                           !valid | !adj_p.value_ftest_r2 < 0.1 | !significant_x~ "ns",
                           T ~ "other.untested"),.groups = "drop") |>
  dplyr::select(sex=model,gene_symbol, feature.y.scrna,class, estimate.qnorm) |>
  distinct()  

# ------------------------------------------------------------------------------
# import supercell object from upstream
# ------------------------------------------------------------------------------

SC_rds <- read_rds("upstream/calderon22_supercells.rds")

# extract the useful bits to their own objects
GE <- SC_rds$SC.GE
SC <- SC_rds$SC

# make informative cell names
colnames(GE) <- sprintf("hr%s.mc%s.%s",SC$window,1:SC$N.SC,SC$lineage)

GE_all <- GE

GE_tf_te_ix <- (rownames(GE) %in% c(tfs,tes))
GE_rowsds0_ix <- rowSds(GE) == 0

GE <- GE[GE_tf_te_ix & !GE_rowsds0_ix,]

# ------------------------------------------------------------------------------
# correlations
# ------------------------------------------------------------------------------


get_weighted_cor <- function(x, n) {
  cm <- t(x) |>
    as.matrix() |>
    psych::cor.wt(w = n)
  
  cm$r
}

get_corrected_cor <- function(x, ave_GE) {
  message("nother")
  ngrp <- 20; sizegrp <- 200; refgrp <- 6
  
  cormat.spqn <- spqn::normalize_correlation(x, 
                                             ave_exp = ave_GE,
                                             ngrp = ngrp,
                                             size_grp = sizegrp,
                                             ref_grp = refgrp)
  
  rownames(cormat.spqn) <- rownames(x)
  colnames(cormat.spqn) <- colnames(x)
  cormat.spqn
}


get_p_from_r <- function(r, n) {
  t_stat <- r *sqrt(n-2) / sqrt(1-r^2)
  # ripped from source code of stats::cor.test
  p <- 2 * min(pt(t_stat, n-2), pt(t_stat, n-2, lower.tail=FALSE))
  p
}

# annotates and reshapes correlatio coef mat into tbl with extra info and pvals
cormat2tbl <- function(x,n) {
  x |> 
    as_tibble(rownames = "feature") |>
    pivot_longer(-feature,names_to = "y", values_to = "coef") |>
    filter(feature %in% tfs & y %in% tes & feature!=y) |>
    mutate(p=map_dbl(coef,get_p_from_r,n=n)) |>
    mutate(padj = p.adjust(p,method="BH")) |>
    #filter(feature %in% c("pan")) |>
    left_join(relationship_df, by=c(feature="gene_symbol", y = "feature.y.scrna")) |>
    mutate(class = replace_na(class,"untested")) |>
    filter(!is.na(class)) |>
    mutate(class=fct_relevel(class,c("untested","ns","coex","corr.gene.fixed.overlap","fixed
                                     .overlap")))
} 

# raw info for calculating results
df <- split(SC$lineage,SC$lineage) |>
  c(list(all_cells = SC$lineage)) |>
  map(names) |>
  map(as.integer) |>
  enframe(name = "lineage", value = "SC_idx") |>
  mutate(n=map_dbl(SC_idx, length)) |>
  mutate(SC_sizes = map(SC_idx, ~{SC$supercell_size[.x]})) |>
  mutate(GE_mat = map(SC_idx, ~{GE[,.x]}))

df <- filter(df, lineage %in% c("all_cells"))

# get weighted correlation
df <- mutate(df, cor.raw = map2(GE_mat, SC_sizes, get_weighted_cor))

# get average expression used for spqn
df <- mutate(df, ave_GE = map(GE_mat, ~log1p(rowMeans(expm1(.x)))))

# get spqn-corrected weighted correlation
df <- mutate(df, cor.spqn = map2(cor.raw, ave_GE, get_corrected_cor))

# get pvals and coefs and annotations in tbl format
df <- mutate(df, res.raw = map2(cor.raw, n, cormat2tbl),
             res.spqn  = map2(cor.spqn, n, cormat2tbl) )


write_rds(df, snakemake@output$df)

#g_pan_highly_corr_with_tes <- 
#df %>%
#  dplyr::select(lineage,res=res.spqn) |>
#  unnest(res) |>
#  filter(sex %in% c("female","male")) |> # only look at TEs/feature pairs that were modelable in male or female DGRP data
#  dplyr::select(feature,y,coef,p,padj) |>
#    distinct() |>
#  mutate(feature2 = if_else(feature%in%c("pan","Unr","CG16779","vvl","NfI"),feature,"other")) |>
#  mutate(feature2 = fct_reorder(feature2,coef)) |>
#  mutate(feature2 = fct_relevel(feature2,"other")) |>
#  ggplot(aes(feature2,coef)) +
#    geom_boxplot()
