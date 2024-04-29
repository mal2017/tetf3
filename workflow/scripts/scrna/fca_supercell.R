# ------------------------------------------------------------------------------
# import supercell object from upstream
# ------------------------------------------------------------------------------

SC_rds <- readRDS("upstream/fca_supercells.rds")

# extract the useful bits to their own objects
GE <- SC_rds$SC.GE
SC <- SC_rds$SC
rm(SC_rds);gc()

library(tidyverse)
library(SuperCell)
library(psych) # for weighted correlation
library(spqn)

lkup <- read_tsv("results/resources/gene_symbol_lookup.tsv.gz") |>
  dplyr::select(gene_ID,gene_symbol) |>
  deframe()

tfs <- read_tsv("resources/Drosophila_melanogaster_TF.txt")

tes <- jsonlite::read_json("upstream/te_element_lookup.json") %>%
  names()

rownames(GE) <- map_chr(rownames(GE),~{if_else(.x %in% names(lkup),lkup[.x],.x)})

GE_all <- GE

#GE_tf_te_ix <- (rownames(GE) %in% c(tfs$Symbol,tes))
GE_rowsds0_ix <- rowSds(GE) == 0

GE <- GE[!GE_rowsds0_ix,]


# ------------------------------------------------------------------------------
# correlations
# ------------------------------------------------------------------------------
get_weighted_cor <- function(x, n) {
  cm <- t(x) |>
    as.matrix() |>
    psych::cor.wt(w = n)
  
  cm$r
}

# fun to get spqn-corrected corr coefs
get_corrected_cor <- function(x, ave_GE) {
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

# annotates and reshapes correlatio coef mat
cormat2tbl <- function(x,n) {
  x |> 
    as_tibble(rownames = "feature") |>
    pivot_longer(-feature,names_to = "y", values_to = "coef") |>
    filter(y %in% tes & feature!=y) |>
    mutate(p=map_dbl(coef,get_p_from_r,n=n)) |>
    mutate(padj = p.adjust(p,method="BH"))
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

print("checkpoint0")

# get weighted correlation
df <- mutate(df, cor.raw = map2(GE_mat, SC_sizes, get_weighted_cor))

print("checkpoint1")

# get average expression used for spqn
df <- mutate(df, ave_GE = map(GE_mat, ~log2(rowMeans(2^(.x)-1)+1))) # scater/scran/scanpy uses log2 by default

# get spqn-corrected weighted correlation
df <- mutate(df, cor.spqn = map2(cor.raw, ave_GE, get_corrected_cor))

print("checkpoint2")

# get pvals and coefs and annotations in tbl format
df <- mutate(df, res.raw = map2(cor.raw, n, cormat2tbl),
             res.spqn  = map2(cor.spqn, n, cormat2tbl) )

# export results
write_rds(df, snakemake@output$df)
write_rds(list(GE=GE_all,SC=SC),snakemake@output$supercell)