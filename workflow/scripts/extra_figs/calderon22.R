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
library(wCorr);library(psych) # for weighted correlation
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

lms <- lms |>
  mutate(signed_prop_x_explained_var = sign(estimate.qnorm) * sumsq_anova_x/total_variance)

relationship_df <- lms |>
  filter(valid & adj_p.value_ftest_r2 < 0.1) |>
  mutate(class = case_when(significant_x ~ "coex",
                           str_detect(overlap_filter_type,"direct") ~"fixed.overlap",
                           str_detect(overlap_filter_type,"correlated") ~ "corr.gene.fixed.overlap",
                           !valid | !adj_p.value_ftest_r2 < 0.1 | !significant_x~ "ns",
                           T ~ "other.untested"),.groups = "drop") |>
  dplyr::select(sex=model,gene_symbol, feature.y.scrna,class, estimate.qnorm, signed_prop_x_explained_var) |>
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

GE_tf_te_ix <- (rownames(GE) %in% c(tfs,tes,"piwi"))
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

# get weighted correlation
df <- mutate(df, cor.raw = map2(GE_mat, SC_sizes, get_weighted_cor))

# get average expression used for spqn
df <- mutate(df, ave_GE = map(GE_mat, ~log1p(rowMeans(expm1(.x)))))

# get spqn-corrected weighted correlation
df <- mutate(df, cor.spqn = map2(cor.raw, ave_GE, get_corrected_cor))

# get pvals and coefs and annotations in tbl format
df <- mutate(df, res.raw = map2(cor.raw, n, cormat2tbl),
             res.spqn  = map2(cor.spqn, n, cormat2tbl) )

df_lineages <- df
df <- filter(df, lineage %in% c("all_cells"))

# ------------------------------------------------------------------------------
# qc plots
# ------------------------------------------------------------------------------

g_supercell_size <- SC$supercell_size |>
  enframe(value = "n_cells") |>
  ggplot(aes(n_cells)) +
  geom_histogram() +
  scale_x_log10()

gg_spqn <- df |>
  filter(lineage=="all_cells") |>
  mutate(gg_coefs_raw = map2(cor.raw, ave_GE, plot_signal_condition_exp, signal=0.01),
         gg_coefs_spqn = map2(cor.spqn, ave_GE, plot_signal_condition_exp, signal=0.01)) |>
  dplyr::select(lineage,gg_coefs_raw, gg_coefs_spqn)


g_all_cell_corr_pvals <- df |>
  filter(lineage=="all_cells") |>
  dplyr::select(lineage, res=res.spqn) |>
  unnest(res) |>
  filter(sex=="female") |> # sex only deals with the LM coefs here, so we can just look at one sex
  ggplot(aes(p)) +
  geom_histogram() +
  scale_fill_grey()


g_all_cell_overlapping_features_coexpressed <- df |> 
  filter(lineage=="all_cells") |>
  dplyr::select(lineage,res=res.spqn) |>
  unnest(res) |>
  drop_na() |>
  ggplot(aes(class,coef)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ggpubr::stat_compare_means(ref.group = "ns",method = "wilcox.test") +
  facet_wrap(~sex) +
  coord_cartesian(ylim = c(-0.1,0.1))

#pan is correlated with most TEs at a high level.
#We can look at enrichment of coexpressed TEs (from DGRP data) in our coexpression here.
#Note that I limit two two categories (ns and coex) and exclude the pairs filtered out because they appear to have fixed insertions that I can't control for.
#And again, consider that I have excluded TE/TF pairs that were unmodelable in the DGRP data (non-valid models or non-significant fit).
#we also don't consider sex here; a te is coexpressed if it is significant in either sex.


get_cont_mat <- function(x) {
  message("nother")
  z <- filter(x, class %in% c("ns","coex")) |>
    group_by(feature,y) |>
    summarise(class=if_else(any(class=="coex"),"coex","ns"),padj=unique(padj), .groups="drop") |>
    group_by(class,sc.class= if_else(padj < 0.1,"sc.coex","sc.ns")) |>
    tally() |>
    ungroup() |>
    pivot_wider(names_from = sc.class, values_from = n, values_fill = 0) |>
    arrange(desc(class)) |>
    column_to_rownames("class")
  
  z[c("coex","ns"),c("sc.coex","sc.ns")]
}

possibly_get_cont_mat <-  possibly(get_cont_mat)
possibly_fisher <-  possibly(\(x) broom::tidy(fisher.test(x)))

sc_lm_enrich_df <- df_lineages |> 
  dplyr::select(lineage,res=res.spqn) |>
  unnest(res) |>
  #filter(feature == "pan") |>
  nest(data=-c(lineage)) |>
  mutate(cont_mat = map(data,possibly_get_cont_mat)) |>
  mutate(fisher_res = map(cont_mat, possibly_fisher)) |>
  unnest(fisher_res)

lineage_enrichments <- sc_lm_enrich_df |> 
  arrange(p.value) |>
  filter(p.value < 0.05 & estimate > 1) |>
  dplyr::select(-data)


pull(lineage_enrichments, cont_mat)[[1]] |>
  rownames_to_column("DGRP class") |>
  gt::gt(rowname_col = "DGRP class")

dplyr::select(lineage_enrichments, -cont_mat) |>
  gt::gt()

g_pan_highly_corr_with_tes <- df %>%
  dplyr::select(lineage,res=res.spqn) |>
  unnest(res) |>
  filter(sex %in% c("female","male")) |> # only look at TEs/feature pairs that were modelable in male or female DGRP data
  mutate(feature2 = if_else(feature=="pan","pan","other")) |>
  ggplot(aes(coef, fill=feature2)) +
  geom_density() +
  facet_wrap(~sex)

# ------------------------------------------------------------------------------
# pan coexpression scatters
# ------------------------------------------------------------------------------

plot_indiv_relationship <- function(x, y, caption="") {
  supercell_GeneGenePlot(GE_all,gene_x = x, gene_y = y,supercell_size = SC$supercell_size,color.use = "black") -> res
  
  g <- res$p$data |> 
    as_tibble() |>
    ggplot(aes(x, y, size=size)) +
    ggdensity::geom_hdr_points(method="mvnorm") +
    xlab(x) +
    ylab(y) +
    labs(title=sprintf("%s vs %s",x,y),subtitle = sprintf("raw weighted corr:%s; p=%s",format.pval(res$w.cor[[1]],digits = 3),format.pval(res$w.pval[[1]],digits = 3)),caption = caption) +
    geom_smooth(method="lm",se=F,color="red",linetype="dashed",linewidth=1)
  
  return(g)
}

g_poscon_and_hits_panel <- plot_indiv_relationship("pan","nkd","strongly Wnt controlled") +
  plot_indiv_relationship("pan","1360","Fixed insertion") +
  plot_indiv_relationship("pan","rooA"," Fixed ins. correlated genes") +
  plot_indiv_relationship("pan","invader2","DGRP coexpressed") +
  plot_indiv_relationship("pan","Stalker4","DGRP coexpressed") +
  plot_indiv_relationship("pan","S","Coexpressed in scRNA, not DGRP") + patchwork::plot_layout(nrow=2, guides = "collect") &
  guides(size="none")


# negcon panel selected from highly expressed embryo genes (Modencode RNA-seq via flybase)
g_negcon_panel <- plot_indiv_relationship("pan","ci") +
  plot_indiv_relationship("pan","HmgD") +
  plot_indiv_relationship("pan","awd") +
  plot_indiv_relationship("pan","EAChm") +
  plot_indiv_relationship("pan","tin") +
  plot_indiv_relationship("pan","Odj") +patchwork::plot_layout(nrow=2,guides = "collect") &
  guides(size="none")

# simplified trio of plots for figures v2 (simplified figs)
g_pan_vs <- c("nkd","1360","invader2") |>
set_names() |>
  map(~plot_indiv_relationship("pan",.x))


write_rds(df_lineages, snakemake@output$df)
write_rds(g_supercell_size, snakemake@output$g_supercell_size)
write_rds(gg_spqn,snakemake@output$gg_spqn)
write_rds(g_all_cell_corr_pvals, snakemake@output$g_all_cell_corr_pvals)
write_rds(g_all_cell_overlapping_features_coexpressed, snakemake@output$g_all_cell_overlapping_features_coexpressed)
#lineage_enrichments # note the sign doesn't always agree
write_rds(g_pan_highly_corr_with_tes, snakemake@output$g_pan_highly_corr_with_tes)
write_rds(g_poscon_and_hits_panel,snakemake@output$g_poscon_and_hits_panel)
write_rds(g_negcon_panel,snakemake@output$g_negcon_panel)
write_rds(g_pan_vs,snakemake@output$g_pan_vs)

