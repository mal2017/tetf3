Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

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

source("workflow/scripts/scrna/supercell_plotting_utils.R")

# get smk inputs ---------------------------------------------------------------

tfs_fl <- ifelse(exists("snakemake"),snakemake@input$tfs, "resources/Drosophila_melanogaster_TF.txt")

tes_fl <- ifelse(exists("snakemake"),snakemake@input$tes, "upstream/te_element_lookup.json")

lms_fl <- ifelse(exists("snakemake"),snakemake@input$mods, "upstream/final-models.collected-info.tsv.gz")

supercell_fl <- ifelse(exists("snakemake"),snakemake@input$supercell, "upstream/supercells.rds") 
# load data --------------------------------------------------------------------

tfs <- read_tsv(tfs_fl)
tfs <- unique(tfs$Symbol)

tes <- jsonlite::read_json(tes_fl) %>%
  names()


#Now I get the lms and categorize them, because we're looking for specific signatures to improve our confidence in this analysis.
#Note that currently (check code below in case I tweak this) I exclude 'untestable' coefs.

lms <- read_tsv(lms_fl) |>
  mutate(feature.y.scrna  = str_remove(feature.y,"-element"))

# might use signed variance explained by each gene's expression as a metric
lms <- lms |>
  mutate(signed_prop_x_explained_var = sign(estimate.qnorm) * sumsq_anova_x/total_variance)

# categorize each te/tf [air]
relationship_df <- lms |>
  filter(valid & adj_p.value_ftest_r2 < 0.1) |>
  mutate(class = case_when(significant_x ~ "coex",
                           str_detect(overlap_filter_type,"direct") ~"fixed.overlap",
                           str_detect(overlap_filter_type,"correlated") ~ "corr.gene.fixed.overlap",
                           !valid | !adj_p.value_ftest_r2 < 0.1 | !significant_x~ "ns",
                           T ~ "other.untested"),.groups = "drop") |>
  dplyr::select(sex=model,gene_symbol, feature.y.scrna,class, estimate.qnorm, signed_prop_x_explained_var) |>
  distinct()  

#Now get the supercell object we generated with upstream workflow tetf_calderon22
SC_rds <- read_rds(supercell_fl)


# process supercell results ----------------------------------------------------

# extract the useful bits to their own objects
GE <- SC_rds$SC.GE
SC <- SC_rds$SC

# make informative cell names
colnames(GE) <- sprintf("hr%s.mc%s.%s",SC$window,1:SC$N.SC,SC$lineage)

GE_all <- GE

GE_tf_te_ix <- (rownames(GE) %in% c(tfs,tes,"piwi"))
GE_rowsds0_ix <- rowSds(GE) == 0

GE <- GE[GE_tf_te_ix & !GE_rowsds0_ix,]

 

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


# plots ------------------------------------------------------------------------
#We want a balance between smoothing (high supercell size) and dropout-induced spurious correlations (small supercell size). Empirically, median supercell size of 10 is as low as is reasonable.


#When considering all cells, and all TE/TF pairs, these calls follow what we'd expect based on our lm correlations.

df |> 
  filter(lineage=="all_cells") |>
  dplyr::select(lineage,res=res.spqn) |>
  unnest(res) |>
  drop_na() |>
  ggplot(aes(class,abs(coef))) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ggpubr::stat_compare_means(ref.group = "ns",method = "wilcox.test") +
  facet_wrap(~sex)

df |> 
  filter(lineage=="all_cells") |>
  dplyr::select(lineage,res=res.spqn) |>
  unnest(res) |>
  filter(class %in% c("ns","coex")) |>
  group_by(feature, y, sex) |>
  summarise(coef=unique(coef),class=if_else(any(class == "coex"),"coex","ns")) |>
  drop_na() |>
  ggplot(aes(class,abs(coef))) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ggpubr::stat_compare_means(ref.group = "ns",method = "wilcox.test") +
  facet_wrap(~sex)

# But this doesn't quite work for pan's coefs, which aren't statistcally different between TEs we expect are coexpressed with pan and TEs we expect aren't.

df |> 
  filter(lineage=="all_cells") |>
  dplyr::select(lineage,res=res.spqn) |>
  unnest(res) |>
  filter(feature == "pan") |>
  ggplot(aes(class,abs(coef))) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ggpubr::stat_compare_means(ref.group = "ns",method = "t.test") +
  facet_wrap(~sex)

# But, it looks like that's because pan is correlated with most TEs at a high level.
# We can look at enrichment of coexpressed TEs (from DGRP data) in our coexpression here.
# Note that I limit two two categories (ns and coex) and exclude the pairs filtered out because they appear to have fixed insertions that I can't control for.
# And again, consider that I have excluded TE/TF pairs that were unmodelable in the DGRP data (non-valid models or non-significant fit).
# we also don't consider sex here; a te is coexpressed if it is significant in either sex.


sc_lm_enrich_df <- df |> 
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

# pan, across the board, has higher correlation coefficients with TEs than other TFs.

df %>%
  dplyr::select(lineage,res=res.spqn) |>
  unnest(res) |>
  #filter(sex == "female") |>
  mutate(feature2 = if_else(feature=="pan","pan","other")) |>
  ggplot(aes(coef, fill=feature2)) +
  geom_density() +
  facet_wrap(~sex)

## effect size vs corr
# agreement not amazing between lm effect sizes and correlation, but we have to keep in mind these are different tissues and different forces may be mediating coexpression at the population level and the tissue level. Additionally

df|>
  dplyr::select(lineage,res.spqn) |>
  unnest(res.spqn) |>
  filter(lineage=="all_cells") |>
  #filter(class =="coex") %>%
  #filter(p<0.05) %>%
  ggplot(aes(coef, signed_prop_x_explained_var)) +
  geom_point() +
  ggpubr::stat_cor(method="spearman") +
  theme(aspect.ratio = 1) +
  facet_wrap(~class) +
  geom_smooth(method="lm")

df|>
  dplyr::select(lineage,res=res.raw) |>
  unnest(res) |>
  filter(lineage=="all_cells") |>
  #filter(class =="coex") %>%
  #filter(p<0.05) %>%
  ggplot(aes(coef, estimate.qnorm)) +
  geom_point() +
  ggpubr::stat_cor(method="spearman") +
  theme(aspect.ratio = 1) +
  facet_wrap(~class) +
  geom_smooth(method="lm")

## pan coex tes
# here, class coex means that it is significantly coexpressed in adult lms. signed_prop_x_explained_variance is the proportion of variance in scaled, transformed, normalized expression explained by pan's expression in DGRP RNA-seq data.



plot_indiv_relationship("pan","nkd","strongly Wnt controlled") +
  plot_indiv_relationship("pan","1360","Fixed insertion") +
  plot_indiv_relationship("pan","rooA"," Fixed ins. correlated genes") +
  plot_indiv_relationship("pan","invader2","DGRP coexpressed") +
  plot_indiv_relationship("pan","Stalker4","DGRP coexpressed") +
  plot_indiv_relationship("pan","S","Coexpressed in scRNA, not DGRP") + patchwork::plot_layout(nrow=2) &
  guides(size="none")


# negcon panel selected from highly expressed embryo genes (Modencode RNA-seq via flybase)
plot_indiv_relationship("pan","ci") +
  plot_indiv_relationship("pan","HmgD") +
  plot_indiv_relationship("pan","awd") +
  plot_indiv_relationship("pan","EAChm") +
  plot_indiv_relationship("pan","tin") +
  plot_indiv_relationship("pan","Odj") +patchwork::plot_layout(nrow=2) &
  guides(size="none")
