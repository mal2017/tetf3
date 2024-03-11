library(patchwork)
library(plotgardener)
library(tidyverse)
library(clusterProfiler)
library(ggpubr)
library(ggdensity)
library(patchwork)
library(DiagrammeRsvg)
library(rsvg)
library(vtree)
library(DiagrammeR)

# read in results of DGRP coexpression analysis
lms <- read_tsv("upstream/final-models.collected-info.tsv.gz") |>
  mutate(feature.y.scrna  = str_remove(feature.y,"-element"))

# categorize the types of relationships we uncovered from that analysis
relationship_df <- lms |>
  #filter(valid) |>
  mutate(class = case_when(significant_x ~ "coex",
                           str_detect(overlap_filter_type,"direct") ~"fixed.overlap",
                           str_detect(overlap_filter_type,"correlated") ~ "corr.gene.fixed.overlap",
                           !valid ~ "not.valid",
                           (!adj_p.value_ftest_r2 < 0.1) | (!significant_x)~ "ns",
                           T ~ "untested"),.groups = "drop") |>
  dplyr::select(sex=model,gene_symbol, feature.y.scrna,class, estimate.qnorm) |>
  distinct()  

# read in results of supercell analysis
supercell <- read_rds("results/calderon22/calderon22_reanalysis_supercell.rds")

# read in TE TF correlations we identifed from the supercell results
tf_te_correlations <- read_rds("results/calderon22/calderon22_reanalysis_correlations.rds") |>
  dplyr::select(lineage,res=res.spqn) |>
  unnest(res)


# ------------------------------------------------------------------------------
# plot supercell sizes to help explain that analysis
# ------------------------------------------------------------------------------

g_supercell_size <- supercell$SC$supercell_size |>
  enframe(value = "n_cells") |>
  ggplot(aes(n_cells)) +
  geom_histogram() +
  xlab("barcodes per supercell") + ylab("N")

# ------------------------------------------------------------------------------
# plot results of supercell correlation analysis
# ------------------------------------------------------------------------------

# asking if pan is highly correlated with TEs in general
g_pan_highly_corr_with_tes <- 
tf_te_correlations %>%
  dplyr::select(feature,y,coef,p,padj) |>
    distinct() |>
  mutate(feature2 = if_else(feature%in%c("pan","Unr"),feature,"other")) |>
  mutate(feature2 = fct_reorder(feature2,coef)) |>
  mutate(feature2 = fct_relevel(feature2,"other")) |>
  ggplot(aes(feature2,coef)) +    
  geom_boxplot() +
  ggpubr::stat_compare_means(ref="other",size=2) +
  xlab("") + ylab("weighted correlation")

# examine pan's coexpressed or not coexpressed TEs specifically
pan_te_tf_correlations <- tf_te_correlations |>
  filter(feature=="pan") |>
  inner_join(relationship_df, by=c(feature="gene_symbol",y="feature.y.scrna")) |>
  mutate(feature2 = if_else(feature%in%c("pan","Unr"),feature,"other")) |>
  mutate(feature2 = fct_reorder(feature2,coef)) |>
  mutate(class=fct_reorder(class,coef)) |>
  mutate(class=fct_relevel(class,"not.valid","ns","coex","corr.gene.fixed.overlap","fixed.overlap")) |>
  group_by(class,sex) |>
  mutate(n=n()) |>
  ungroup() |>
  mutate(lab = sprintf("%s (n=%s)",class,n)) |>
  mutate(lab=fct_reorder(lab,rank(class)))

g_pan_coex_tes_m <- pan_te_tf_correlations |>
  filter(sex=="male") |>
  ggplot(aes(lab,coef)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(comparisons = list(c(2,3)),size=2,label.y = 0.2) +
  #facet_wrap(~sex,scales="free") +
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  xlab("") + ylab("weighted correlation")


g_pan_coex_tes_f <- pan_te_tf_correlations |>
  filter(sex=="female") |>
  ggplot(aes(lab,coef)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(comparisons = list(c(2,3)),size=2,label.y = 0.2) +
  #facet_wrap(~sex,scales="free") +
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  xlab("") + ylab("weighted correlation")


# ------------------------------------------------------------------------------
# plot exemplary and control results from supercell analysis
# ------------------------------------------------------------------------------

plot_indiv_relationship <- function(x, y, caption="") {
  require(SuperCell)
  supercell_GeneGenePlot(supercell$GE,
                         gene_x = x, 
                         gene_y = y,
                         supercell_size = supercell$SC$supercell_size,
                         color.use = "black") -> res
  
  g <- res$p$data |> 
    as_tibble() |>
    ggplot(aes(x, y)) +
    ggdensity::geom_hdr_points(method="mvnorm",aes(size=size)) +
    xlab(x) +
    ylab(y) +
    labs(title=sprintf("%s vs %s",x,y),subtitle = sprintf("raw weighted corr:%s",round(res$w.cor[[1]],digits = 3)),caption = caption) +
    geom_smooth(method="lm",se=F,color="red",linetype="dashed",linewidth=1) +
    scale_size_continuous(name="barcodes/supercell") +
    scale_color_viridis_d(name="density")
  
  return(g)
}
# https://elifesciences.org/articles/68573
# summarizes refs for deposition of piwi in embryo - espec piwi
# which is also deposited in embryonic soma
# see franz et al 2017 (probing canonicity of wnt pathway) for selection
# of control genes
g_negcorr_poscon_pirna <- plot_indiv_relationship("piwi","1360")
g_pan_te_exemplary <- plot_indiv_relationship("pan","1360")
g_pan_poscon <- plot_indiv_relationship("pan","nkd")
g_pan_negcon <- plot_indiv_relationship("pan","CG4115")

g_controls <- g_pan_poscon + g_pan_negcon +
  g_negcorr_poscon_pirna + g_pan_te_exemplary + 
  plot_layout(nrow=1,guides = "collect") & 
  theme(legend.position = "bottom")


# ------------------------------------------------------------------------------
# overview of approach
# ------------------------------------------------------------------------------


g_a <- grViz("
digraph boxes_and_circles {

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]

  # several 'node' statements
  node [shape = box,
        fontname = Arial]
  data [label='processed embryo scRNA-seq data\n(Calderon et al. 2022)'];
  merge [label='collapse per-insertion UMI counts to TE family'];
  metacell [label='create metacells with SuperCell'];
  corr [label='calculate weighted TE/gene correlations'];

  data -> merge;
  merge -> metacell;
  metacell -> corr;
}
")

cartoon_temp <- tempfile()
#grVizToPNG(g_a, filename = "project_overview.png")
export_svg(g_a) |>
  charToRaw() |>
  rsvg_svg(cartoon_temp,width = 2000,height = 2000)

g_a_cartoon <- magick::image_read_svg(cartoon_temp) |> magick::image_ggplot(interpolate = T)

# plotting page 1 --------------------------------------------------------------------
theme_set(theme_classic() + theme(text=element_text(size=unit(7,"pt"))))

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 8.5, height = 11)

pageCreate(width =8.5, height = 11, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = g_a_cartoon, x = 1, y=0.3, width = 3, height=2.25)
plotText(label = "A", x = 1, y = 0.5)

pb <- plotGG(plot = g_supercell_size, x = 4.5, y=0.5, width = 3, height=2.25)
plotText(label = "B", x = 4.5, y = 0.5)

pc <- plotGG(plot = g_pan_highly_corr_with_tes, x = 0.5,  y=3, width = 2.75, height=2.25)
plotText(label = "C", x = 0.5, y = 3)

pd <- plotGG(g_pan_coex_tes_f,x=3.5,y=3,width=2, height=2.5)
plotText(label = "D", x = 3.5, y = 3)

pe <- plotGG(g_pan_coex_tes_m,x=5.75,y=3,width=2, height=2.5)
plotText(label = "E", x = 5.75, y = 3)

pf_i <- plotGG(g_controls,x=0.4,y=5.5,width=7.6,height=2.75)
plotText(label = "F", x = 0.5, y = 5.5)
plotText(label = "G", x = 2.4, y = 5.5)
plotText(label = "H", x = 4.15, y = 5.5)
plotText(label = "I", x = 6.15, y = 5.5)

dev.off()
