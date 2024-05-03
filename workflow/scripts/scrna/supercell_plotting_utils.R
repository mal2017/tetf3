Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

# get a correlation weighted by the number of cells represented in each supercell/metacell
get_weighted_cor <- function(x, n) {
  cm <- t(x) |>
    as.matrix() |>
    psych::cor.wt(w = n)
  
  cm$r
}

# use spqn to correct correlation values
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

# get p value from a correlation value
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

# make contingency table from 
get_cont_mat <- function(x) {
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

# allow failure when contingency table would have missing columns
possibly_get_cont_mat <-  possibly(get_cont_mat)
possibly_fisher <-  possibly(\(x) broom::tidy(fisher.test(x)))


# make scatterplot of 1 feature vs another
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
