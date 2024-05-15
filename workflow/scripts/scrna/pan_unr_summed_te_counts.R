Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(scater)
library(scran)
library(scuttle)
library(jsonlite)
library(psych)

coex <- read_json(ifelse(exists("snakemake"),snakemake@input$json,"results/resources/coexpressed_tes.json"))

detected_tes <- coex$all |> unlist() |> unique()

bound <- read_tsv(ifelse(exists("snakemake"),snakemake@input$rip,"results/ripseq/unr_ripseq.tsv.gz")) |>
  filter(type=="TE" & status == "bound") |>
  pull(feature)

# ------------------------------------------------------------------------------
sce <- read_rds(ifelse(exists("snakemake"),snakemake@input$sce,'resources/all.usa.fca_annot.supercell.sce.rds'))

# see https://support.bioconductor.org/p/9154078/ and 
# description of 'log' argument of lognormcounts (called upstream via batchelor::multiBatchNorm)
normcounts(sce) <- (2^logcounts(sce))-1

metadata_df <- normcounts(sce) |>
  as.matrix() |>
  colSums() |>
  enframe(name = "metacell",value = "colsums")

metadata_df <- left_join(metadata_df,as_tibble(colData(sce),rownames="metacell"), by="metacell")

# check colsum dist
#ggplot(metadata_df,aes(colsums)) +
#  geom_histogram()

# check relationship of size to sum
#ggplot(metadata_df,aes(size,colsums)) +
#  geom_point() +
#  scale_x_log10() +
#  scale_y_log10()

# ------------------------------------------------------------------------------
# get aggregated (summed) counts for gene sets and individual genes of interest

# generate sce with only TEs and genes of interest to make 'invert' arg work as expected
sce.sub <- sce[c("FBgn0085432","FBgn0263352",rownames(sce)[!str_detect(rownames(sce),"FBgn")])]

# fun to get
get_summed_normcounts <- function(gene_set,gsname="x", invert=F) {
  
  if (invert) {
    sx <- aggregateAcrossFeatures(sce.sub,if_else(!rownames(sce.sub) %in% gene_set,"summed.normcounts","other"), average=F, use.assay.type = "normcounts")  
  } else {
    sx <- aggregateAcrossFeatures(sce.sub,if_else(rownames(sce.sub) %in% gene_set,"summed.normcounts","other"), average=F, use.assay.type = "normcounts")
  }
  
  sx["summed.normcounts",] |>
    normcounts() |>
    as.matrix() |>
    t() |>
    as_tibble(rownames="metacell") |>
    mutate(gene_set = gsname) |>
    dplyr::relocate(gene_set,.after="metacell")
}



summed_normcounts <- bind_rows(get_summed_normcounts(rownames(sce)[!str_detect(rownames(sce),"FBgn")],"TEs"),
                               get_summed_normcounts("FBgn0085432","pan"),
          get_summed_normcounts(coex$all$pan,"pan.coexpressed"),
          get_summed_normcounts(coex$all$pan,"pan.noncoexpressed",invert = T),
          get_summed_normcounts("FBgn0263352","Unr"),
          get_summed_normcounts(coex$all$Unr,"Unr.coexpressed"),
          get_summed_normcounts(coex$all$Unr,"Unr.noncoexpressed",invert = T),
          get_summed_normcounts(bound,"Unr.bound"),
          get_summed_normcounts(bound,"Unr.nonbound",invert = T)) |>
  pivot_wider(names_from = "gene_set", values_from = "summed.normcounts") |>
  left_join(metadata_df,y=_,by="metacell")


# test the function
stopifnot(all(dplyr::select(summed_normcounts,pan)$pan == normcounts(sce)["FBgn0085432",]))


write_tsv(summed_normcounts,snakemake@output$tsv)
