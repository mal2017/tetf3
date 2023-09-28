library(tidyverse)

# this script returns a list of two lists - one for males and one for female
# each of these has a dist object for TEs and a dist object for piRNA genes
# and a dataframe of the values used to generate these
# dist objects are returned so that it's easy to compare clustering methods 
# downstream, if desired

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]
mods <- read_tsv(mods_path) %>% filter(valid)

#pirna_path <- "results/resources/pirna_pathway.tsv"
pirna_path <- snakemake@input[["pirna"]]
pirna_tbl <- read_tsv(pirna_path)

# runs corr based dist clustering even with missing values
get_dist <- function(x) x |>
  scale() |>
  cor(use = "pairwise.complete.obs") |>
  as.data.frame() |> 
  mutate(across(everything(), replace_na, 0)) |> 
  as.matrix() %>%
  {sqrt(1-abs(.))} |>
  as.dist()

# takes the mod df and a sex argument, returns everything needed for a nice heatmap
make_hm_bundle <-  function(sex="female", x=mods) {
  res <- x %>%
    filter(feature.x %in% pirna_tbl$gene_ID) %>%
    filter(model == "female") |>
    dplyr::select(gene_symbol, feature.y, estimate.qnorm) |>
    pivot_wider(names_from = gene_symbol, values_from = estimate.qnorm) |>
    column_to_rownames("feature.y")
  
  # get rid of TEs missing coefs
  res <- res[which(!is.na(matrixStats::rowVars(as.matrix(res),na.rm = T))),]
  
  dist.pirnagenes <- res |> get_dist()
  dist.tes <- res |> t() |> get_dist()
  
  return(list(values = res, dist.pirnagenes = dist.pirnagenes, dist.tes=dist.tes))
}


c("male", "female") |>
  set_names() |>
  lapply(make_hm_bundle) |>
  saveRDS(snakemake@output[["rds"]])
