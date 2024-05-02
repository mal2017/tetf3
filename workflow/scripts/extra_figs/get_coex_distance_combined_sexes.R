library(tidyverse)

mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]
mods <- read_tsv(mods_path)

tfs_fl <- "resources/Drosophila_melanogaster_TF.txt"
tfs_fl <- snakemake@input$tfs
tfs <- unique(read_tsv(tfs_fl) |> pull(Ensembl))

# noted that animaltfdb misses some zad genes
zads_fl <- "results/resources/zad_genes.tsv"
zads_fl <- snakemake@input$zads
zads <- read_tsv(zads_fl)$ensembl_gene_id

coex_mat <- mods |>
  filter(feature.x %in% tfs | feature.x %in% zads) |>
  nest(data=-model) |>
  mutate(dist = map(data,.f=function(x) {
      dplyr::select(x, feature.x, feature.y, score=estimate.qnorm) %>%
        pivot_wider(names_from = c(feature.x), values_from = score, values_fill = 0) |>
        column_to_rownames("feature.y") |>
        as.matrix() |>
        dist(method="euclidean")
  })) |>
  mutate(tab = map(dist, .f=function(x) {
    as.matrix(x) |>
      as_tibble(rownames="te1") |>
      pivot_longer(-te1, names_to = "te2",values_to = "dissimilarity") |>
      filter(te1!=te2)
  })) |>
  dplyr::select(-data)


saveRDS(coex_mat, snakemake@output$rds)
saveRDS(filter(coex_mat,model=="female")$dist[[1]], snakemake@output$female_dist)
saveRDS(filter(coex_mat,model=="male")$dist[[1]], snakemake@output$male_dist)