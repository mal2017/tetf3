Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)

mods <- "upstream/final-models.collected-info.tsv.gz"
mods <- snakemake@input$mods

tes <- mods |> read_tsv() |> pull(feature.y) |> unique()

#ifl <- "results/te_sequence_similarity/te_sketch_dist.txt"
ifl <- snakemake@input$txt

# if outputting  a similarity from dashing2, use this chunk
#x <- read_table(ifl, skip=2) |>
#  dplyr::rename(te1 = `#Sources`) |>
#  pivot_longer(-te1,names_to = "te2", values_to = "dissimilarity") |>
# mutate(dissimilarity = 1 - dissimilarity) # dash2 emits similarity by default

dat<- read_table(ifl, skip = 2) |>
  dplyr::rename(te1 = `#Sources`) %>%
  mutate(across(-te1,as.numeric))

#dashing_type <- "distance"
dashing_type <- snakemake@params$dashing_type

if (dashing_type == "distance") {
  
  dat <- dat |> # stalker3t and INE-1 don't like each other
    mutate(across(-te1,\(x){if_else(is.infinite(x),1.0,x)}))
  
} else if (dashing_type == "similarity") {
  dat <- dat |>
    mutate(across(-te1,\(x){1-x}))

} else {
  stop()
}

d <- dat |> column_to_rownames("te1") |> as.matrix()

d <- d[rownames(d) %in% tes,]
d <- d[,colnames(d) %in% tes]

d <-as.dist(d)

write_rds(d,snakemake@output$dist)