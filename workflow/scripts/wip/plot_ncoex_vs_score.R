library(tidyverse)

mods2 <- read_tsv("upstream/final-models.collected-info.tsv.gz")

ncoex <- mods2 |> 
  filter(significant_x) |> 
  dplyr::select(gene_symbol,feature.y) |>
  distinct() |>
  count(gene_symbol,sort = T) |>
  mutate(n_coex_rank = row_number())


ggplot(ncoex, aes(n)) +
  geom_histogram()


ncoex |> filter(gene_symbol == "pan")  

mods2 |> 
  filter(gene_symbol == "pan") |>
  filter(significant_x)

  
head(ncoex, n=100) |> print(n=Inf)
