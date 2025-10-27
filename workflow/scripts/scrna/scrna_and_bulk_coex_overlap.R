library(tidyverse)


library(plotgardener)
library(tidyverse)
library(DiagrammeRsvg)
library(rsvg)
library(DiagrammeR)
library(SuperCell)
library(patchwork)

tfs <- read_tsv("resources/Drosophila_melanogaster_TF.txt")

tes <- jsonlite::read_json("upstream/te_element_lookup.json") %>%
  names()

mods <- read_tsv("upstream/final-models.collected-info.tsv.gz")

# read in TE TF correlations we identifed from the supercell results
feature_correlations <- read_rds("results/calderon22/fca_reanalysis_correlations.rds") |>
  dplyr::select(lineage,res=res.spqn) |>
  unnest(res)

scrna <- feature_correlations |>
  group_by(feature,y) |>
  summarise(scrna_coex = padj < 0.1,.groups = "drop")

bulk <- mods |>
  #filter(model == "female") |>
  group_by(gene_symbol,feature.y) |>
  summarise(bulk_coex = all(significant_x),.groups = "drop")


inner_join(scrna,bulk,by=c(feature="gene_symbol",y="feature.y"))  |>
  group_by(scrna_coex,bulk_coex) |>
  summarise(n=n()) |>
  mutate(scrna_coex = if_else(scrna_coex,"scRNA coex.","n.s"),
         bulk_coex = if_else(bulk_coex,"bulk coex.","n.s")) |>
  pivot_wider(names_from = bulk_coex,values_from = n) |>
  dplyr::relocate(`bulk coex.`,.after=scrna_coex) |>
  arrange(rev(scrna_coex)) |>
  column_to_rownames("scrna_coex") |>
  fisher.test()

