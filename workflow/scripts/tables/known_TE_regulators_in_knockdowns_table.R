library(tidyverse)
library(GenomicRanges)
library(gt)
library(writexl)

te_regulators <- read_tsv("results/resources/pirna_pathway.tsv")

de <- read_rds("results/deg/ourKD.de.grs.rds")

de <- de$adjusted |>
  map_df(as_tibble, .id="comparison")

de <- de |>
  filter(feature %in% te_regulators$gene_ID) |>
  left_join(te_regulators,by=c(feature="gene_ID")) |>
  separate(comparison,into=c("z","target","sex","tissue","driver"), remove = F,extra = "drop",sep="_") |>
  dplyr::select(-z) |>
  dplyr::select(target,sex,tissue,driver,feature,gene_symbol,padj,log2FoldChange) |>
  mutate(target=if_else(target=="NFI","NfI",target))


writexl::write_xlsx(de,"results/tables/table_known_te_regulators_in_knockdown.xlsx")