library(tidyverse)


okd <- read_rds("results/deg/ourKD.de.grs.rds")

okd <- okd$adjusted |>
  map_df(as_tibble,.id="experiment") |>
  dplyr::select(experiment, feature, baseMean, lfc=log2FoldChange, padj) |>
  mutate(knockdown= str_extract(experiment,"(?<=knockdown2_).+?(?=_)")) |>
  dplyr::relocate(knockdown, feature,.after="experiment") |>
  mutate(knockdown = str_replace(knockdown,"NFI","NfI"))

s2r <- read_tsv("results/deg/s2rplus.res.tsv.gz")
s2r <- dplyr::select(s2r, knockdown=comparison, feature, baseMean=AveExpr, lfc=logFC, padj = adj.P.Val)

dat <- inner_join(s2r, okd, by=c("knockdown","feature"), suffix=c(".s2r",".whole_tissue"))

dat |>
  filter(!str_detect(feature,"FBgn")) |>
  ggplot(aes(lfc.s2r, lfc.whole_tissue)) +
  geom_point() +
  #ggdensity::geom_hdr_points() +
  ggpubr::stat_cor() +
  facet_wrap(~experiment) +
  geom_smooth(method="lm")
