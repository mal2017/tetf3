library(tidyverse)
library(GenomicRanges)
library(gt)

te_regulators <- read_tsv("results/resources/pirna_pathway.tsv")

de <- read_rds("results/deg/ourKD.de.grs.rds")

de <- de$adjusted |>
  map_df(as_tibble, .id="comparison")

de <- de |>
  filter(padj < 0.1 & feature %in% te_regulators$gene_ID) |>
  left_join(te_regulators,by=c(feature="gene_ID")) |>
  separate(comparison,into=c("z","target","sex","tissue","driver"), remove = F,extra = "drop",sep="_") |>
  dplyr::select(-z) |>
  filter(driver %in% c("Mef2.R","tj","aTub")) |>
  dplyr::select(target,sex,tissue,driver,feature,gene_symbol,padj,log2FoldChange) |>
  mutate(target=if_else(target=="NFI","NfI",target))

d <- read_rds("results/pirna/encode_peaks_dist_to_pirna.gr.rds")

# Unr is rna binding
d |>
  as_tibble() |>
  filter(is.piRNA.pathway) |>
  left_join(te_regulators, by=c(gene_id="gene_ID")) |>
  dplyr::select(feature="gene_id",gene_symbol,ChIP,distance) |>
  distinct() |>
  left_join(de,y=_,by=c(feature="feature",gene_symbol="gene_symbol",target="ChIP")) |>
  filter(target!="Unr") |>
  filter(distance < 350) 


d2 <- d |>
  as_tibble() |>
  filter(is.piRNA.pathway) |>
  left_join(te_regulators, by=c(gene_id="gene_ID")) |>
  dplyr::select(feature="gene_id",gene_symbol,ChIP,distance) |>
  distinct() |>
  mutate(bound = distance < 350)


tab <- left_join(de,y=d2,by=c(feature="feature",gene_symbol="gene_symbol",target="ChIP")) |>
  filter(target!="Unr") |>
  mutate(gene_symbol = if_else(bound,paste0(gene_symbol,"<sup>†</sup>"),gene_symbol)) |>
  group_by(RNAi=target,sex,tissue,driver,direction=if_else(sign(log2FoldChange)==1,"increase","decrease")) |>
  summarise(gene=paste(gene_symbol,collapse=", "),.groups = "drop") |>
  pivot_wider(names_from = direction,values_from = gene,values_fill = "") |>
  gt() |>
  gt::tab_spanner(label="expression change",columns=c("increase","decrease")) |>
  tab_style(style = cell_text(style="italic"),
            locations = cells_body(columns=c("increase","decrease","RNAi"))) |>
  gt::tab_footnote("†: evidence that factor targeted for RNAi binds within 350bp of 5' end of gene by ENCODE ChIP") |>
  fmt_markdown(columns = everything(), rows = everything())

dir.create("results/tables/")
gtsave(tab,"results/tables/known_TE_regulators_in_knockdown.table.docx")