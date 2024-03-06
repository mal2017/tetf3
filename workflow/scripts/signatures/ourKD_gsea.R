library(clusterProfiler)
library(tidyverse)
library(GenomicRanges)

# ------------------------------------------------------------------------------
# read in data
# ------------------------------------------------------------------------------
gene_universe <- read_tsv("resources/fbgn_fbtr_fbpp_expanded_fb_2021_04.tsv.gz", skip = 4)

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()

lms_path <- "upstream/final-models.collected-info.tsv.gz"
lms_path <-snakemake@input[["mods"]]
lms <- read_tsv(lms_path) |> filter(significant_x)

pirna_path <- "results/resources/pirna_pathway.tsv"
pirna_path <- snakemake@input$pirna

res_path <- "results/deg/ourKD.de.grs.rds"
res_path <- snakemake@input[["res"]]
res0 <- read_rds(res_path) %>%
  map_df(~map_df(.x,as_tibble,.id="comparison"),.id="adjustment")

# ------------------------------------------------------------------------------
# make gene sets
# ------------------------------------------------------------------------------

# this is nested because we will combine it with other nested tibble below
pirna <-  read_tsv(pirna_path) |>
  dplyr::select(ensembl_gene=gene_ID) |>
  distinct() |>
  mutate(gs_name="TE.regulators") |>
  dplyr::select(gs_name, ensembl_gene) |>
  mutate(signature_name = gs_name) |>
  nest(data=-signature_name) |>
  dplyr::rename(signature = data)

# tbl holds a TE signature for each host gene
t2g <- lms %>% 
  dplyr::select(gs_name=gene_symbol,ensembl_gene=feature.y) %>%
  distinct() %>%
  arrange(gs_name)

# add a signature containing all TEs
at2g <- unnest(enframe(list(all_tes=unique(t2g$ensembl_gene)),name = "gs_name", value = "ensembl_gene"),"ensembl_gene")


signatures <- bind_rows(t2g,y=at2g) |>
  mutate(signature_name = gs_name) |> 
  nest(data=-signature_name) |>
  dplyr::rename(signature = data) |>
  #filter(signature_name == "all_tes" | signature_name == "TE.regulators") |>
  bind_rows(pirna)

# ------------------------------------------------------------------------------
# set up tibble with each results ranking and paired gene sets
# ------------------------------------------------------------------------------
res1 <- res0 %>% 
  filter(adjustment %in% c('adjusted')) %>%
  left_join(lkup, by=c(feature="gene_ID")) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol),feature, gene_symbol)) %>%
  dplyr::select(comparison,adjustment,stat, feature, gene_symbol) |>
  mutate(kd = str_extract(comparison,"(?<=knockdown2_).+?(?=_)"))

# make df with a row for every comparison we want
# per kd, one row for the knocked-down gene's TE signature
# plus one for all TEs
res2 <- res1 %>%
  dplyr::select(comparison,adjustment,kd,feature,score=stat) %>%
  arrange(-score) |>
  nest(data=c(feature,score)) %>%
  mutate(data=map(data,drop_na)) |>
  cross_join(signatures) |>
  filter(kd == signature_name | signature_name == "all_tes" | signature_name == "TE.regulators")

# ------------------------------------------------------------------------------
# run gsea
# ------------------------------------------------------------------------------
# gsea func
possibly_gsea <- possibly(function(.x,.y) {
  GSEA(deframe(.x), 
       TERM2GENE = .y,
       seed=2022,
       pvalueCutoff = 1,
       minGSSize = 3, 
       eps=0)
  },
  otherwise = NULL)

set.seed(2)

res.gsea <- res2 |>
  mutate(gsea = map2(data,signature,possibly_gsea)) %>%
  mutate(gsea.tidy = map(gsea,as_tibble)) %>%
  dplyr::select(-data) %>%
  unnest(gsea.tidy)

# internal p value adjustment doesn't make sense here, because each gsea
# run is done with 1 gene set at a time , so I do adjustment
# after myself within each class of question we're asking:
# 'all_tes', individual factors', TE regulators
gsea.tbl <- res.gsea %>%
  group_by(signature_type = if_else(signature_name == "kd","factor-specific",signature_name)) |>
  mutate(padj = p.adjust(pvalue, method="BH")) |>
  ungroup() |>
  dplyr::select(-qvalue, -p.adjust) # these were made by clusterprofiler/fgsea, but don't apply as stated in comment above

# ------------------------------------------------------------------------------
# run gsea
# ------------------------------------------------------------------------------

# "results/ourKD.gsea.tbl.rds"
saveRDS(gsea.tbl, snakemake@output[["rds"]])
