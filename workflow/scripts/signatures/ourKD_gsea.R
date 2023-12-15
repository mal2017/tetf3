library(clusterProfiler)
library(tidyverse)
library(GenomicRanges)

gene_universe <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2022_04.tsv.gz", skip = 4)

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()

#lms_path <- "upstream/final-models.collected-info.tsv.gz"
lms_path <-snakemake@input[["mods"]]

lms <- read_tsv(lms_path)
lms <- lms %>% filter(significant_x)

res_path <- "results/deg/ourKD.de.grs.rds"
res_path <- snakemake@input[["res"]]

res <- read_rds(res_path) %>%
  map_df(~map_df(.x,as_tibble,.id="comparison"),.id="adjustment")

# add gene symbol name for each feature
res <- res %>% 
  filter(adjustment == 'adjusted') %>%
  left_join(lkup, by=c(feature="gene_ID")) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol),feature, gene_symbol)) %>%
  dplyr::select(comparison,stat, feature, gene_symbol)

res <- mutate(res, kd = str_extract(comparison,"(?<=knockdown2_).+?(?=_)"))

res <- res %>% mutate(kd=if_else(kd == "NFI","NfI",kd))

#pirna_path <- "results/resources/pirna_pathway.tsv"
pirna_path <- snakemake@input$pirna

pirna <-  read_tsv(pirna_path) |>
  dplyr::select(ensembl_gene=gene_ID) |>
  distinct() |>
  mutate(gs_name="TE.regulators") |>
  dplyr::select(gs_name, ensembl_gene)


# tbl holds a TE signature for each host gene
t2g <- lms %>% 
  dplyr::select(gs_name=gene_symbol,ensembl_gene=feature.y) %>%
  distinct() %>%
  arrange(gs_name) |>
  bind_rows(pirna)





# add a signature containing all TEs
at2g <- unnest(enframe(list(all_tes=unique(t2g$ensembl_gene)),name = "gs_name", value = "ensembl_gene"),"ensembl_gene") |>
  bind_rows(t2g,y=_) |>
  mutate(signature_name = gs_name) |> nest(data=-signature_name) |>
  dplyr::rename(signature = data) |>
  filter(signature_name %in% res$kd | signature_name == "all_tes" | signature_name == "TE.regulators")

# make df with a row for every comparison we want
# per kd, one row for the knocked-down gene's TE signature
# plus one for all TEs
res <- res %>%
  dplyr::select(comparison,kd,feature,score=stat) %>%
  arrange(-score) %>% 
  nest(data=c(feature,score)) %>%
  cross_join(at2g) |>
  filter(kd == signature_name | signature_name == "all_tes" | signature_name == "TE.regulators")

# run gsea
possibly_gsea <- possibly(function(.x,.y){GSEA(deframe(.x), TERM2GENE = .y,seed=2022,pvalueCutoff = 1,minGSSize = 5, eps=0)},otherwise = NULL)

set.seed(2)

res <- res |>
  mutate(gsea = map2(data,signature,possibly_gsea)) %>%
  mutate(gsea.tidy = map(gsea,as_tibble)) %>%
  dplyr::select(-data) %>%
  unnest(gsea.tidy)

# internal p value adjustment doesn;t make sense here, because the way I ran this
# doesn't run fgsea with multiple sets at the same time, so I do adjustment
# after this myself
# additionally, we're asking 2 different questions with 'all_tes' and individual factors'
# TE target sets, so these are adjusted separately
gsea.tbl <- res %>%
  #group_by(signature_type = signature_name == "all_tes") |>
  mutate(padj = p.adjust(pvalue, method="BH")) |>
  ungroup() |>
  dplyr::select(-qvalue, -p.adjust) # these were made by clusterprofiler/fgsea, but don't apply as stated in comment above

# "results/ourKD.gsea.tbl.rds"
saveRDS(gsea.tbl, snakemake@output[["rds"]])