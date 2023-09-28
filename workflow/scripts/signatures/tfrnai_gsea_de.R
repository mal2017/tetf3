library(clusterProfiler)
library(tidyverse)

#coex_path <- "upstream/final-models.collected-info.tsv.gz"
#res_path <- "results/deg/s2rplus.res.tsv.gz"

coex_path <- snakemake@input[["coex"]]
res_path <- snakemake@input[["deg"]]

lms <- read_tsv(coex_path) %>% filter(significant_x)
res <- read_tsv(res_path)

t2g <- lms %>% 
  dplyr::select(gs_name=gene_symbol,ensembl_gene=feature.y) %>%
  distinct() %>%
  group_by(gs_name) %>%
  ungroup() %>%
  arrange(gs_name)

# add a signature containing all TEs
at2g <- unnest(enframe(list(all_tes=unique(t2g$ensembl_gene)),name = "gs_name", value = "ensembl_gene"),"ensembl_gene") |>
  bind_rows(t2g,y=_) |>
  mutate(signature_name = gs_name) |> nest(data=-signature_name) |>
  dplyr::rename(signature = data) |>
  filter(signature_name %in% res$comparison | signature_name == "all_tes")

# make df with a row for every comparison we want
# per kd, one row for the knocked-down gene's TE signature
# plus one for all TEs
res <- res %>%
  dplyr::select(comparison,feature,score=t) %>%
  arrange(-score) %>% 
  nest(data=c(feature,score)) %>%
  cross_join(at2g) |>
  filter(comparison == signature_name | signature_name == "all_tes")

# perform gsea
possibly_gsea <-  possibly(function(.x,.y) {
  set.seed(2022)
  GSEA(deframe(.x), TERM2GENE = .y, seed=2022,pvalueCutoff = 2, minGSSize = 5, eps = 0)
}, NULL)

res <- res |>
  mutate(gsea = map2(data,signature,possibly_gsea)) %>%
  mutate(gsea.tidy = map(gsea,as_tibble))

# we have many different DE results, but the p adjustment that occurs internally
# only would work if we're dealing with multiple gene sets per DE,
# so I recompute the rankings and readjust the p values in the broader context of 
# TF-specific rankings (one family of test)
# or separately in the context of all TE regulation
gsea.res <- res %>%
  dplyr::select(comparison,gsea,gsea.tidy) %>%
  mutate(kd = comparison) %>%
  unnest(gsea.tidy) %>%
  group_by(comparison) %>%
  mutate(pvRnk = dense_rank(-log10(pvalue)),
         nesRnk = dense_rank(abs(NES))) %>%
  ungroup() %>%
  relocate(pvRnk,nesRnk) %>%
  mutate(padj = p.adjust(pvalue, method="BH")) %>%
  ungroup() |>
  dplyr::select(RNAi=comparison,TE.set=ID,setSize,NES,pvalue, padj, gsea)

write_rds(gsea.res,snakemake@output[["rds"]])

