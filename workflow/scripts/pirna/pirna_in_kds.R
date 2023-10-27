library(clusterProfiler)
library(tidyverse)
library(GenomicRanges)
library(tidyverse)
library(clusterProfiler)


# ------------------------------------------------------------------------------
# piRNA genes
# ------------------------------------------------------------------------------

gene_universe <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2022_04.tsv.gz", skip = 4)

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()

pirna_tbl <- ifelse(exists("snakemake"),snakemake@input$pirna,"results/resources/pirna_pathway.tsv") %>% read_tsv()

# ------------------------------------------------------------------------------
# get rankings
# ------------------------------------------------------------------------------

res_path <- "results/deg/ourKD.de.grs.rds"
res_path <- snakemake@input[["res"]]

# import reslts
res <- read_rds(res_path) %>%
  map_df(~map_df(.x,as_tibble,.id="comparison"),.id="adjustment")

# add gene symbol name for each feature
res <- res %>% 
  filter(adjustment == 'adjusted') %>%
  left_join(lkup, by=c(feature="gene_ID")) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol),feature, gene_symbol)) %>%
  dplyr::select(comparison,log2FoldChange,padj, feature, gene_symbol)

res <- mutate(res, kd = str_extract(comparison,"(?<=knockdown2_).+?(?=_)"))

res <- res %>% mutate(kd=if_else(kd == "NFI","NfI",kd))

res <- res |> left_join(pirna_tbl, by=c("feature"="gene_ID")) |> filter(str_detect(feature,"FBgn"))

res <- res |>
  mutate(across(c(in.Handler13,in.Czech13),replace_na, F)) |>
  mutate(across(c(padj),replace_na, 1))


res <- res |>
  group_by(piRNA = in.Handler13 | in.Czech13, sig=padj < 0.1, kd, comparison)  |>
  tally() |>
  ungroup() |>
  arrange(comparison) |>
  mutate(piRNA=if_else(piRNA,"piRNA","other"),sig=if_else(sig,"sig","ns"))

fisher_res <- res |>
  pivot_wider(names_from = sig,values_from = n,values_fill = 0) |>
  dplyr::relocate(sig, .before=ns) |>
  mutate(piRNA=fct_relevel(as_factor(piRNA),"piRNA")) |>
  arrange(comparison,piRNA) |>
  nest(data=c(piRNA,sig,ns)) |>
  mutate(data=map(data,column_to_rownames,"piRNA")) |>
  mutate(fish = map(data,fisher.test)) |>
  mutate(tidy.fish = map(fish,broom::tidy)) |>
  unnest(tidy.fish)

write_rds(res, snakemake@output$de_pirna_tbl)
write_rds(fisher_res, snakemake@output$de_pirna_fisher)