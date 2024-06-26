Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(ggtext)
library(patchwork)
library(plotgardener)

# ------------------------------------------------------------------------------
# get de results
# ------------------------------------------------------------------------------

gene_universe <- read_tsv("resources/fbgn_fbtr_fbpp_expanded_fb_2021_04.tsv.gz", skip = 4)

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()

pirna_genes <- read_tsv("results/resources/pirna_pathway.tsv")
pirna_genes <- read_tsv(snakemake@input$pirna)

res_path <- "results/deg/ourKD.de.grs.rds"
res_path <- snakemake@input$grs

res <- read_rds(res_path) %>%
  map_df(~map_df(.x,as_tibble,.id="comparison"),.id="adjustment")

# add gene symbol name for each feature
res <- res %>% 
  filter(adjustment == 'adjusted') %>%
  left_join(lkup, by=c(feature="gene_ID")) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol),feature, gene_symbol)) %>%
  dplyr::select(comparison,feature, gene_symbol, log2FoldChange, pvalue, padj) |>
  mutate(comparison = str_replace_all(str_extract(comparison,"(?<=knockdown2_).+(?=_control)"), "_", " / "))

res <- res |>
  mutate(feature.type = if_else(str_detect(feature,"FBgn"),"gene","TE")) |>
  mutate(feature.type = if_else(feature %in% pirna_genes$gene_ID, "TE silencer", feature.type)) |>
  #filter(feature.type %in% c("Czech 2013/Handler 2013","TE")) |>
  mutate(feature.type = fct_relevel(feature.type, c("gene","TE silencer","TE")))

labs <- res |>
  filter(padj < 0.1) |>
  count(comparison,feature.type) |>
  pivot_wider(names_from = feature.type, values_from = n,values_fill = 0) |>
  mutate(label = sprintf("DE features\ngenes: %s\nTE silencers: %s\nTEs: %s",gene,`TE silencer`,TE))


run_phyper <- function(universe,te_regs, hits) {
  te_regs <- intersect(universe,te_regs)
  M <- length(unique(universe))
  N <- length(unique(te_regs))
  n <- length(unique(hits))
  k <- length(intersect(unique(hits),unique(te_regs)))
  
  phyper(k,N,M-N,n, lower.tail = F)
}


res |>
  filter(feature.type!="TE") |>
  nest(-comparison) |>
  mutate(ph = map(data, ~{run_phyper(.x$feature,pirna_genes$gene_ID,filter(.x,padj < 0.1)$feature)})) |>
  mutate(pht = map(ph, broom::tidy)) |>
  unnest(pht) |>
  mutate(padj = p.adjust(x,method="BH"))



# ------------------------------------------------------------------------------
# exemplary all-te, piRNA volc
# ------------------------------------------------------------------------------

g_a <- res |> 
  #filter(comparison == "knockdown2_CG16779_male_gonad_aTub_control_male_gonad_aTub") |>
  arrange(-str_detect(feature,"FBgn")) |>
  arrange(feature.type, -log2FoldChange) |>
  ggplot(aes(log2FoldChange, -log10(pvalue))) +
  geom_point(size=rel(0.5),aes(color=feature.type)) +
  scale_color_manual(values = c("TE"="red","TE silencer"="blue", gene="lightgray")) +
  geom_text(data=labs, aes(x=-10,y=Inf,label=label),hjust=0, vjust=1) +
  theme(legend.position = "bottom",
        legend.justification="right",
        legend.box.spacing = unit(0,"pt"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.key.height = unit(0.2,"in"), legend.key.width = unit(0.05,"in")) +
  facet_wrap(~comparison, ncol=3,scales = "free")

write_rds(g_a,snakemake@output$gg)