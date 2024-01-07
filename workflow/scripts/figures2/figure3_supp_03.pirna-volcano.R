library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(ggtext)
library(patchwork)
library(plotgardener)
library(tidyverse)

# ------------------------------------------------------------------------------
# get de results
# ------------------------------------------------------------------------------

gene_universe <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2022_04.tsv.gz", skip = 4)

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()

pirna_genes <- read_tsv("results/resources/pirna_pathway.tsv")

res_path <- "results/deg/ourKD.de.grs.rds"

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
  mutate(label = sprintf("genes: %s\nTE silencers: %s\nTEs: %s",gene,`TE silencer`,TE))

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
  geom_text(data=labs, aes(x=-10,y=75,label=label),hjust=0, vjust=1) +
  theme(legend.position = "bottom",
        legend.justification="right",
        legend.box.spacing = unit(0,"pt"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.key.height = unit(0.2,"in"), legend.key.width = unit(0.05,"in")) +
  facet_wrap(~comparison, ncol=3)


# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=7))
)

dir.create("results/figures2/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_a, x = 0.5, y=0.5, width = 7.5,height = 10)

dev.off()

