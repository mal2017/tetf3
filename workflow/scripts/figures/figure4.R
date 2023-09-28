library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(ggtext)
library(patchwork)
library(plotgardener)
library(tidyverse)

# ------------------------------------------------------------------------------
# get signature enrichment results from our kd
# ------------------------------------------------------------------------------
x <- read_rds("results/signatures/ourKD_gsea.rds")

x <- x |> mutate(lab = str_replace_all(str_extract(comparison,"(?<=knockdown2_).+(?=_control)"), "_", " / "))
  

gg <- x |> filter(padj < 0.1) |>
  mutate(gg = pmap(list(lab, ID, gsea),
                   .f = function(lab,gs, obj) {
                     enrichplot::gseaplot2(obj, geneSetID=gs, title = lab)
                   }
  )) |>
  dplyr::select(comparison, signature_name, gg)


# ------------------------------------------------------------------------------
# get de results
# ------------------------------------------------------------------------------

gene_universe <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2022_04.tsv.gz", skip = 4)

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()

res_path <- "results/deg/ourKD.de.grs.rds"

res <- read_rds(res_path) %>%
  map_df(~map_df(.x,as_tibble,.id="comparison"),.id="adjustment")

# add gene symbol name for each feature
res <- res %>% 
  filter(adjustment == 'adjusted') %>%
  left_join(lkup, by=c(feature="gene_ID")) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol),feature, gene_symbol)) %>%
  dplyr::select(comparison,feature, gene_symbol, log2FoldChange, pvalue, padj)

# ------------------------------------------------------------------------------
# significant all-te hit barchart
# ------------------------------------------------------------------------------

g_a <- x |>
  filter(signature_name == "all_tes") |>
  mutate(lab = fct_reorder(lab, pvalue)) |>
  ggplot(aes(-log10(pvalue), lab)) +
  geom_col() +
  geom_vline(xintercept = -log10(0.05),color="red",linetype="dashed") +
  ylab("RNAi / sex / sample / driver") +
  ggtitle("all TE enrichment")


# ------------------------------------------------------------------------------
# exemplary all-te random walk
# ------------------------------------------------------------------------------


g_b <- gg |>
  filter(comparison == "knockdown2_CG16779_male_gonad_aTub_control_male_gonad_aTub" &
           signature_name == "all_tes") |>
  pull(gg) |>
  pluck(1) & theme(axis.title = element_text(size=5), axis.text = element_text(size=5), plot.title = element_text(size=7, hjust=0.5))

# ------------------------------------------------------------------------------
# exemplary all-te volc
# ------------------------------------------------------------------------------

pirna_genes <- read_tsv("results/resources/pirna_pathway.tsv")

g_c <- res |> 
  filter(comparison == "knockdown2_CG16779_male_gonad_aTub_control_male_gonad_aTub") |>
  arrange(-str_detect(feature,"FBgn")) |>
  mutate(feature.type = if_else(str_detect(feature,"FBgn"),"gene","TE")) |>
  mutate(feature.type = if_else(feature %in% pirna_genes$gene_ID, "Czech 2013/Handler 2013", feature.type)) |>
  #filter(feature.type %in% c("Czech 2013/Handler 2013","TE")) |>
  mutate(feature.type = fct_relevel(feature.type, c("gene","Czech 2013/Handler 2013","TE"))) |>
  arrange(feature.type, -log2FoldChange) |>
  ggplot(aes(log2FoldChange, -log10(pvalue), color=feature.type)) +
  geom_point(size=rel(0.5)) +
  scale_color_manual(values = c("TE"="red","Czech 2013/Handler 2013"="blue", gene="lightgray")) +
  theme(legend.position = "bottom",
        legend.justification="right",
        legend.box.spacing = unit(0,"pt"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.key.height = unit(0.2,"in"), legend.key.width = unit(0.05,"in"))

# ------------------------------------------------------------------------------
# target-specific all-te hit barchart
# ------------------------------------------------------------------------------

g_d <- x |>
  filter(signature_name != "all_tes") |>
  mutate(lab = fct_reorder(lab, pvalue)) |>
  ggplot(aes(-log10(pvalue), lab)) +
  geom_col() +
  geom_vline(xintercept = -log10(0.05),color="red",linetype="dashed") +
  ylab("RNAi / sex / sample / driver") +
  ggtitle("target-specific coexpressed TE signature")


# ------------------------------------------------------------------------------
# exemplary target-specific random walk
# ------------------------------------------------------------------------------

g_e <- gg |>
  filter(comparison == "knockdown2_pan_female_head_Mef2.R_control_female_head_Mef2.R" &
           signature_name != "all_tes") |>
  pull(gg) |>
  pluck(1) & theme(axis.title = element_text(size=5), axis.text = element_text(size=5), plot.title = element_text(size=7, hjust=0.5))



# ------------------------------------------------------------------------------
# exemplary volc #2
# ------------------------------------------------------------------------------

g_f <- res |> 
  filter(comparison == "knockdown2_pan_female_head_Mef2.R_control_female_head_Mef2.R") |>
  arrange(-str_detect(feature,"FBgn")) |>
  mutate(feature.type = if_else(str_detect(feature,"FBgn"),"gene","TE")) |>
  mutate(feature.type = if_else(feature %in% pirna_genes$gene_ID, "Czech 2013/Handler 2013", feature.type)) |>
  #filter(feature.type %in% c("Czech 2013/Handler 2013","TE")) |>
  mutate(feature.type = fct_relevel(feature.type, c("gene","Czech 2013/Handler 2013","TE"))) |>
  arrange(feature.type, -log2FoldChange) |>
  ggplot(aes(log2FoldChange, -log10(pvalue), color=feature.type)) +
  geom_point(size=rel(0.5)) +
  scale_color_manual(values = c("TE"="red","Czech 2013/Handler 2013"="blue", gene="lightgray")) +
  theme(legend.position = "bottom",
        legend.justification="right",
        legend.box.spacing = unit(0,"pt"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.key.height = unit(0.2,"in"), legend.key.width = unit(0.05,"in"))


# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=5))
)

dir.create("results/figures/")

pdf("results/figures/figure4.pdf",width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_a, x = 0.5, y=0.5, width = 2.5,height = 2.1)
plotText("A", x = 0.5, y=0.5)

plotGG(g_b, x = 3, y=0.5, width = 3,height = 2.2)
plotText("B", x = 3.5, y=0.5)

plotText("C",  x = 5.95, y=0.5)
plotGG(g_c, x = 6, y=0.5, width = 2,height = 2.1)

plotGG(g_d, x = 0.5, y=2.85, width = 2.5,height = 2.1)
plotText("D",  x = 0.5, y=2.85)


plotGG(g_e, x = 3, y=2.85, width = 3,height = 2.2)
plotText("E" , x = 3.25, y=2.85)

plotGG(g_f, x = 6, y=2.85, width = 2,height = 2.1)
plotText("F", x = 6, y=2.85)

dev.off()

