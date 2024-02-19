library(tidyverse)
library(patchwork)
library(plotgardener)

gsea_df <- read_rds('results/signatures/s2rplus_te_gsea.rds')
lkup_path <- "results/resources/gene_symbol_lookup.tsv.gz"

lkup <- read_tsv(lkup_path)

# ------------------------------------------------------------------------------
# barplot
# ------------------------------------------------------------------------------

g_b <- gsea_df |>
  dplyr::select(-gsea) |>
  filter(TE.set == "all_tes") |>
  mutate(rnk = dense_rank(NES)) |>
  ggplot(aes(rnk,NES)) +
  geom_col(width=1.01) +
  geom_col(data =  . %>% filter(RNAi %in% c("pan","NfI")), width=1.01, color="red") +
  geom_text(data = . %>% filter(RNAi %in% c("pan","NfI")), 
            aes(label=RNAi, y=NES + sign(NES)*0.3),nudge_x = 10,
            color="red", fontface="italic") +
  annotate(label=paste0("n=",length(unique(gsea_df$RNAi))),x=Inf,y=-Inf, geom="text", hjust=1, vjust=0) +
  xlab("rank") + ylab("NES (all TEs)")


g_cs <- gsea_df |>
  filter(RNAi %in% c("pan","NfI") & TE.set == "all_tes") |>
  mutate(lab = paste0(RNAi,"; NES=", round(NES,2),"; BH padj=",format.pval(padj))) |>
  mutate(plt = map2(gsea,lab,~enrichplot::gseaplot2(.x,geneSetID="all_tes", title=.y))) |>
  pull(plt) |>
  map( ~{ .x & theme(axis.title = element_text(size=5), axis.text = element_text(size=5), plot.title = element_text(size=7, hjust=0.5))})


# ------------------------------------------------------------------------------
# All TE leading edge bar chart
# ------------------------------------------------------------------------------
x <- read_rds("results/signatures/ourKD_gsea.rds")

x <- x |> mutate(lab = str_replace_all(str_extract(comparison,"(?<=knockdown2_).+(?=_control)"), "_", " / "))

plot_barchart <- function(dat) {
  dat |>
    filter(signature_name == kd) |>
    mutate(lab = fct_reorder(lab, padj)) |>
    ggplot(aes(-log10(padj), lab)) +
    geom_col() +
    geom_vline(xintercept = -log10(0.1),color="red",linetype="dashed") +
    ylab("RNAi / sex / sample / driver")
}

g_a <- plot_barchart(x)



# ------------------------------------------------------------------------------
# create page

theme_set(theme_classic() + 
            theme(text = element_text(size=7)) +
            theme(plot.title = element_text(hjust = 0.5))
)

dir.create("results/figures2/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_a, x=0.5, y=0.5, width=3.75, height=2.5)
plotText("A", x = 0.5, y=0.5)


plotGG(g_b, x = 4.25, y=0.7, width = 3.85,height = 2)
plotText("B", x = 4.5, y=0.5)


plotGG(g_cs[[1]], x = .5, y=3.5, width = 3.7,height = 3)
plotText("C", x = .5, y=3.5)

plotGG(g_cs[[2]], x = 4.5, y=3.5, width = 3.7,height = 3)
plotText("D", x = 4.5, y=3.5)

dev.off()


