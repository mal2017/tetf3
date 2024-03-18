library(tidyverse)

gsea_fl <- "results/signatures/ourKD_gsea.rds"
gsea_fl <- snakemake@input$rds
x <- read_rds(gsea_fl)

x <- x |> mutate(lab = str_replace_all(str_extract(comparison,"(?<=knockdown2_).+(?=_control)"), "_", " / "))

plot_bar <- function(obj) {
  mutate(obj, lab = fct_reorder(lab, padj)) |>
  ggplot(aes(-log10(padj), lab)) +
  geom_col() +
  geom_vline(xintercept = -log10(0.1),color="red",linetype="dashed") +
  ylab("RNAi / sex / sample / driver")
}

g_all <- x |>
  filter(signature_name == "all_tes") |>
  plot_bar()


g_factor.specific <- x |>
  filter(signature_name == kd) |>
  plot_bar()

g_teregs <- x |>
  filter(signature_name == "TE.regulators") |>
  plot_bar()

g_sirna <- x |>
  filter(signature_name == "siRNA") |>
  plot_bar()


res <- list(all_tes = g_all, 
            factor.specific=g_factor.specific,
            TE.regulators = g_teregs,
            siRNA=g_sirna)

write_rds(res,snakemake@output$gg_list)