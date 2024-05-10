Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(writexl)

main_ps_ctrl <- read_rds("results/phylosignal/phylosignal_df.rds") |> filter(score_type == "control")

unr_ps_ctrl <- read_rds("results/ripseq/unr_ripseq_phylosignal.tbl.rds") |> filter(coef %in% c("bm","random"))

fully_filtered <- read_tsv("results/phylosignal/phylosignal_filtered_hits.tsv.gz") |>
  filter(score_type == "score" & n_tests_sig >=1) |>
  dplyr::select(TF) |>
  distinct()

te_regs <- read_tsv("results/resources/pirna_pathway.tsv")

zads <- read_tsv("results/resources/zad_genes.tsv")

te_regs_intersection <- fully_filtered |>
  filter(TF %in%te_regs$gene_symbol)

zad_intersection <- fully_filtered |>
  filter(TF %in%zads$gene_symbol)

zad_cont <- read_rds("results/phylosignal/phylosignal_df.rds") |>
  filter(metric!="Lambda" & score_type=="score") |>
  group_by(sex,TF) |>
  summarise(significance=if_else(any(padj < 0.1),"sig","n.s.")) |>
  mutate(is_zad = if_else(TF %in% zads$gene_symbol,"ZAD","other")) |>
  group_by(is_zad,significance) |>
  summarise(n=n(),.groups = "drop") |>
  arrange(desc(significance)) |>
  pivot_wider(names_from = is_zad, values_from = n)

zad_fish <- column_to_rownames(zad_cont,"significance") |>
  fisher.test() |>
  broom::tidy()

x <- list(`coexpression score phyloSignal results`=read_rds("results/phylosignal/phylosignal_df.rds") |> filter(metric!="Lambda"),
     `Unr RIP-seq phyloSignal`=read_rds("results/ripseq/unr_ripseq_phylosignal.tbl.rds"),
     `>=1 significant tests M or F`=fully_filtered,
     `controls for main tree`=main_ps_ctrl,
     `controls for Unr specific tree`=unr_ps_ctrl,
     `intersection with piRNA regulators`=te_regs_intersection,
     `intersection with ZADs`=zad_intersection,
     `ZAD enrichment test`=zad_fish)


write_xlsx(x,snakemake@output$xlsx)