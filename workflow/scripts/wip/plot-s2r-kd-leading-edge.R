library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(ggtext)
library(enrichplot)

x <- read_rds("results/signatures/s2rplus_te_gsea.rds")

x |>
  mutate(sig2 = if_else(TE.set == "all_tes", "all tes","gene-specific")) |>
  arrange(NES) |>
  group_by(sig2) |>
  mutate(rnk = row_number()) |>
  #filter(RNAi=="pan")
  ggplot(aes(rnk,NES)) +
  geom_col(width = 1.01) +
  geom_col(data=\(x)filter(x, RNAi=="pan"),fill="red", width = 1.01) +
  facet_grid(.~sig2, scales="free",space = "free") +
  ggrepel::geom_label_repel(data= \(x)filter(x, RNAi=="pan"),aes(label=RNAi))

x |> filter(TE.set!="all_tes") |> arrange(pvalue)


# kipferl
# Clark et al. 2017 (refs Handler's finding that TfIIS suppresses TEs)
x |> 
  filter(RNAi %in% c("CG2678","TfIIS","key")) |>
  #pull(gsea)
  mutate(gg = pmap(list(RNAi, TE.set, gsea),
                 .f = function(cmp,gs, obj) {
                   enrichplot::gseaplot2(obj, geneSetID=gs, title = paste(cmp, "\n", gs))
                 }
  )) |>
  pull(gg)
