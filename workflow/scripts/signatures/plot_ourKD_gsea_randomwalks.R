library(tidyverse)
library(clusterProfiler)
library(enrichplot)

gsea_fl <- "results/signatures/ourKD_gsea.rds"
gsea_fl <- snakemake@input$rds
x <- read_rds(gsea_fl)

x <- x |> mutate(lab = str_replace_all(str_extract(comparison,"(?<=knockdown2_).+(?=_control)"), "_", " / "))

gg <- x |> 
  filter(padj < 0.1) |>
  mutate(gg = pmap(list(lab, ID, gsea),
                   .f = function(lab,gs, obj) {
                     enrichplot::gseaplot2(obj, geneSetID=gs, title = lab)
                   }
  )) |>
  dplyr::select(comparison, signature_name, gg)

write_rds(gg,snakemake@output$gg_df)
