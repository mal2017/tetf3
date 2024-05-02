library(SuperCell)

fl <-  ifelse(exists("snakemake"),snakemake@input$supercell, "upstream/fca_supercells.rds")
x <- readRDS(fl)

membership <- x$SC$membership |> tibble::enframe(value = "metacell",name="cell")

library(scater)
#library(scran)
#library(scuttle)

sce <- supercell_2_sce(SC = x$SC, SC.GE = x$SC.GE)
sce$label <- x$SC$lineage[colnames(sce)]

rm(x);gc()

library(tidyverse)
library(writexl)


cor_fl <- ifelse(exists("snakemake"),snakemake@input$corr, "results/calderon22/fca_reanalysis_correlations.rds") 

x <- read_rds(cor_fl)
x <- x[1,"res.spqn"] |> unnest(res.spqn)
x <- x |> filter(padj < 0.1)

o <- list(`metacell metadata`=as_tibble(colData(sce),rownames="metacell"),
          `metacell membership`=membership,
          `significant gene vs TE correlations`=x)

write_xlsx(o,snakemake@output$xlsx)
write_rds(sce,snakemake@output$sce)
