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


o <- list(`metacell metadata`=as_tibble(colData(sce),rownames="metacell"),`metacell membership`=membership)

write_xlsx(o,snakemake@output$xlsx)
write_rds(sce,snakemake@output$sce)
