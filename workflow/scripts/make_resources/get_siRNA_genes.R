library(tidyverse)
library(AnnotationHub)


fl <- "resources/Participating_Molecules_R-DME-211000.tsv"
fl <- snakemake@input$reactome_sirna
uniprot <-read_tsv(fl, col_select = 1)
ah = AnnotationHub()
library(org.Dm.eg.db)
orgdb <- query(ah, c("OrgDb","dm6"))

keytypes(org.Dm.eg.db)

# format just like my txt file with putative TE silencers
# and remove histone genes and RnaPII genes, as these seem highly non-specific
# to siRNA pathway
res <- select(org.Dm.eg.db, 
       uniprot$Identifier,  c("SYMBOL", "GENENAME","FLYBASE"),"UNIPROT") |>
  filter(!str_detect(SYMBOL,"^His|^Rpb|^RpII")) |>
  dplyr::select(gene_ID = FLYBASE,gene_symbol=SYMBOL,gene_symbol.Reactome=UNIPROT) |>
  mutate(`in.R-DME-211000`=T)

write_tsv(res, snakemake@output$tsv)
