library(tidyverse)
library(plyranges)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(AnnotationDbi)


# import saved txdb
txdb <- ifelse(exists("snakemake"),
    snakemake@input$txdb,
    "results/resources/txdb") %>%
    loadDb()

# import pirna gene ids
pirna_gene_ids <- ifelse(exists("snakemake"),
    snakemake@input$pirna,
    "results/resources/pirna_pathway.tsv") %>%
    read_tsv()

# import remap peaks as gr -  these are my 'simplified peaks',
# meaning that peaks are merged when multiple experiments exist in 
# remap22
remap <- ifelse(exists("snakemake"),
    snakemake@input$remap,
    "results/resources/remap.gr.rds") %>%
    readRDS()

# get the remap peaks we're interested in
remap <- remap[names(remap) %in% c("pan","CG16779","NfI","vvl","Unr")]

remap <- remap %>%
  unlist() %>%
  mutate(.,ChIP = names(.)) 

# get a single 'tss' for each gene
all_genes <- genes(txdb) %>% 
  anchor_5p() %>%
  mutate(width=0)

# annotate these gene tss' by membership in set of 
# previously id'd TE regulators
all_genes <- all_genes |>
  mutate(is.piRNA.pathway = gene_id %in% pirna_gene_ids$gene_ID)

# make sure all seqlevels are shared between ranges to avoid warnings/errors
shared_seqs <- intersect(seqlevelsInUse(remap),seqlevelsInUse(all_genes))
seqlevels(all_genes, pruning.mode="coarse") <- shared_seqs
seqlevels(remap, pruning.mode="coarse") <- shared_seqs

# split peaks by ChIP, then annotate all_genes with the nearest
gr <- remap %>% 
  split(.,.$ChIP) %>%
  as.list() %>%
  map(~join_nearest(all_genes, .x, suffix = c(".piRNA",".ChIP"), distance = T)) %>%
  GRangesList() %>%
  unlist() |>
  mutate(ChIP2  = ifelse(ChIP %in% c("NfI","CG16779","pan"),ChIP,"other"))


write_rds(gr, snakemake@output$rds)

