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

# these are narrowPeak format with multiple summits per peak, so need to reduce
pks0 <- list(pan = "https://www.encodeproject.org/files/ENCFF118OVG/@@download/ENCFF118OVG.bed.gz",
            NfI = "https://www.encodeproject.org/files/ENCFF919JDR/@@download/ENCFF919JDR.bed.gz",
            vvl = "https://www.encodeproject.org/files/ENCFF741HLT/@@download/ENCFF741HLT.bed.gz",
            CG16779 = "https://www.encodeproject.org/files/ENCFF789WVF/@@download/ENCFF789WVF.bed.gz") |>
  map(import,format="narrowPeak") |>
  map(GenomicRanges::reduce)

pks <- pks0 |>
  GRangesList() |>
  unlist() %>%
  plyranges::mutate(.,ChIP = names(.)) 

seqlevelsStyle(pks) <- "NCBI"


# get a single 'tss' for each gene
all_genes <- genes(txdb) %>% 
  anchor_5p() %>%
  mutate(width=0)

# annotate these gene tss' by membership in set of 
# previously id'd TE regulators
all_genes <- all_genes |>
  mutate(is.piRNA.pathway = gene_id %in% pirna_gene_ids$gene_ID)

# make sure all seqlevels are shared between ranges to avoid warnings/errors
shared_seqs <- intersect(seqlevelsInUse(pks),seqlevelsInUse(all_genes))
seqlevels(all_genes, pruning.mode="coarse") <- shared_seqs
seqlevels(pks, pruning.mode="coarse") <- shared_seqs

# split peaks by ChIP, then annotate all_genes with the nearest
gr <- pks %>% 
  split(.,.$ChIP) %>%
  as.list() %>%
  map(~join_nearest(all_genes, .x, suffix = c(".piRNA",".ChIP"), distance = T)) %>%
  GRangesList() %>%
  unlist() |>
  mutate(ChIP2  = ifelse(ChIP %in% c("NfI","CG16779","pan","vvl"),ChIP,"other"))

write_rds(gr, snakemake@output$rds)
