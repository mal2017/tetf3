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

# embryo only for consistency
# these are narrowPeak format with multiple summits per peak, so need to reduce
pks0 <- list(pan = c("https://www.encodeproject.org/files/ENCFF118OVG/@@download/ENCFF118OVG.bed.gz",
                     "https://www.encodeproject.org/files/ENCFF878WBZ/@@download/ENCFF878WBZ.bed.gz",
                     "https://www.encodeproject.org/files/ENCFF432EKX/@@download/ENCFF432EKX.bed.gz",
                     "https://www.encodeproject.org/files/ENCFF004BZT/@@download/ENCFF004BZT.bed.gz"),
             gro = "https://www.encodeproject.org/files/ENCFF724OZI/@@download/ENCFF724OZI.bed.gz",
            NfI = "https://www.encodeproject.org/files/ENCFF919JDR/@@download/ENCFF919JDR.bed.gz",
            vvl = "https://www.encodeproject.org/files/ENCFF249YBA/@@download/ENCFF249YBA.bed.gz",
            CG16779 = c("https://www.encodeproject.org/files/ENCFF472HIX/@@download/ENCFF472HIX.bed.gz")) |>
  map(~map(.x,import,format="narrowPeak"))

pks <- map(pks0, GRangesList) |>
  map(unlist) |>
  map(GenomicRanges::reduce)

pks <- pks |>
  GRangesList() |>
  unlist() %>%
  plyranges::mutate(.,ChIP = names(.)) 

seqlevelsStyle(pks) <- "NCBI"

expressed_in_embryo <- read_tsv("resources/embryo_expressed.FlyBase_IDs.txt",col_names = "gene_id")

# get a single 'tss' for each gene
all_genes <- genes(txdb) %>% 
  anchor_5p() %>%
  mutate(width=0)

# annotate these gene tss' by membership in set of 
# previously id'd TE regulators
all_genes <- all_genes |>
  mutate(is.piRNA.pathway = gene_id %in% pirna_gene_ids$gene_ID,
         embryo.expressed =  gene_id %in% expressed_in_embryo$gene_id)

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
  mutate(ChIP2  = ifelse(ChIP %in% c("NfI","CG16779","pan","vvl","gro"),ChIP,"other"))

write_rds(gr, snakemake@output$rds)
