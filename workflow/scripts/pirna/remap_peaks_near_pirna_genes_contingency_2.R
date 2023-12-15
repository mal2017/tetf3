library(tidyverse)
library(plyranges)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(AnnotationDbi)
library(nullranges)

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


# import ovary gene ids
ovary_gene_ids <- ifelse(exists("snakemake"),
                         snakemake@input$ovary,
                         "resources/ovary_expressed.FlyBase_IDs.txt") %>%
  read_tsv(col_names = "gene_ID")

# import remap peaks as gr
remap0 <- ifelse(exists("snakemake"),
    snakemake@input$remap,
    "results/resources/remap.gr.rds") %>%
    readRDS() %>%
  unlist()

seqlevels(remap0,pruning.mode="coarse") <- seqlevels(txdb)

seqinfo(remap0) <- seqinfo(txdb)

# unexpected dropping of elements unless I change seqlevels when a GRanges
# instead of GRanges list
remap0 <- split(remap0,names(remap0))

# get the promoter regions we'll be testing for 
pirna_genes <- genes(txdb) %>% promoters(2000) %>% filter(gene_id %in% pirna_gene_ids$gene_ID)
ovary_genes <- genes(txdb) %>% promoters(2000) %>% filter(gene_id %in% ovary_gene_ids$gene_ID)
other_genes <- genes(txdb) %>% promoters(2000) %>% filter(!gene_id %in% pirna_gene_ids$gene_ID)

# generate segmentation by ovary-expressed gene density
ovary_genes <- sortSeqlevels(ovary_genes)
ovary_genes <- sort(ovary_genes)
names(ovary_genes) <- NULL

set.seed(2)

L_s <- 5e4

seg_cbs <- segmentDensity(ovary_genes, n = 5, L_s = L_s, type = "cbs")

blockLength <- 1e4 # size of blocks to bootstrap
R <- 100

# https://nullranges.github.io/nullranges/articles/bootRanges.html
grl <- remap0[c("CG16779","NfI")] %>% as.list()

# get Granges for each TF we're looking at. expression
# yields a single Granges with columns iter and ChIP
# iter==0 is the original
boots <- map(grl,
             .f = function(x) {
               boots <- bootRanges(x, blockLength=blockLength, R=R,
                          seg=seg_cbs)
  
               combined <- x %>%
                 mutate(iter=0) %>%
                 bind_ranges(boots) %>%
                 plyranges::select(iter)
               return(combined)
             }) %>%
  imap(~mutate(.x, ChIP=.y)) %>%
  GRangesList() %>%
  unlist()

res <- pirna_genes %>%
  join_overlap_inner(boots) %>%
  group_by(ChIP, iter) %>%
  summarise(n_overlaps = n()) %>%
  as_tibble()
  
res %>% ggplot(aes(n_overlaps)) +
  geom_histogram(data = . %>% filter(iter!=0)) +
  geom_vline(data = . %>% filter(iter==0), aes(xintercept=n_overlaps), color="red") +
  facet_wrap(~ChIP)

# NFI performed in embryo

genes(txdb) %>%
  



split(boots,boots$ChIP) %>%
  map_df(.f= function(x) {
    split(x,x$iter) %>%
      map(~add_nearest_distance(pirna_genes,.x)) %>%
      map_df(as_tibble,.id="iter")
  },.id="ChIP") %>%
  ggplot(aes(iter,distance)) +
  facet_wrap(~ChIP) +
  geom_boxplot(outlier.shape = NA)
