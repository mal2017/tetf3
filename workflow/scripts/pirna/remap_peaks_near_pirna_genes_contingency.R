library(tidyverse)
library(plyranges)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(AnnotationDbi)
library(furrr)

deg <- ifelse(exists("snakemake"),snakemake@input[["deseq_gr"]],
              "results/analysis/deg/ourKD.de.grs.rds") %>%
                read_rds() %>%
  .$adjusted %>%
  map_df(as_tibble,.id="RNAi")  %>%
  mutate(feature.x = str_extract(RNAi,"NFI|CG16779")) %>%
  mutate(feature.x=if_else(feature.x == "NFI","NfI",feature.x)) %>%
  filter(padj < 0.1) %>%
  filter(!is.na(feature.x))

# import saved txdb
threads <- ifelse(exists("snakemake"),
               snakemake@threads,4)

plan(multisession, workers = threads)

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

# import stage-expressed gene ids
wpp_gene_ids <- ifelse(exists("snakemake"),
                         snakemake@input$wpp,
                         "resources/wpp_expressed.FlyBase_IDs.txt") %>%
  read_tsv(col_names = "gene_ID")

embryo_gene_ids <- ifelse(exists("snakemake"),
                       snakemake@input$embryo,
                       "resources/embryo_expressed.FlyBase_IDs.txt") %>%
  read_tsv(col_names = "gene_ID")

# import remap peaks as gr
remap0 <- ifelse(exists("snakemake"),
    snakemake@input$remap,
    "results/resources/remap.gr.rds") %>%
    readRDS()

remap0 <- remap0 %>%
  unlist() %>%
  mutate(.,ChIP = names(.)) 

remap <- remap0 %>%
  filter(ChIP %in% c("CG16779","NfI"))

# get first bp of all genes so I can set upstream distances on the fly
all_genes <- genes(txdb) %>% 
  anchor_5p() %>%
  mutate(width=0)

shared_seqs <- intersect(seqlevelsInUse(remap),seqlevelsInUse(all_genes))
seqlevels(all_genes, pruning.mode="coarse") <- shared_seqs

pirna_genes <- all_genes %>% filter(gene_id %in% pirna_gene_ids$gene_ID)
pirna_genes_wpp <- pirna_genes %>% filter(gene_id %in% wpp_gene_ids$gene_ID)
pirna_genes_embryo <- pirna_genes %>% filter(gene_id %in%embryo_gene_ids$gene_ID)
other_genes <- all_genes %>% filter(!gene_id %in% pirna_gene_ids$gene_ID)
wpp_genes <- other_genes %>% filter(gene_id %in% wpp_gene_ids$gene_ID)
embryo_genes <- other_genes %>% filter(gene_id %in% embryo_gene_ids$gene_ID)

bg_tbl <- list(NfI = embryo_genes, CG16779 = wpp_genes) %>%
  enframe(name="fac", value = "bg")

pirna_tbl <-  list(NfI = pirna_genes_embryo, CG16779 = pirna_genes_wpp) %>%
  enframe(name="fac", value="piRNA")

# NfI performed in embryo
# https://www.encodeproject.org/experiments/ENCSR085FSS/
# CG16779 performed in wpp
# https://www.encodeproject.org/experiments/ENCSR846PTH/

get_cont_tbl <- function(fac, p, bg, d=1000, fix="end") {
  foc <- countOverlaps(resize(p, fix=fix,width=d), fac)
  bg2 <- countOverlaps( resize(bg, fix=fix,width=d), fac)
  
  c(sum(foc > 0), sum(bg2 > 0), sum(foc == 0), sum(bg2 == 0)) %>%
    matrix(nrow=2,dimnames = list(c("piRNA","bg"),c("bound","unbound")), byrow = F)
}

# sanity check for my matrix population/dimnaming:
# c(1,2,3,4) %>% matrix(nrow=2, byrow=F, dimnames=list(c("A","B"),c("C","D")))

# sanity check for the actual number
#subsetByOverlaps(resize(pirna_genes, 1000, fix="end"), plyranges::filter(remap,ChIP=="NfI"))

d <- tibble(d =  c(100,250,500,1000))
res <- remap %>%
  split(.,.$ChIP) %>%
  as.list() %>%
  enframe(name="fac",value = "ChIP") %>%
  left_join(pirna_tbl) %>%
  left_join(bg_tbl) %>%
  cross_join(d) %>%
  mutate(contingency_table = pmap(list(ChIP,piRNA, bg,d), function(w,x,y,z) get_cont_tbl(w,x,y,d=z))) %>%
  #pull(contingency_table)
  mutate(fish.test = map(contingency_table,~broom::tidy(fisher.test(.)))) %>%
  unnest(fish.test)

# list(NfI = pirna_genes_embryo %>% add_nearest_distance(filter(remap,ChIP=="NfI")),
#      CG16779 = pirna_genes_wpp %>% add_nearest_distance(filter(remap,ChIP=="CG16779")),
#      NfI.background = embryo_genes %>% add_nearest_distance(filter(remap,ChIP=="NfI")),
#      CG16779.background = wpp_genes %>% add_nearest_distance(filter(remap,ChIP=="CG16779"))) %>%
#   imap(~mutate(.x, ChIP = .y)) %>%
#   GRangesList() %>%
#   unlist() %>%
#   as_tibble() %>%
#   #filter(is.na(distance))
#   separate(ChIP, into=c("ChIP","set"),sep="\\.") %>%
#   mutate(set = replace_na(set,"observed")) %>%
#   ggplot(aes(distance+1,fill=set)) +
#   geom_histogram() +
#   facet_wrap(~ChIP + set, scales="free")
  
saveRDS(res,snakemake@output[["rds"]])


# -------------------------------------------------------------------
# get list of piRNA pathway genes that are DE in
# -------------------------------------------------------------------

poi <- c("aub","piwi","rhi","del","AGO3","arx","tej","egg","vas")

# now find differentially expressed pirna genes nearby those peaks
pirna_genes %>%
  {x <- .; x@elementMetadata$distance <- NULL; x} %>%
  join_overlap_left(remap, maxgap=1000) %>%
  #filter(distance == 0) %>%
  as_tibble() %>%
  inner_join(deg, by=c(ChIP="feature.x", gene_id="feature")) %>% 
  left_join(pirna_gene_ids, by=c(gene_id="gene_ID")) %>%
  mutate(tissue = str_extract(RNAi,"male_gonad|female_head|female_gonad")) %>%
  dplyr::select(`ChIP/RNAi`=ChIP,tissue,gene_symbol,log2FoldChange, padj) %>%
  distinct() %>%
  arrange(`ChIP/RNAi`, tissue, -log2FoldChange) %>%
  #print(n=Inf)
  saveRDS(snakemake@output[["kd_chip_intersect_rds"]])
  

