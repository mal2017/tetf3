library(tidyverse)
library(tximport)
library(DESeq2)
library(plyranges)

json_fl <- "results/resources/coexpressed_tes.json"
json <-snakemake@input$te_json
coex <- jsonlite::read_json(json_fl)

alpha <- 0.1
alpha <- snakemake@params$alpha

te_names <- Biostrings::readDNAStringSet("~/work/tetf3/resources/Tidalbase_transposon_sequence.fasta") |>
  names()

gtf <- snakemake@input$gtf
gtf <-"upstream/transcripts-plus-tes.gtf"
txome <- plyranges::read_gff(gtf)

allowed.types <- c("mRNA")

allowed.features <- unique(filter(txome,type %in% allowed.types)$gene_id)

fl <- snakemake@input$tsv
fl <- "upstream/unr_ripseq_counts.tsv.gz"

star.dat <- read_tsv(fl,col_names = c("feature","unstranded","first.strand","second.strand","sample")) |>
    filter(!str_detect(feature,"^N_"))

# confirm strand

#star.dat |>
#  filter(str_detect(feature,"FBgn")) |>
# pivot_longer(-c(sample,feature,unstranded),names_to = "library.mode", values_to = "count") |>
#  ggplot(aes(sample,count+1,fill=library.mode)) +
#  geom_boxplot() +
#  scale_y_log10() +
#  theme(axis.text.x=element_text(angle=45,hjust=1))

# make into matrix
star.mat <- star.dat |> 
  dplyr::select(-unstranded,-first.strand) |>
  pivot_wider(names_from = "sample",values_from = "second.strand") |>
  filter(feature %in% allowed.features) |>
  column_to_rownames("feature") |>
  as.matrix()


# get samples as DF
coldat <- colnames(star.mat) |>
  tibble(sample = _) |>
  mutate(condition= str_extract(sample,"WT|CTRL|274|DM")) |>
  mutate(possible.batch = paste0("grp",str_extract(sample,"47|48|49"))) |>
  column_to_rownames("sample") |>
  DataFrame()

dds <- DESeqDataSetFromMatrix(star.mat,coldat, design = ~condition)


smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$condition <- relevel(dds$condition,"WT")

dds <- DESeq(dds)

res <- list(SM_vs_WT=c("condition","274","WT"),
            DM_vs_WT=c("condition","DM","WT"),
            SM_vs_CTRL=c("condition","274","CTRL"),
            DM_vs_CTRL=c("condition","DM","CTRL"),
            WT_vs_CTRL=c("condition","WT","CTRL")) |>
  enframe(name="comparison",value = "contrast") |>
  mutate(deseq.res = map(contrast,~as_tibble(results(dds,contrast = .x, alpha=alpha,altHypothesis="greater"),rownames="feature")))

res <- res |>
  dplyr::select(comparison,deseq.res) |>
  unnest(deseq.res)


res <- txome |>
  as_tibble() |>
  dplyr::select(gene_id,gene_symbol) |>
  distinct() |>
  left_join(res,y=_, by=c("feature"="gene_id")) |>
  dplyr::relocate(gene_symbol,.before="feature")

res <- res |> 
  mutate(feature.type = if_else(feature %in% te_names,"TE","gene"))

# check vs figure rom hollmann et al - 
# note I'm using an ma plot here bc I do a 1 way test so the p-values
# don't make sense to plot as a volcano (because I'm interested in
# Unr-bound transcripts, but the genes they call out as outlers are also outliers here)
#fromHollmann2020 <- c("Dro","CecB","CecA","CecA1","pain","pirk","Mtk","PGRP-LB")
#res |>
#  filter(comparison %in% c("DM_vs_WT","SM_vs_WT")) |>
#  ggplot(aes(baseMean,log2FoldChange)) +
# facet_wrap(~comparison) +
# #geom_point()+
#  geom_point(data=\(x)filter(x,str_detect(feature,"FBgn")),size=rel(0.5)) +
#  geom_point(data=\(x)filter(x,gene_symbol %in% fromHollmann2020),color="red") +
# geom_text(data=\(x)filter(x,gene_symbol %in% fromHollmann2020),aes(label=gene_symbol),color="red") +
# scale_x_log10()

# final calls
unr_rip <- res |>
  dplyr::select(comparison,feature=gene_symbol,rip_baseMean=baseMean, rip_log2FoldChange=log2FoldChange,rip_padj=padj,rip_padj=padj, rip_stat=stat) |>
  filter(comparison == "WT_vs_CTRL") |>
  mutate(type = if_else(feature %in% te_names,"TE","gene")) |>
  mutate(status = if_else(rip_padj < alpha & rip_log2FoldChange > 0,"bound","not_bound")) |>
  mutate(coex = feature %in% coex$all$Unr) |>
  mutate(status = replace_na(status,"not_bound")) |>
  mutate(group=case_when(type == "gene" & status == "bound" ~ "bound gene",
                       type == "gene" & status == "not_bound" ~ "non-bound gene",
                       type == "TE" & status == "bound" & coex ~ "coexpressed and bound TE",
                       type == "TE" & status == "bound" & !coex ~ "noncoex. and bound TE",
                       type == "TE" & status == "not_bound" & coex ~ "coexpressed and non-bound TE",
                       type == "TE" & status == "not_bound" & !coex ~ "noncoex. and non-bound TE"))

write_tsv(unr_rip,snakemake@output$tsv)
write_rds(dds,snakemake@output$dds)
