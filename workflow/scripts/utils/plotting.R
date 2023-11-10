library(ggbio)
library(plyranges)
library(GenomicRanges)
library(GenomicAlignments)
library(tidyverse)
library(AnnotationDbi)
library(GenomicFeatures)

# func to summarize by fixed length windows from a bigwig file
import_and_summarize <- function(x, tiles, which) {
  
  gr <- rtracklayer::import.bw(x,which=which)
  
  gr$score <- replace_na(gr$score, replace = 0)
  
  ba <- binnedAverage(bins=tiles,
                      numvar=mcolAsRleList(gr, varname = "score")[seqlevels(tiles)],
                      varname="score", na.rm=T)
}


# func to make plot of binned granges genome wide, uses ggbio. It returns an s3 object, one slot of which
# holds a plain ggplot2 / grob objec
plot_genome_signal <-  function(x) {
  g <- autoplot(x, aes(y=score), 
                geom="area", 
                chr.weight=c("2L"=1/6,"2R"=1/6,"3R"=1/6,"3L"=1/6,"4"=1/6,"X"=1/6)) +
    #theme_clear() +
    theme_genome() +
    facet_grid(grl_name~seqnames, scales="free",space = "free_x") +
    ylab("log2(IP/input)") +
    xlab("chromosomal position")
  
  return(g)
}

# plots the csem log2 IP/input signal for a given region
# currently hardcoded for our flagship encode experiments
# intended for relatively broad views (>=5kb)
plot_region_signal <- function(roi) {
  txdb <- loadDb("results/resources/txdb")
  
  tes <- read_bed("upstream/reference_insertions.bed") |> 
    group_by(name) |>
    reduce_ranges() |>
    ungroup() |>
    filter_by_overlaps(roi) |> sort()
  
  tes_grl <- split(tes, tes$name)
  
  # subset bw by region a bit biger than what we want to plot
  l2r_rep1 <- read_bigwig("upstream/pan_ENCSR058DSI_rep1.bw") |> filter_by_overlaps(mutate(anchor_center(roi),width=width*1.5))
  l2r_rep2 <- read_bigwig("upstream/pan_ENCSR058DSI_rep2.bw") |> filter_by_overlaps(mutate(anchor_center(roi),width=width*1.5))
  
  peaks_gr <- Sys.glob("upstream/pan_ENCSR058DSI_rep*.mosaics.bed") |>
    map(read_bed) |>
    GRangesList() |>
    unlist() |>
    reduce_ranges() |>
    filter_by_overlaps(roi)
  
  col <- "darkslategrey"
  peaks <- autoplot(peaks_gr, geom="rect",color=col,fill=col)
  repeats <- autoplot(tes,geom="rect" ,fill=col,color=col)
  signal_rep1  <- autoplot(l2r_rep1, geom="area", aes(y=score),ylab="",fill=col,color=col)
  signal_rep2  <- autoplot(l2r_rep2, geom="area", aes(y=score), ylab="",fill=col,color=col)
  genes <- autoplot(txdb, which=resize(roi,width = 5e6,fix = "center"), 
                    fill=col,color=col,
                    geom="rect", 
                    label=T, 
                    names.expr="gene_id",
                    truncate.gaps = F,
                    mode="reduce")
  
  if(length(filter_by_overlaps(genes(txdb),roi)) == 1) {
    genes <- genes + annotate(geom="text",label=filter_by_overlaps(genes(txdb),roi)$gene_id,x=start(roi),y=Inf, hjust=0, vjust=1.1)
  }
  
  x <- tracks(peaks = peaks,
              repeats = repeats,
              `CSEM signal (rep 1)`=signal_rep1,
              `CSEM signal (rep 2)`=signal_rep2,
              genes = genes,
              title = as.character(roi),
              heights = c(0.25,0.35,1,1,0.5),
              label.bg.fill = "white",
              label.bg.color = "black",
              label.width = unit(15,"pt"), label.text.cex = 0.5)
  
  xlim(x) <- roi  
  
  return(x)
}


# plots BT1 alignments from my CSEM/MOSAICS pipeline
# currently hardcoded for our flagship encode experiments
# intended for relatively narrow views (25bp-500bp)
plot_region_reads <- function(roi) {
  
  txdb <- loadDb("results/resources/txdb")
  
  tes <- read_bed("upstream/reference_insertions.bed") |> 
    group_by(name) |>
    reduce_ranges() |>
    ungroup() |>
    filter_by_overlaps(roi) |> sort()
  
  tes_grl <- split(tes, tes$name)
  
  bfq0_rep1 <- readGAlignments("upstream/pan_ENCSR058DSI_rep1.sorted.bam", param=ScanBamParam(which=roi, mapqFilter = 0), use.names = T)
  bfq3_rep1 <- readGAlignments("upstream/pan_ENCSR058DSI_rep1.sorted.bam", param=ScanBamParam(which=roi, mapqFilter = 3), use.names = T)
  
  bfq0_rep2 <- readGAlignments("upstream/pan_ENCSR058DSI_rep2.sorted.bam", param=ScanBamParam(which=roi, mapqFilter = 0), use.names = T)
  bfq3_rep2 <- readGAlignments("upstream/pan_ENCSR058DSI_rep2.sorted.bam", param=ScanBamParam(which=roi, mapqFilter = 3), use.names = T)  
  
  readsq0_rep1 <- autoplot(bfq0_rep1, which=roi, geom="rect", aes(fill=strand), color=NA)
  readsq3_rep1 <- autoplot(bfq3_rep1, which=roi, geom="rect", aes(fill=strand), color=NA, legend=F) + guides(fill=F)
  readsq0_rep2 <- autoplot(bfq0_rep2, which=roi, geom="rect", aes(fill=strand),color=NA, legend=F) + guides(fill=F)
  readsq3_rep2 <- autoplot(bfq3_rep2, which=roi, geom="rect", aes(fill=strand), color=NA, legend=F) + guides(fill=F)
  repeats <- autoplot(tes_grl, group.selfish=F)
  
  x <- tracks(`rep 1 mapq>=0`=readsq0_rep1, 
              `rep 1 mapq>=3`=readsq3_rep1,
              `rep 2 mapq>=0`=readsq0_rep2, 
              `rep 2 mapq>=3`=readsq3_rep2,
              repeats=repeats, 
              title = as.character(roi),
              heights = c(1,1,1,1,0.25),
              label.bg.fill = "white",
              label.bg.color = "black",
              label.width = unit(15,"pt"), label.text.cex = 0.5)
  
  xlim(x) <- roi
  
  return(x)
}

possibly_plot_region_reads <- possibly(plot_region_reads,otherwise = NULL)


