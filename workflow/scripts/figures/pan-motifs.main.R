library(tidyverse)
library(plotgardener)
library(patchwork)
library(GenomeInfoDb)
library(AnnotationDbi)
library(rtracklayer)

# ------------------------------------------------------------------------------
# get info for plotting tracks
# ------------------------------------------------------------------------------
bws <- Sys.glob("upstream/csem_mosaics/bigwigs/pan_ENCSR636FRF_rep*.log2ratio.bw") |>
  c(Sys.glob("upstream/csem_mosaics/bigwigs/gro_ENCSR981DLO_rep1.log2ratio.bw")) |>
  c(Sys.glob("upstream/csem_mosaics/bigwigs/E0.4_H3K9Me3_ChIPSeq_1.log2ratio.bw"))

bws <- as.list(bws)

names(bws) <- str_extract(bws,"(?<=bigwigs\\/).+(?=\\.log2ratio)")

bws <- bws |> map(import)

bws <- map(bws, ~{seqlevelsStyle(.x) <- "UCSC";return(.x)})

yrng <- calcSignalRange(bws,negData = T) |> round(1)
chroi <- c("chr2L","chr2R","chr3L","chr3R","chrX")
sl <- getChromInfoFromNCBI("GCF_000001215.4") |>
  filter(UCSCStyleName%in% chroi)

chrend <- sl |> dplyr::select(UCSCStyleName,SequenceLength) |>
  deframe()

sl <- sl |>
  mutate(maxw = ((7.5-(5*0.05))/5)*SequenceLength/mean(SequenceLength)) |> # calc size of each plot, normalized to the mean chrlen = 1.5in
  dplyr::select(UCSCStyleName,maxw) |>
  deframe()

tracklabs <- sprintf("%s",str_extract(names(bws),"pan|gro|H3K9Me3"))
enclabs <- sprintf("%s",replace_na(str_extract(names(bws),"ENCSR.+(?=_rep)"),""))

# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + theme(text=element_text(size=unit(7,"pt"))))

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

# ------------------------------------------------------------------------------
# plot motif figs
# ------------------------------------------------------------------------------

motif_fig_df <- read_rds("results/motifs/comparison/pan_denovo_comparison.meme.gg_df.rds") |>
  mutate(g_rnk = map2(denovo,g_rnk,~{.y + labs(title=.x)}))

g_a <- plotGG(motif_fig_df$g_aln[[1]], x = 0.5, y=0.5, width = 2,height = 1.5)
plotText("A", x = 0.5, y=0.5)

g_b <- plotGG(motif_fig_df$g_aln[[2]], x = 3, y=0.5, width = 2*(23/9),height = 1.5)
plotText("B", x = 3, y=.5)


g_c <- plotGG(Reduce(`+`,motif_fig_df$g_rnk) + plot_layout(nrow=1,guides = "collect") & theme(legend.position = "right") & aes(color=class),
       x=0.75,y=2.25,width = 7,height=2)
plotText("C", x = 0.75, y=2.25)

# ------------------------------------------------------------------------------
# plot tracks
# ------------------------------------------------------------------------------
plotText("D", x = 0.5, y=4.25)

params_c <- pgParams(assembly = "dm6", default.units = "inches",range=yrng,fontsize=5)

# loop to plot multitracks for each chr
xi <- 0.5
for (chr in chroi) {
  
  params_n <- pgParams(title=chr,assembly = "dm6",chrom=chr,chromstart=0,chromend=chrend[[chr]])
  
  h <- 2.75
  yn <- 4.5
  w <- sl[[chr]]
  
  # only put tracklab on left side of overall figure
  if (chr == "chr2L") {
    lab <- tracklabs
  } else if (chr == "chrX"){
    lab <- enclabs
  } else {
    lab <- NULL
  }
  
  # add stacked signal tracks for 1 chrom, all bigwigs
  sig_track <- plotMultiSignal(data = bws, 
                          just = c("left","top"),
             params = c(params_c,params_n),
             range=yrng,
             label = lab,
             x=xi, 
             y=yn,
             width = w,
             height = h,
             start=0,end=chrend[[chr]],
             gapdistance = 0.01,
             fill = "#253494", linecolor = "#253494",negData = T,fontsize=5,cex=0.5)
  
  # add genomic distance scale and chr label on bottom
  annoGenomeLabel(plot = sig_track[[1]], 
                  params = c(params_c,params_n),
                  scale = "Mb", fontsize = 5, digits = 0,
                  start=0,
                  end=chrend[[chr]],
                  x=xi,
                  y = yn+h,
                  just = c("left","top"))
  
  # add signal min/max scale
  if (chr=="chrX") {
    plotText(sprintf("[ %s - %s ]",yrng[1],yrng[2]),
             fontsize = 9,
             x = xi+w,
             y=yn-0.2,
             just="right")
  }
  
  # highlight rough pericentromeres
  annoHighlight(
      plot = sig_track[[1]],
      chrom = chr,
      chromstart = if_else(str_detect(chr,"R$"),0,0.8*chrend[[chr]]),
      chromend = if_else(str_detect(chr,"R$"),0.2*chrend[[chr]],chrend[[chr]]),
      y=yn,height = h, just = c("left", "top"),
      default.units = "inches")

  xi <- xi + w + 0.05
}



dev.off()