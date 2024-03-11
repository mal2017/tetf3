library(tidyverse)
library(ggnewscale)
library(plotgardener)
library(patchwork)
library(GenomeInfoDb)
library(AnnotationDbi)
library(ggbio)
library(org.Dm.eg.db)
library(GenomicFeatures)
library(rtracklayer)
library(plyranges)
library(GenomicRanges)

# ------------------------------------------------------------------------------
# get info for plotting tracks
# ------------------------------------------------------------------------------
bws <- Sys.glob("upstream/csem_mosaics/bigwigs/pan*.log2ratio.bw")

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
  mutate(maxw = ((7.5-(5*0.02))/5)*SequenceLength/mean(SequenceLength)) |>
  dplyr::select(UCSCStyleName,maxw) |>
  deframe()

# ------------------------------------------------------------------------------
# colors for whole genome tracks
# ------------------------------------------------------------------------------
supporting <- jsonlite::read_json("resources/pericent_enriched_pan_chips_by_inspection.json",simplifyVector = T)

colrtbl <- tibble(Sample =  names(bws)) |>
  mutate(pericent = Sample %in% supporting) |>
  mutate(experiment=str_extract(names(bws),"ENCSR.+(?=_)")) |> 
  mutate(experiment =fct_relevel(experiment,
                                 "ENCSR058DSI","ENCSR636FRF","ENCSR033IIP","ENCSR074LKQ")) |>
  arrange(experiment) |>
  mutate(color=if_else(pericent,"red","darkgray"))
         
# reorder via the manual fct order i specifed above
stopifnot(length(bws)==nrow(colrtbl))
bws <- bws[colrtbl$Sample]

tracklabs <- sprintf("%s",str_extract(names(bws),"ENCSR.+(?=_rep)"))

# ------------------------------------------------------------------------------
# create page 1
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=7))
)

dir.create("results/figures/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

params_c <- pgParams(assembly = "dm6", default.units = "inches",range=yrng)

# loop to plot multitracks for each chr
xi <- 0.5
h <- 6.5
yn <- 0.75
for (chr in chroi) {
  
  params_n <- pgParams(title=chr,assembly = "dm6",
                       chrom=chr,chromstart=0,chromend=chrend[[chr]])
  
  w <- sl[[chr]]
  
  # only put tracklab on left side of overall figure
  if (chr == "chr2L") {
    lab <- tracklabs
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
                               fill = colrtbl$color, 
                               linecolor = colrtbl$color,
                               negData = T,
                               fontsize=2)
  
  # add genomic distance scale and chr label on bottom
  annoGenomeLabel(plot = sig_track[[names(bws)[1]]], 
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
             y=yn,
             just="right")
  }
  
  # highlight rough pericentromeres
  annoHighlight(
    plot = sig_track[[names(bws)[1]]],
    chrom = chr,
    chromstart = if_else(str_detect(chr,"R$"),0,0.8*chrend[[chr]]),
    chromend = if_else(str_detect(chr,"R$"),0.2*chrend[[chr]],chrend[[chr]]),
    y=yn,height = h, just = c("left", "top"),
    default.units = "inches")
  
  xi <- xi + w + 0.02
}

plotText("A", x = 0.5, y=0.5)

# ------------------------------------------------------------------------------
# get other plots
# ------------------------------------------------------------------------------

# how repetitiveness vs manually id'd supporting chips
g_rep <- read_rds("results/repetitiveness/repetitiveness_by_visual_pericent_inspection_status.gg.rds") +
    ylab("(IP TE/IP genomic)/(input TE/input genomic)") +
  aes(color=visual.pericentromeric.enrichment) +
  scale_color_manual(values=c(`TRUE`="red",`FALSE`="darkgray")) +
  guides(color="none")


# show qc of manually id'd supporting chips
g_qual <- read_rds("results/repetitiveness/quality_by_visual_pericent_inspection_status.gg.rds") +
  aes(color=visual.pericentromeric.enrichment) +
  scale_color_manual(values=c(`TRUE`="red",`FALSE`="darkgray")) +
  guides(color="none")


plotGG(g_rep, x = 2, y=7.75, width = 2,height = 1.5)

plotGG(g_qual, x = 4.5, y=7.75, width = 2,height = 1.5)

plotText("B", x = 2, y=7.75)
plotText("C", x = 4.5, y=7.75)

dev.off()