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

source("workflow/scripts/utils/plotting.R")

# ------------------------------------------------------------------------------
# whole genome tracks
# ------------------------------------------------------------------------------
sl <- getChromInfoFromNCBI("GCF_000001215.4")

chroi <-  c("2L","2R","3L","3R","X","4")

# get full range of main chroms
wh <-  sl |>
  as_tibble() |>
  dplyr::filter(SequenceName %in%chroi) |>
  dplyr::select(chr=SequenceName,end=SequenceLength) |>
  mutate(start=1) |>
  GRanges()

# get names/paths of bigwigs

# get tiles and average within each tile - plotting takes forever with default 50bp windows
tiles <- tileGenome(deframe(sl[,c("SequenceName","SequenceLength")])[chroi],tilewidth = 50000, cut.last.tile.in.chrom=TRUE)

pan_bws <- Sys.glob("~/amarel-matt/tetf/subworkflows/tetf_csem_mosaics/results/csem_mosaics/viz/pan*.log2ratio.bw")
names(pan_bws) <- str_extract(pan_bws,"(?<=viz\\/).+(?=\\.log2)")

# specify the set of colors I'll use later, as well as a manual ordering of pan samples
colrtbl <- tibble(experiment=unique(str_extract(names(pan_bws),"ENCSR.+(?=_)")),
                  cl=paletteer::paletteer_d("khroma::dark"))

colrtbl <- tibble(Sample =  names(pan_bws)) |>
  mutate(experiment=str_extract(names(pan_bws),"ENCSR.+(?=_)")) |> inner_join(x=colrtbl,y=_) |>
  mutate(experiment =fct_relevel(experiment,
                                 "ENCSR058DSI","ENCSR636FRF","ENCSR033IIP","ENCSR074LKQ")) |>
  arrange(experiment)
         
# reorder via the manual fct order i specifed above
pan_bws <- pan_bws[order(match(str_extract(names(pan_bws),"ENCSR.+(?=_)"), colrtbl$experiment))]

bws <- list(pan = pan_bws)

grsl <- map(bws, ~GRangesList(map(.x, ~{import_and_summarize(.x, tiles, wh)} ))) # this func from workflow/scripts/utils/plotting.R

gs <- map(grsl, plot_genome_signal)

g_pan_profile <- gs$pan@ggplot + aes(color=grl_name, fill=grl_name) + 
  scale_color_manual(values = deframe(mutate(colrtbl, cl=as.character(cl)) |> dplyr::select(Sample,cl))) +
  scale_fill_manual(values = deframe(mutate(colrtbl, cl=as.character(cl)) |> dplyr::select(Sample,cl))) +
  guides(fill="none",color="none") +
  theme(strip.text.y=element_blank())

# ------------------------------------------------------------------------------
# qc for pan chips
# ------------------------------------------------------------------------------

qc_df0 <- Sys.glob("~/amarel-matt/tetf/subworkflows/tetf_basic_chip/results/basic_chip/qc/masked/pan*_rep*.fingerprint.metrics.txt") |>
  map_df(read_tsv)

qc_df <- filter(qc_df0, !str_detect(Sample,"input"))

qc_df <- qc_df |> mutate(experiment = str_extract(Sample,"ENCSR.+(?=_rep)"))

qc_df <- mutate(qc_df, library = str_extract(Sample,"pan_.+_rep\\d"))
  
qc_df <- qc_df |>
  pivot_longer(-c(Sample,experiment, library), names_to = "metric",values_to = "score") |>
  filter(metric == "JS Distance") 

qc_df <- qc_df |> mutate(experiment=factor(experiment, levels=levels(colrtbl$experiment)))

qc_df <- qc_df |> left_join(x=colrtbl, y=_, by=c("Sample"="library",experiment="experiment"))

qc_df <- qc_df |> arrange(experiment, Sample) |> mutate(Sample= factor(Sample))

g_js <- qc_df |> mutate(Sample=fct_reorder(Sample,-as.numeric(experiment))) |>
  arrange(Sample) |>
ggplot(aes(score,Sample, fill=as.character(`cl`), color=as.character(`cl`))) +
    geom_point(size=rel(2)) +
  geom_path(aes(group=experiment)) +
  #geom_col() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
   xlab("JS Distance") +
  scale_color_identity() +
  scale_fill_identity() +
  xlim(c(0,1)) +
  scale_x_sqrt()

# ------------------------------------------------------------------------------
# create page 1
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=5))
)

dir.create("results/figures2/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_pan_profile, x = 3, y=0.5, width = 4.5,height = 5)

plotGG(g_js, x = 0.25, y=0.5, width = 2.75,height = 5.2)

plotText("A", x = 0.5, y=0.5)
plotText("B", x = 3, y=0.5)

dev.off()