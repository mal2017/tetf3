library(plotgardener)
library(rtracklayer)
library(tidyverse)


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

pan_bws <- Sys.glob("upstream/csem_mosaics/bigwigs/pan*.log2ratio.bw")
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
#
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

bws <- Sys.glob("upstream/csem_mosaics/bigwigs/pan_ENCSR058DSI_rep*.log2ratio.bw") |>
  c(Sys.glob("upstream/csem_mosaics/bigwigs/gro_ENCSR981DLO_rep1.log2ratio.bw")) |>
  c(Sys.glob("upstream/csem_mosaics/bigwigs/E0.4_H3K9Me3_ChIPSeq_1.log2ratio.bw"))

names(bws) <- str_extract(bws,"(?<=viz\\/).+(?=\\.log2)")

# get tiles and average within each tile - plotting takes forever with default 50bp windows
tiles <- tileGenome(deframe(sl[,c("SequenceName","SequenceLength")])[chroi],tilewidth = 50000, cut.last.tile.in.chrom=TRUE)

grs <- map(bws, import_and_summarize, tiles, wh) # this func from workflow/scripts/utils/plotting.R

gs <- grs |> GRangesList() |> plot_genome_signal()

g_tracks <- gs@ggplot 
