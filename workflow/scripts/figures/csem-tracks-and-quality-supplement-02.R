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
# repetitivess scires for pan and select others
# ------------------------------------------------------------------------------
repet <- "results/repetitiveness/chip_repetitiveness.rds"
repet <- read_rds(repet)

g_e_repetitiveness <- repet |>
  mutate(target = fct_reorder(target,estimate)) |>
  ggplot(aes(target,estimate)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ggpubr::stat_compare_means(size=rel(2),label.x.npc = "center",label.y.npc = 0.9) +
  ylab("mapped read ratio:\n(IP TE/IP genomic) / (WCE TE/WCE genomic)")

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

h3k9me3_bws <- Sys.glob("~/amarel-matt/tetf/subworkflows/tetf_csem_mosaics/results/csem_mosaics/viz/*H3K9Me3*.log2ratio.bw")
names(h3k9me3_bws) <- str_extract(h3k9me3_bws,"(?<=viz\\/).+(?=\\.log2)")

# reorder by developmental timepoint
h3k9me3_bws <- h3k9me3_bws[names(h3k9me3_bws) |> str_extract("(?<=E)[\\d\\.]+(?=_)") |> as.numeric() |> order()]

bws <- list(h3k9me3 = h3k9me3_bws)

grsl <- map(bws, ~GRangesList(map(.x, ~{import_and_summarize(.x, tiles, wh)} ))) # this func from workflow/scripts/utils/plotting.R

gs <- map(grsl, plot_genome_signal)


# ------------------------------------------------------------------------------
# create page 1
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=5))
)

dir.create("results/figures/")

pdf("results/figures/csem-tracks-h3k9me3-profile.pdf",width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(gs$h3k9me3@ggplot, x = 0.5, y=0.5, width = 7.5,height = 3)

plotGG(g_e_repetitiveness + theme(axis.text.y = element_text(size=5)), x = 0.25, y=3.75, width = 7.5,height = 2)

plotText("A", x = 0.5, y=0.5)
plotText("B", x = 0.5, y=3.75)

dev.off()
