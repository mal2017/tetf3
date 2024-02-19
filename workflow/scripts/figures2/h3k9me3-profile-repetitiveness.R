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
library(ggsignif)
source("workflow/scripts/utils/plotting.R")


# ------------------------------------------------------------------------------
# repetitivess scires for pan and select others
# ------------------------------------------------------------------------------
repet <- "results/repetitiveness/chip_repetitiveness.rds"
repet <- read_rds(repet) |>
  filter(str_detect(target,"H3K|pan|CG16779|vvl|NfI")) |>
  mutate(target = fct_reorder(target,ratio.te, .fun=mean))

pw <- \(x) {crossing(x,set_names(x,paste0(colnames(x),".2")),.name_repair = "universal")}

cmps <- crossing(x=unique(repet$target),y=unique(repet$target)) |>
  filter(x == "pan" & y!= "pan") |>
  mutate(l = map2(x,y,~as.character(c(.x,.y)))) |>
  pull(l)

# check normality for t test
repet_test <- repet |>
  dplyr::select(target,ratio.te) |>
  group_by(target) |>
  summarise(data=list(ratio.te)) |>
  #filter(map_dbl(data,length) > 2) |>
  #mutate(normality = map(data, ~broom::tidy(shapiro.test(.x)))) |>
  #unnest(normality,names_sep = "_") |>
  #filter(normality_p.value>0.05) |>
  pw() |>
  filter(target == "pan" &target.2!="pan") |>
  mutate(ct = map2(data,data.2,~broom::tidy(t.test(.x,.y)))) |>
  unnest(ct) |>
  arrange(p.value) |>
  ungroup() |>
  mutate(padj = p.adjust(p.value,method="BH")) |>
  mutate(cmp=map2(target,target.2,.f=~(c(.x,.y))))

g_e_repetitiveness <- repet |>
  ggplot(aes(target,ratio.te)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ggsignif::geom_signif(comparisons = map(repet_test$cmp,as.character),annotations = paste0("BH padj=",format(repet_test$padj,scientific = T,digits = 2)),step_increase = 0.5) +
  #ggpubr::stat_compare_means(size=rel(2),label.x.npc = "center",label.y.npc = 0.9,ref.group = "pan") +
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

h3k9me3_bws <- Sys.glob("upstream/csem_mosaics/bigwigs/*H3K9Me3*.log2ratio.bw")
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

dir.create("results/figures2/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(gs$h3k9me3@ggplot, x = 0.5, y=0.5, width = 7.5,height = 3)

plotGG(g_e_repetitiveness + theme(axis.text.y = element_text(size=5)), x = 0.5, y=3.75, width = 4.1,height = 3.5)

plotText("A", x = 0.5, y=0.5)
plotText("B", x = 0.5, y=3.75)

dev.off()
