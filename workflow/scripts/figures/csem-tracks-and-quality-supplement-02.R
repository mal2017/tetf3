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
  filter(str_detect(target,"H3K|pan")) |>
  mutate(target = fct_reorder(target,estimate))

pw <- \(x) {crossing(x,set_names(x,paste0(colnames(x),".2")),.name_repair = "universal")}

# performed t tests for pairs of normally distributed
# values
repet_test <- repet |>
  dplyr::select(target,ratio.te) |>
  group_by(target) |>
  summarise(data=list(ratio.te)) |>
  filter(map_dbl(data,length) > 2) |>
  mutate(normality = map(data, ~broom::tidy(shapiro.test(.x)))) |>
  unnest(normality,names_sep = "_") |>
  filter(normality_p.value>0.05) |>
  pw() |>
  filter(target == "pan" &target.2!="pan") |>
  mutate(ct = map2(data,data.2,~broom::tidy(t.test(.x,.y)))) |>
  unnest(ct) |>
  arrange(p.value) |>
  mutate(padj = p.adjust(p.value,method="BH")) |>
  mutate(cmp=map2(target,target.2,.f=~(c(.x,.y)))) |>
  filter(padj < 0.1)

g_e_repetitiveness <- repet |>
  ggplot(aes(target,estimate)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ggsignif::geom_signif(comparisons = map(repet_test$cmp,as.character),annotations = format(repet_test$padj,scientific = T,digits = 2),step_increase = 0.5) +
  #ggpubr::stat_compare_means(size=rel(2),label.x.npc = "center",label.y.npc = 0.9,ref.group = "pan") +
  ylab("mapped read ratio:\n(IP TE/IP genomic) / (WCE TE/WCE genomic)")


# ------------------------------------------------------------------------------
# quality vs repetitiveness
# ------------------------------------------------------------------------------

qc_df0 <- Sys.glob("~/amarel-matt/tetf/subworkflows/tetf_basic_chip/results/basic_chip/qc/masked/pan*_rep*.fingerprint.metrics.txt") |>
  map_df(read_tsv)

qc_df <- filter(qc_df0, !str_detect(Sample,"input")) |> 
  mutate(experiment = str_extract(Sample,"ENCSR.+(?=_rep)")) |>
  mutate(library = str_extract(Sample,"pan_.+_rep\\d")) |>
  dplyr::select(library,c("JS Distance","diff. enrichment","CHANCE divergence"))

supporting <- c('pan_ENCSR058DSI_rep1',
                'pan_ENCSR058DSI_rep2',
                'pan_ENCSR636FRF_rep1',
                'pan_ENCSR636FRF_rep2',
                'pan_ENCSR636FRF_rep3',
                'pan_ENCSR033IIP_rep1',
                'pan_ENCSR074LKQ_rep2',
                'pan_ENCSR455AWG_rep2',
                'pan_ENCSR455AWG_rep3')

toplot <- inner_join(dplyr::select(repet,library=sample,repetitiveness_index=estimate),
                     qc_df, by="library") |>
  mutate(supporting = library %in% supporting) |>
  pivot_longer(-c(library,supporting),names_to = "metric", values_to = "score")

# all normal
repet_normality <- toplot |>
  group_by(supporting,metric) |>
  summarise(data=list(score)) |>
  mutate(normality = map(data, ~broom::tidy(shapiro.test(.x)))) |>
  unnest(normality,names_sep = "_")
  
toplot |> 
ggplot(aes(supporting, score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.3) +
  facet_wrap(~metric, scales="free") +
  ggpubr::stat_compare_means(method='t.test')

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

pdf("results/figures/csem-tracks-h3k9me3-profile-repetitiveness-supplement.pdf",width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(gs$h3k9me3@ggplot, x = 0.5, y=0.5, width = 7.5,height = 3)

plotGG(g_e_repetitiveness + theme(axis.text.y = element_text(size=5)), x = 0.25, y=3.75, width = 7.5,height = 2)

plotText("A", x = 0.5, y=0.5)
plotText("B", x = 0.5, y=3.75)

dev.off()
