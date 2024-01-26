library(tidyverse)
library(plotgardener)

df <- read_tsv("results/motifs/csem_peak_sea.pan.tsv.gz")

gg <- df |>
  unite(motif, ID,CONSENSUS,sep = "\n") |>
  dplyr::select(library=peak_set, motif,QVALUE) |>
  #mutate(hit = QVALUE < 0.1) |>
  #filter(QVALUE < 0.1) |>
  ggplot(aes(library,motif)) +
  geom_tile(aes(fill=QVALUE)) +
  #geom_point(color="red",size=rel(10)) +
  geom_text(data=\(x) filter(x,QVALUE < 0.1), label="*", size=10) +
  scale_fill_distiller(palette = 2) +
  theme(axis.text.x = element_text(angle=45, hjust=1))



# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

theme_set(theme_minimal() + 
            theme(text = element_text(size=7))
)

dir.create("results/figures2/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(gg, x = 1, y=0.5, width = 6.5,height = 3)

dev.off()

