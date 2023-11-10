library(tidyverse)
library(universalmotif)
library(patchwork)
library(plotgardener)

sils <- read_rds("results/motifs/known_motifs/archbold_degenerate_clustering.rds")

motifs <- read_meme("results/motifs/known_motifs/archbold_degenerate.meme")

names(motifs) <- map_chr(motifs,~.x@name)

g_memes <-  motifs |>
  imap(~{view_motifs(.x, type="PPM", relative_entropy=F, normalize.score=T) + labs(title = .y)}) |>
  Reduce(`+`,x=_) & theme(text=element_text(size=7))

g_sil_comparison <- sils |>
  ggplot(aes(as.factor(klet_size),as.factor(centers), size=mean_width, color=n_mislabeled)) +
  geom_point() +
  scale_color_viridis_b(option = 'C', name='n mislabeled') +
  scale_size(name="mean sil. width") +
  xlab("klet size") + ylab("n clusters")


best_sil <- filter(sils, n_mislabeled == 0) |>
  slice_max(mean_width, n=1, with_ties = T) |>
  slice_max(centers, n=1, with_ties = F)

g_best_sil <- best_sil |> pull(gg) |> pluck(1)


# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=7))
)

dir.create("results/figures/")

pdf("results/figures/archbold_merging_supplement-01.pdf",width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_sil_comparison, x = 0.5, y=0.5, width = 4,height = 2.5)
plotText("A", x = 0.5, y=0.5)

plotGG(g_best_sil, x = 4.75, y=0.5, width = 3.25,height = 2.5)
plotText("B", x = 4.75, y=.5)

plotGG(g_memes, x=0.25, y=3.25, width = 7.75, height = 1.25)
plotText("C", x = 0.5, y=3.1)


dev.off()