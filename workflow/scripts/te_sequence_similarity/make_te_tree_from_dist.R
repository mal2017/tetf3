Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(tidytree)
library(ggtree) # required for conversion to treedata, oddly enough

# get te classes
te.classes.fl <- "~/work/tetf3/resources/Tidalbase_Dmel_TE_classifications_2015.txt"
te.classes.fl <- snakemake@input$te_classes

te.classes <- read_tsv(te.classes.fl) %>%
  dplyr::select(feature.y = Flybase_name, Class, repClass, repFamily) %>%
  distinct()

# fix some missing info
te.classes <- mutate(te.classes, repFamily = ifelse(feature.y == "TART-C","Jockey",repFamily)) %>%
  mutate(repFamily = ifelse(feature.y == "Tc3","TcMar-Tc1",repFamily)) %>%
  mutate(repFamily = ifelse(feature.y %in% c("TLD2_LTR","Stalker3T"),"Gypsy",repFamily)) %>%
  mutate(repClass = ifelse(feature.y == "TART-C","LINE",repClass)) %>%
  mutate(repClass = ifelse(feature.y == "Tc3","DNA",repClass)) |>
  mutate(repClass = ifelse(feature.y %in% c("TLD2_LTR","Stalker3T"),"LTR",repClass)) 

# get dist obj
dist_fl <- "results/te_sequence_similarity/te_female-coex_dist.rds"
dist_fl <- snakemake@input$dist

d <- read_rds(dist_fl)

phylo <- ape::bionjs(d)

tr <- as.treedata(phylo) |>
  tidytree::left_join(te.classes, by=c(label="feature.y"))

#ggtree(tr, layout = "radial") +
#  ggtree::geom_tiplab(aes(label=label), size=rel(3)) +
#  ggtree::geom_tippoint(aes(color=repClass), size=rel(3)) +
#  paletteer::scale_color_paletteer_d("ggsci::default_igv")


write_rds(tr,snakemake@output$tidytree)

# for future reference - hears an easy conversion of
# an unrooted non-ultrametric tree to an hclust obj
#phylo |>
#  phytools::reroot(which(phylo$tip.label == "INE-1")) |>
# phytools::force.ultrametric(method="extend") |>
# ape::as.hclust.phylo() |>
#  plot()
