library(tidyverse)
library(tidytree)
library(ggtree)
library(treeio)
library(plyranges)
library(ggtreeExtra)
library(ggnewscale)
library(memes)

# tree to be plotted
tre_fl <- "~/work/tetf3/results/te_sequence_similarity/te_sketch_tidytree.rds"
tre_fl <- snakemake@input$tidytree
tre <- read_rds(tre_fl)

# get info about pan coexpression
mods_fl <- "~/work/tetf3/upstream/final-models.collected-info.tsv.gz"
mods_fl <- snakemake@input$coex
mods <- read_tsv(mods_fl)

# get info on the tf we're working with - only going to be pan unless we need to 
# examine others later
tf_symbol <- "pan"
tf_symbol <- snakemake@params$tf
coex <- filter(mods, gene_symbol == tf_symbol)

coex <- coex |>
  dplyr::select(label=feature.y, sex=model, significant_x, coex.score = estimate.qnorm)

# get info about motif presence from fimo
fimo_fl <- "results/motifs/fimo_on_tes/denovo/pan/fimo.tsv"
fimo_fl <- paste0(snakemake@input$fimo,"/fimo.tsv")
motif_instances <- read_tsv( fimo_fl,comment="#")

motif_hits <- motif_instances |>
  #dplyr::select(label=sequence_name, denovo.motif=motif_alt_id) |>
  #mutate(has.motif=T) |>
  group_by(motif_id, motif_alt_id, sequence_name) |>
  summarise(n=as.numeric(n()), .groups = "drop") |> # cfor continuous scale later
  mutate(denovo.motif=motif_alt_id) |>
  mutate(id = str_extract(motif_alt_id,"\\d+"))

# make tree
g <- tre |> 
  ggtree(layout = "circular", branch.length = "none")

# add tiles with coexpression score for each sex
g <-  g +
  geom_fruit(data = filter(coex, sex=="male"), geom="geom_col", aes(y=label,fill=abs(coex.score), x=abs(coex.score)),
             axis.params=list(
               axis       = "x",
               text.size  = 1.8,
               hjust      = 1,
               vjust      = 0.5,
               nbreak     = 3,
             )) +
  scale_fill_gradient2(low="darkgray",high="red", name="abs. male coex. score")  +
  geom_fruit(data = filter(coex, sex=="male" & significant_x), geom="geom_point", aes(y=label), color="red")
  
g <- g + new_scale_fill() +
  geom_fruit(data = filter(coex, sex=="female"), geom="geom_col", aes(y=label,fill=abs(coex.score), x=abs(coex.score)),
             axis.params=list(
               axis       = "x",
               text.size  = 1.8,
               hjust      = 1,
               vjust      = 0.5,
               nbreak     = 3,
             )) +
  scale_fill_gradient2(low="darkgray",high="blue",name="abs. female coex. score") +
  geom_fruit(data = filter(coex, sex=="female" & significant_x), geom="geom_point", aes(y=label), color="blue")

g <- g + new_scale_fill() + 
  geom_fruit(data = filter(motif_hits, motif_alt_id == "MEME-6"), aes(y=sequence_name,x=n, fill=n), geom="geom_col",
             axis.params=list(
               axis       = "x",
               text.size  = 1.8,
               hjust      = 1,
               vjust      = 0.5,
               nbreak     = 3,
             )) +
  scale_fill_distiller(name="n MEME-6 motif hits", palette = 6, direction=1)
  scale_fill_manual()
  

#g + geom_tiplab(offset = 17) 

saveRDS(g, snakemake@output$rds)

