library(tidyverse)
library(tidytree)
library(ggtree)
library(treeio)
library(plyranges)
library(ggtreeExtra)
library(ggnewscale)

# tree to be plotted
tre_fl <- "~/work/tetf3/results/te_sequence_similarity/te_sketch_tidytree.rds"
tre_fl <- snakemake@input$tidytree
tre <- read_rds(tre_fl)

# get info about pan coexpression
mods_fl <- "~/work/tetf3/upstream/final-models.collected-info.tsv.gz"
mods_fl <- snakemake@input$coex
mods <- read_tsv(mods_fl)


tf_symbol <- "pan"
tf_symbol <- snakemake@params$tf
coex <- filter(mods, gene_symbol == tf_symbol)

coex <- coex |>
  dplyr::select(label=feature.y, sex=model, significant_x, coex.score = estimate.qnorm)

# get info about motif presence from sea
#motif_instances <- read_tsv("~/work/tetf3/results/motifs/xstreme_per_tf/pan/sea_out/sequences.tsv",comment = "#")
fimo_fl <- "results/motifs/fimo_on_tes/pan/fimo.tsv"
fimo_fl <- paste0(snakemake@input$fimo,"/fimo.tsv")
motif_instances <- read_tsv( fimo_fl,comment="#")


motif_hits <- motif_instances |>  
  dplyr::select(label=sequence_name, denovo.motif=motif_alt_id) |>
  mutate(has.motif=T) |>
  mutate(id = str_extract(denovo.motif,"\\d+"))

# make tree
g <- tre |> 
  ggtree(layout = "fan", branch.length = "none", open.angle = 10)

g <- g + geom_tippoint(aes(color=repFamily),size=rel(4))

# add tiles with coexpression score for each sex
#g <-  g +
#  geom_fruit(data = filter(coex, sex=="female"), geom="geom_tile", aes(y=label,fill=coex.score),pwidth = 2, offset=0.1) +
#  scale_fill_gradient2(low="darkgray",high="blue",name="female coex. score") +
#  new_scale_fill() +
#  geom_fruit(data = filter(coex, sex=="male"), geom="geom_tile", aes(y=label,fill=coex.score),pwidth = 2,offset = 0.1) +
# scale_fill_gradient2(low="darkgray",high="red", name="male coex. score")

# add point indicating signficant coexpression signal
g <- g +
  geom_fruit(data = filter(coex, sex=="female" & significant_x), geom="geom_point", aes(y=label), color="blue", offset = 0.1, size=rel(3.5)) +
  geom_fruit(data = filter(coex, sex=="male" & significant_x), geom="geom_point", aes(y=label), color="red", offset = 0.1, size=rel(3.5))

# add more tiles indicating motif presence/absence for each TE
fruitlist <- geom_fruit_list(
  geom_fruit(data=filter(motif_hits, has.motif), geom="geom_tile", aes(y=label, x=denovo.motif),color="darkgray", fill="darkgray",size=rel(0.5)),
  geom_fruit(data=filter(motif_hits, has.motif), geom="geom_text", aes(y=label, x=denovo.motif,label=id),color="black",fill=NA, size=rel(2))
)


g <- g + fruitlist #+ geom_tiplab(aes(label=label),offset = 15.2)

saveRDS(g, snakemake@output$rds)

