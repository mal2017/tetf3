library(tidyverse)

motif_cmp_fl <- "results/motifs/comparison/pan_denovo_comparison.meme.rds"
motif_cmp_fl <- snakemake@input$comparison
motif_comp <- read_rds(motif_cmp_fl)

motif_comp <- motif_comp |>
  group_by(denovo) |>
  arrange(Pval) |>
  mutate(rank = row_number()) |>
  ungroup() |>
  mutate(class = case_when(str_detect(known,"Degenerate")~"Archbold 2014 degenerate",
                           str_detect(known,"::MA0") & name=="pan" ~"jaspar (pan)",
                           str_detect(known,"known::MA")~"jaspar (other)",
                           str_detect(known,"known::")~"Archbold 2014 HMG/helper",
                           T~"wut")) |>
  mutate(label = if_else(is.na(name),class,name))

# ------------------------------------------------------------------------------
# get known motif matches to highlight
# if jaspar is significant but not the top hit, favor the jaspar motif
# ------------------------------------------------------------------------------

# highest p at bh adjusted p of 0.1 - used for cutoff line in plots


# used for plotting the rank of matches
motif_comp2plot <- motif_comp |>
  dplyr::select(rank, denovo, class, label, known, name, Pval, padj)

which.motifs <- motif_comp2plot |> 
  filter(padj < 0.1) |>
  filter(class %in% c("jaspar (pan)")|str_detect(class,"Archbold")) |>
  dplyr::select(denovo,known)

g_rnks <- which.motifs |>
  mutate(g_rnk = map2(denovo,known,
                     ~{
                       dat <- filter(motif_comp2plot,denovo==.x)
                       max_p <- filter(dat, padj < 0.1) |> pull(Pval) |> max()
                       ggplot(dat,aes(rank,-log10(Pval))) +
                         geom_point() +
                         ggrepel::geom_text_repel(data= \(x) filter(x,known==.y), aes(label=label),color="black") +
                         geom_hline(yintercept = -log10(max_p), linetype="dashed", color="darkgray")
                       
                       }))


# --------------------------------
# select only motifs we're highlighting and simplify the results table
# this now has info about the significance of the comparison,
# alignment of the compared motifs, and a ranked dotplot with best known pan motif
# highlighted
# --------------------------------
res <- motif_comp |> 
  inner_join(which.motifs) |>
  dplyr::select(denovo,known,name,class,label,Pval,Eval,padj,motifs,g_aln=gg) |>
  left_join(g_rnks,by=c("denovo","known"))

write_rds(res, snakemake@output$gg)
