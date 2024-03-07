library(tidyverse)

motif_cmp_fl <- "results/motifs/comparison/pan_denovo_comparison.meme.rds"
motif_cmp_fl <- snakemake@input$comparison
motif_comp <- read_rds(motif_cmp_fl)

motif_comp <- motif_comp |>
  group_by(denovo) |>
  arrange(Pval) |>
  mutate(rank = row_number()) |>
  ungroup() |>
  mutate(class = case_when(str_detect(known,"degenerate")~"Archbold 2014 degenerate",
                           str_detect(known,"::MA0") & name=="pan" ~"jaspar (pan)",
                           str_detect(known,"known::MA")~"jaspar (other)",
                           str_detect(known,"known::")~"Archbold 2014 HMG/helper",
                           T~"wut")) |>
  mutate(label = if_else(is.na(name),class,name))

# highest p at bh adjusted p of 0.1 - used for cutoff line in plots
max_p <- filter(motif_comp, padj < 0.1) |> pull(Pval) |> max()

motif_comp2plot <- motif_comp |>
  dplyr::select(rank, denovo, class, label, known, name, Pval, padj)
#filter(class == "pan (jaspar)")

# ------------------------------------------------------------------------------
# get known motif matches to highlight
# if jaspar is significant but not the top hit, favor the jaspar motif
# ------------------------------------------------------------------------------
jaspar_which.motifs <- motif_comp2plot |> filter(padj < 0.1 & class == "jaspar (pan)")

which.motifs <- motif_comp2plot |>
  #filter(str_length(denovo) > str_length("denovo::")+6) |>
  group_by(denovo) |>
  slice_min(padj, n=1) |>
  ungroup() |>
  #slice_max(str_length(denovo)) |>
  filter(padj < 0.1)

which.motifs <- rbind(jaspar_which.motifs,which.motifs) |>
  group_by(denovo) |>
  slice_max(class=="jaspar (pan)",n=1) |>
  dplyr::select(denovo,known) |>
  mutate(highlight=T)


# ------------------------------------------------------------------------------
# plot motif al
# ------------------------------------------------------------------------------
g_a <- motif_comp2plot |>
  filter(denovo %in% which.motifs$denovo) |>
  left_join(which.motifs) |>
  nest(data=-c(denovo)) |>
  mutate(g_rnk=map(data,~{
    g <- ggplot(.x,aes(rank,-log10(Pval))) +
      #geom_point(data=\(x){filter(x, padj >= 0.1)}, size=rel(0.75),color="gray") +
      geom_point(size=rel(1)) +
      ggrepel::geom_text_repel(data= \(x) filter(x,highlight), aes(label=label),color="black") +
      scale_color_brewer(type = "qual", palette = 2,name="") +
      geom_hline(yintercept = -log10(max_p), linetype="dashed", color="darkgray") +
      #facet_wrap(~denovo, scales = "free", nrow=1) +
      theme(legend.position = "bottom")
    return(g)
  }))
  
g_motif_alns <- filter(motif_comp, denovo %in% g_a$denovo) |>
  right_join(which.motifs) |>
  group_by(denovo) |>
  dplyr::select(denovo,class,g_aln=gg) #|>
  #Reduce(`+`,x=_ ) & theme_bw() & 
  #guides(color="none", fill="none") & 
  #theme(text = element_text(size=5), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank())

res <- left_join(g_a,g_motif_alns,by="denovo")

write_rds(res, snakemake@output$gg)
