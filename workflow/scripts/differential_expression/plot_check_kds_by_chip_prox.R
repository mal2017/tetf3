library(tidyverse)


d_fl <- "results/pirna/encode_peaks_dist_to_pirna.gr.rds"
d_fl <- snakemake@input$gr
d <- read_rds(d_fl)


res_fl <- "results/deg/ourKD.de.df.rds"
res_fl <- snakemake@input$res
res_us <- read_rds(res_fl)

res_us <- res_us |>
  #filter(group %in% c("head")) |>
  mutate(target = str_extract(comparison,"pan|NfI|vvl|CG16779"))

df <- d |>
  as_tibble() |>
  right_join(res_us, by=c(gene_id = "feature",ChIP="target")) |>
  filter(!is.na(distance)) # tes, as well as genes on some scaffolds/chroms without peaks have na distances

g <- df |>
  filter(embryo.expressed) |>
  mutate(lab = str_replace_all(str_remove(comparison,"knockdown2_"),"_","/")) |>
  ggplot(aes(distance,-log10(pvalue))) +
  geom_point(size=1) +
  facet_wrap(~lab) +
  xlab("distance to RNAi target's nearest modENCODE ChIP peak (embryo)")



write_rds(g,snakemake@output$gg)
  