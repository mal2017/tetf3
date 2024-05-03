Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)


mods_path <- "results/pirna/te_silencers_in_lms.rds"
mods_path <- snakemake@input[["mods"]]
mods <- read_rds(mods_path)

plot_te_silencer_nhits <- \(x) group_by(x,gene_symbol,model,class) |>
  summarise(n_hits = sum(significant_x),.groups = "drop") |>
  ggplot(aes(class,n_hits)) +
  geom_boxplot(outlier.shape = NA) +
  ggpubr::stat_compare_means(label.y = 9) +
  xlab("known TE regulator") +
  ylab("coexpressed TEs") +
  coord_cartesian(ylim=c(0,10))

g_f <- mods |>
  filter(model == "female") |>
  plot_te_silencer_nhits()  

g_m <- mods |>
  filter(model == "male") |>
  plot_te_silencer_nhits()

write_rds(g_f,snakemake@output$female)
write_rds(g_m,snakemake@output$male)