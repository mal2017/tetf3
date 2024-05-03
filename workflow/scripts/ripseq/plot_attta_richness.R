Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)


relpos <- 50
relpos <- snakemake@params$relpos

rip_fl <- "results/ripseq/unr_ripseq.tsv.gz"
rip_fl <-  snakemake@input$rip
rip <- read_tsv(rip_fl)

attta_fl <- "results/ripseq/unr_ripseq_features_attta_sites.tsv.gz"
attta_fl <- snakemake@input$attta

attta <- read_tsv(attta_fl) |>
  inner_join(dplyr::select(rip,feature,group,coex,type,status),by="feature")

attta_grp_n <- attta |>
  group_by(group) |>
  summarise(n=n())

attta <- attta |>
  left_join(attta_grp_n) |>
  mutate(group=sprintf("%s (n=%s)",group,n)) |>
  mutate(sites.per.100bp = sites/(seqlen/100))

plot_attta_sites <- \(x) {mutate(x, group=fct_reorder(group,sites)) |>
    ggplot(aes(group,sites)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle=25, hjust=1)) +
    ylab(sprintf("ATTTA sites\n(last %s%% of tx)",relpos)) +
    ggpubr::stat_compare_means(size=rel(2)) +
    xlab("")
}


g_attta_sites <- split(attta,attta$type) |>
  map(plot_attta_sites)

write_rds(g_attta_sites,snakemake@output$boxplot)
