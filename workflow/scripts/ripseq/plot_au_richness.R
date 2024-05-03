Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)

relpos <- 50
relpos <- snakemake@params$relpos

rip_fl <- "results/ripseq/unr_ripseq.tsv.gz"
rip_fl <-  snakemake@input$rip
rip <- read_tsv(rip_fl)

au_fl <- "results/ripseq/unr_ripseq_features_au_content.tsv.gz"
au_fl <- snakemake@input$au
au <- read_tsv(au_fl) |>
  inner_join(dplyr::select(rip,feature,group,coex,type,status),by="feature")

# au content of specific region we are quantifying and performing stats on - controlled in config.yaml
# in general is last 50% of each tx
au_in_region_fl <- "results/ripseq/unr_ripseq_features_au_content_in_region.tsv.gz"
au_in_region_fl <- snakemake@input$au_in_region
au_in_region <- read_tsv(au_in_region_fl) |>
  inner_join(dplyr::select(rip,feature,group,coex,type,status),by="feature")

au_grp_n <- au_in_region |>
  group_by(group) |>
  summarise(n=n())

# summarize per position
au_summarized <- au |> 
  group_by(position.scaled,group,type) |>
  summarise(nt.content=mean(nt.content,na.rm=T),.groups = "drop") |>
  left_join(au_grp_n, by="group") |>
  mutate(group=sprintf("%s (n=%s)",group,n))

plot_tx_at_content <- \(x) ggplot(x, aes(position.scaled,nt.content,color=group)) +
  geom_line() +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1), 
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.box="vertical") +
  guides(color = guide_legend(nrow = 2)) +
  xlab("relative transcript position") + 
  ylab(sprintf("AT content\n(last %s%% of tx)",relpos)) +
  scale_color_brewer(type="qual", palette = 2) +
  theme(legend.background=element_blank(),legend.title = element_blank()) +
  annotate("rect",xmin=relpos,xmax=100,ymin=min(au_summarized$nt.content),ymax=max(au_summarized$nt.content),fill="red",alpha=0.1)

g_au_content_over_tx_list <- split(au_summarized,au_summarized$type) |>
  map(plot_tx_at_content)


plot_au_content_in_region <- \(x) {
  left_join(x,au_grp_n, by="group") |>
    mutate(group=sprintf("%s (n=%s)",group,n)) |>
    mutate(group=fct_reorder(group,nt.content)) |>
    ggplot(aes(group,nt.content)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle=25,hjust=1)) +
    ggpubr::stat_compare_means(size=rel(2)) +
    xlab("") +
    ylab(sprintf("AT content (last %s%% of tx)",relpos))
}

g_au_content_in_region_list <- split(au_in_region,au_in_region$type) |>
  map(plot_au_content_in_region)

write_rds(g_au_content_over_tx_list,snakemake@output$lineplot)
write_rds(g_au_content_in_region_list,snakemake@output$boxplot)