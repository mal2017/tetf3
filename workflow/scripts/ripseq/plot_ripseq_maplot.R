Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)

rip_fl <- "results/ripseq/unr_ripseq.tsv.gz"
rip_fl <- snakemake@input$rip
rip <- read_tsv(rip_fl)

malab <- rip |>
  filter(status=="bound") |>
  group_by(type,status) |>
  tally() |>
  mutate(lab = paste0(n," ",status," ", type,"s")) |>
  pull(lab) |>
  paste(collapse="\n")

g_ma_plot <- rip |>
  filter(comparison=="WT_vs_CTRL") |>
  arrange(desc(type)) |>
  ggplot(aes(rip_baseMean+1,rip_log2FoldChange,color=status)) +
  geom_point(size=rel(0.5)) +
  geom_point(data=\(x)filter(x,type=="TE"),size=rel(1),color="red",fill=NA,shape=21) +
  scale_x_log10() +
  scale_color_manual(values=c(bound="antiquewhite4",not_bound="lightgray")) +
  #scale_fill_manual(values=c(bound="red",not_bound=NA)) +
  xlab("mean normalized counts") + ylab("log2(IP/control)") +
  geom_hline(yintercept = 0, linetype="dashed",color="black") +
  annotate("text",label=malab,y=max(rip$rip_log2FoldChange)*0.9,x=median(rip$rip_baseMean)*10,size=rel(2)) +
  guides(color="none")

write_rds(g_ma_plot,snakemake@output$gg)
