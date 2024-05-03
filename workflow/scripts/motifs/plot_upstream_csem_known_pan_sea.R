Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)

fl <- "results/motifs/csem_peak_sea.known.pan.tsv.gz"
fl <- snakemake@input$tsv

g <- read_tsv(fl) |>
  filter(str_detect(ID,"MA0237")) |>
  ggplot(aes(peak_set,-log10(PVALUE),fill=PVALUE<0.05)) +
  geom_col() +
  geom_hline(yintercept = -log10(0.05),linetype="dashed",color="red") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  facet_wrap(~ID,scales = "free") +
  scale_fill_manual(values=c(`TRUE`="black",`FALSE`="darkgray")) +
  xlab("modENCODE library")

write_rds(g,snakemake@output$gg)