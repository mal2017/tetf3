library(ggtree)
library(ape)
library(tidyverse)
library(ggtreeExtra)
library(phylosignal)

rip_fl <- "results/ripseq/unr_ripseq.tsv.gz"
rip_fl <- snakemake@input$rip
rip <- read_tsv(rip_fl)

p4d_fl <- "results/ripseq/unr_ripseq_phylosignal.p4d.rds"
p4d_fl <- snakemake@input$p4d
p4d <- read_rds(p4d_fl)

tree_fl <- "results/ripseq/unr_ripseq_phylosignal.tree.rds"
tree_fl <- snakemake@input$tree
tree <- read_rds(tree_fl)

coex_df <- rip |>
  filter(type=="TE") |>
  dplyr::select(label=feature,coex) |>
  mutate(coex = if_else(coex,"Unr coex.","not coex."))

tree <- tree |> left_join(coex_df) 

gtr <- ggtree(tree,ladderize = T) +
  geom_tiplab(size=rel(0.45),as_ylab = T)

gtr <- gtr + 
  geom_fruit(aes(y=label,x=nt.content,fill=coex),geom=geom_col,
             axis.params = list(axis="x",title="AT content", title.size=2.1, line.size=0.5,text.size=1.5,vjust=0.6,title.height=0.01),
             grid.params = list(size=0.25,vline=F),
             width=1.1,
             pwidth = 0.4)

gtr <- gtr + 
  geom_fruit(aes(y=label,x=ARE_sites,fill=coex),geom=geom_col,
             axis.params = list(axis="x",title="ARE sites", title.size=2.1, line.size=0.5,text.size=1.5,vjust=0.6,title.height=0.01),
             grid.params = list(size=0.25,vline=F),
             width=1.1,
             pwidth = 0.4)
gtr <- gtr + 
  geom_fruit(aes(y=label,x=rip_log2FoldChange,fill=coex),geom=geom_col,offset = 0.5,
             axis.params = list(axis="x",title="RIP-seq log2(IP/Control)", title.size=2.1, line.size=0.5,text.size=1.5,vjust=0.6,title.height=0.01),
             grid.params = list(size=0.25,vline=F),
             width=1.1,
             pwidth = 0.4)

gtr <- gtr + 
  scale_fill_manual(values=c("not coex."="gray","Unr coex."="red"))

# ------------------------------------------------------------------------------
# phylosignal correlation
# ------------------------------------------------------------------------------

pc <- phyloCorrelogram(p4d,trait="rip_log2FoldChange")
ps_test_fl <- "results/ripseq/unr_ripseq_phylosignal.tbl.rds"
ps_test_fl <- snakemake@input$ps_test
phylosignal_tests <- read_rds(ps_test_fl)

lab <- filter(phylosignal_tests, coef == "rip_log2FoldChange" & metric == "I") |>
  pull(pval) |>
  sprintf("Moran's test p<%s",x=_)

source("workflow/scripts/utils/phylocorrelogram.R")

g_crlg <- plot_crlg(pc) + annotate("text",label=lab,y = -0.1,x=0,hjust=0,vjust=0)

write_rds(g_crlg,snakemake@output$crlg)
write_rds(gtr,snakemake@output$tree)