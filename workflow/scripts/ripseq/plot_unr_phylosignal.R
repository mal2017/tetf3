library(ggtree)
library(ape)
library(tidyverse)
library(ggtreeExtra)
library(phylosignal)
library(ggpp)

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


labs <- phylosignal_tests |>
  filter(coef == "rip_log2FoldChange") |>
  mutate(label = sprintf("metric: %s\np=%s\nadj. p=%s",metric,round(pval,3),round(padj,3))) |>
  filter(metric == "Cmean")
# ------------------------------------------------------------------------------
# get expected value -  see phylosignal package for calculation
# ------------------------------------------------------------------------------
h0 <- -1/(pc$n - 1) # see the phylosignal paper/code for the expected null value calculation

# ------------------------------------------------------------------------------
# get phylosignal's correlogram data as df for ggplot plotting
# ------------------------------------------------------------------------------
crlg_df <- list(rip_log2FoldChange = pc) |>
  map(`$`,"res") |>
  map_df(as_tibble,.id="coef") |>
  set_names(c("coef","x","ci.upper","ci.lower","y"))

crlg_to_plot_df <- crlg_df |>
  ungroup() |>
  mutate(color="darkgreen") |>
  mutate(linetype = "solid") |>
  mutate(type = if_else(linetype == "dashed","control",coef)) |>
  left_join(labs)


# ------------------------------------------------------------------------------
# plotting func
# ------------------------------------------------------------------------------
plot_crlg <- function(df) {
  key <- dplyr::select(df, coef, label) |> distinct()
  dfnpc <- tibble(x = 1, y = 1, tb = list(key))
  g <- ggplot(df,aes(x,color=color,label=coef, y=y)) +
    geom_path(aes(group=label, linetype=linetype)) +
    geom_path(aes(y=ci.upper), linetype="dotted") +
    geom_path(aes(y=ci.lower), linetype="dotted") +
    ggrepel::geom_text_repel(data= \(dat) {slice_min(group_by(dat,coef),x)},seed = 1,force_pull = 1000,force = 2.5,direction = "both",size=rel(3), 
                             position = position_nudge_keep(x = -0.03), max.iter = 11,color="black") +
    xlab("patristic distance") + ylab("coex. score autocorrelation") +
    geom_hline(yintercept = h0, color="gray") +
    theme_classic() +
    theme(text=element_text(size=5)) +
    scale_color_identity() +
    #geom_table_npc(mapping=aes(npcx=x, npcy=y,label=tb), data=dfnpc, table.theme = ttheme_gtminimal(base_size = 5), table.colnames = F,size=0.1) +
    scale_linetype_identity() +
    ylim(c(-0.1,0.125))
  
  if (!unique(df$label) %in% c("bm","random")) {
    g <- g + annotate("text",x=Inf,y=Inf,label=unique(df$label),hjust=1,vjust=1,size=2)
  }
  
  return(g)
}


g_crlg <- plot_crlg(crlg_to_plot_df)

write_rds(g_crlg,snakemake@output$crlg)
write_rds(gtr,snakemake@output$tree)