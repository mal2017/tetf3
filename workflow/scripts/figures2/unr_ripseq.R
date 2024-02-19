library(tidyverse)
library(patchwork)
library(plotgardener)
library(clusterProfiler)
library(enrichplot)
library(phylosignal,quietly = T)
library(phylobase,quietly = T)
library(ggtree)
library(tidytree)
library(ape)

relpos <- 50
relpos <- snakemake@params$relpos

rip <- read_tsv("results/ripseq/unr_ripseq.tsv.gz")

au <- read_tsv("results/ripseq/unr_ripseq_features_au_content.tsv.gz") |>
  inner_join(dplyr::select(rip,feature,group,coex,type,status),by="feature")

au_in_region <- read_tsv("results/ripseq/unr_ripseq_features_au_content_in_region.tsv.gz") |>
  inner_join(dplyr::select(rip,feature,group,coex,type,status),by="feature")

attta <- read_tsv("results/ripseq/unr_ripseq_features_attta_sites.tsv.gz") |>
  inner_join(dplyr::select(rip,feature,group,coex,type,status),by="feature")


# ------------------------------------------------------------------------------
# ma-plot from RIP-seq data
# ------------------------------------------------------------------------------

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


# ------------------------------------------------------------------------------
# leading edge in knockdown
# ------------------------------------------------------------------------------

gsea <- read_rds("results/ripseq/unr_bound_tx_in_kd.gsea.rds")

g_bound_te_gsea <- gsea |> enrichplot::gseaplot2("bound_TE") & 
  theme(axis.title = element_text(size=5),axis.text = element_text(size=5))
g_bound_gene_gsea <- gsea |> enrichplot::gseaplot2("bound_gene") &
  theme(axis.title = element_text(size=5),axis.text = element_text(size=5))

# ------------------------------------------------------------------------------
# au richness
# ------------------------------------------------------------------------------

au_grp_n <- au_in_region |>
  group_by(group) |>
  summarise(n=n())

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
  xlab("relative transcript position") + ylab("AT content") +
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

# ------------------------------------------------------------------------------
# attta sites
# ------------------------------------------------------------------------------

plot_attta_sites <- \(x) {mutate(x, group=fct_reorder(group,sites)) |>
  ggplot(aes(group,sites)) +
    geom_boxplot() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ylab(sprintf("ATTTA sites (last %s%% of tx)",relpos)) +
  ggpubr::stat_compare_means(size=rel(2)) +
  xlab("")
}


g_attta_sites <- split(attta,attta$type) |>
  map(plot_attta_sites)


# ------------------------------------------------------------------------------
# phylosignal tree
# ------------------------------------------------------------------------------
library(ggtreeExtra)
p4d <- read_rds("results/ripseq/unr_ripseq_phylosignal.p4d.rds")
tree <- read_rds("results/ripseq/unr_ripseq_phylosignal.tree.rds")

coex_df <- rip |>
  filter(type=="TE") |>
  dplyr::select(label=feature,coex) |>
  mutate(coex = if_else(coex,"Unr coex.","not coex."))

tree <- tree |> left_join(coex_df) 


#phylosignal::dotplot.phylo4d(p4d,tree.ladderize = T,
#                             trait = c("rip_stat","ARE_sites","nt.content","male","female"),
#                             scale = T, center=T,tip.cex = 0.5)

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
phylosignal_tests <- read_rds("results/ripseq/unr_ripseq_phylosignal.tbl.rds")

lab <- filter(phylosignal_tests, coef == "rip_log2FoldChange" & metric == "I") |>
  pull(pval) |>
  sprintf("Moran's test p<%s",x=_)

source("workflow/scripts/utils/phylocorrelogram.R")

g_crlg <- plot_crlg(pc) + annotate("text",label=lab,y = -0.1,x=0,hjust=0,vjust=0)


# ------------------------------------------------------------------------------
# create page 1
# ------------------------------------------------------------------------------
theme_set(theme_classic() + 
            theme(text = element_text(size=7)) +
            theme(plot.title = element_text(hjust = 0.5))
          )

dir.create("results/figures2/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_ma_plot, x = .5, y=0.5, width = 2.5,height = 2)
plotText("A", x = .5, y=0.5)


plotGG(g_bound_gene_gsea + xlab("rank in Unr knockdown"), x = 3, y=0.5, width = 2.5,height = 2.25)
plotText("Unr-bound genes",x=4.75,y=0.75,fontsize = 5)
plotText("B", x = 3, y=0.5)

plotGG(g_bound_te_gsea + xlab("rank in Unr knockdown"), x = 5.5, y=0.5, width = 2.5,height = 2.25)
plotText("Unr-bound TEs",x=7.25,y=0.75,fontsize = 5)
plotText("C", x = 5.5, y=0.5)

plotGG(g_au_content_over_tx_list$gene, x = 0.5, y=2.75, width = 5.25,height = 2)
plotText("D",  x = 0.5, y=2.75)

plotGG(g_au_content_over_tx_list$TE + theme(legend.position = c(0,0.4)),
       x = 0.5, y=5, width = 5.25,height = 2)
plotText("F",  x = 0.5, y=5)

plotGG(g_au_content_in_region_list$gene + theme(legend.position = c(0,0.4)),
       x = 5.75, y=2.75, width = 2.25,height = 2.25)
plotText("E",  x = 5.7, y=2.75)


plotGG(g_au_content_in_region_list$TE + theme(legend.position = c(0,0.4)),
       x = 5.75, y=5.1, width = 2.25,height = 2.25)
plotText("G",  x = 5.7, y=5.1)


plotGG(g_attta_sites$gene,
       x = 0.5+0.78125, y=7.5, width = 2.25*0.75,height = 2.25)
plotText("H",  x = 0.5+0.78125, y=7.5)

plotGG(g_attta_sites$TE,
       x = 0.5+0.78125 + 2.25*0.75 + 0.25, y=7.5, width = 5*0.75,height = 2.61)
plotText("I",  x =0.5+0.78125 + 2.25*0.75 + 0.25, y=7.5)

dev.off()

# ------------------------------------------------------------------------------
# create page 2
# ------------------------------------------------------------------------------


pdf(snakemake@output$pdf2,width = 8.5, height = 11)
pageCreate(height = 11, showGuides=interactive())

plotGG(gtr+guides(fill="none"), x = .5, y=0.5, width = 7,height = 5.5)
plotText("A", x = .5, y=0.5)

plotGG(g_crlg, x = 2, y=6.5, width = 4,height = 2)
plotText("B", x = 2, y=6.5)

dev.off()