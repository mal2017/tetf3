library(patchwork)
library(plotgardener)
library(tidyverse)
library(ggpubr)
library(ggrastr)

# ------------------------------------------------------------------------------
# read in plots/data used for multiple figs or requireing some ggplot aesthetic tweaks
# ------------------------------------------------------------------------------

repcor <- read_rds("results/coexpression_replication/replicate_dataset_correlation.rds")
mf <- read_rds("results/coexpression_replication/mf_dataset_correlation.rds")
var_exp <- read_rds("results/exploratory_and_descriptive/variance_explained.gg.rds") 

# ------------------------------------------------------------------------------
# plotting
# ------------------------------------------------------------------------------

# replicate correlation
g_repcor <- filter(repcor, result_set=="main_data") |>
  dplyr::select(model,gg) |>
  deframe() |>
  imap(~{.x + labs(title=.y)})

# male vs female
g_mf <- mf |>
  filter(result_set %in% c("unfiltered","replicated")) |>
  dplyr::select(result_set,gg) |>
  deframe() |>
  imap(~{.x + labs(title=.y)})

# variance explained
g_var_exp <- var_exp + 
  #theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle=25, hjust=1)) +
  scale_fill_grey(start = 0.4, end=0.8)

g_nhits_prev_rep_teregs <- read_rds("results/pirna/te_silencer_n_hits_boxplot.males.gg.rds")

g_score_prev_rep_teregs <- read_rds("results/pirna/te_silencer_scores_boxplot.males.gg.rds")

g_sharedness <- read_rds("results/exploratory_and_descriptive/mf_shared_coex_hits.gg.rds") + 
  theme(axis.text.x=element_text(angle=45, hjust=1))

# plotting --------------------------------------------------------------------
theme_set(theme_classic() + theme(text=element_text(size=unit(7,"pt"))))

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 8.5, height = 11)

pageCreate(width =8.5, height = 11, default.units = "inches", showGuides = interactive())
figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

pa <- plotGG(plot = g_var_exp, x = 0.5, y=0.75, width = 7.5, height=2.75)
plotText(label = "A", x = 0.5, y = 0.5)

pb <- plotGG(plot = g_mf$unfiltered+ scale_fill_distiller(palette = 6) + guides(fill="none"), x = 0.4,  y=3.75, width = 1.75, height=1.75)
pc <- plotGG(plot = g_mf$replicated+ scale_fill_distiller(palette = 6) + guides(fill="none"), x = 2.3,  y=3.75, width = 1.75, height=1.75)
pd <- plotGG(plot = g_repcor$male + scale_fill_distiller(palette = 6) + guides(fill="none"), x = 4.3,  y=3.75, width = 1.75, height=1.75)
pe <- plotGG(plot = g_repcor$female+ scale_fill_distiller(palette = 6) + guides(fill="none"), x = 6.3,  y=3.75, width = 1.75, height=1.75)


plotText(label = "B", x = 0.5, y = 3.75)
plotText(label = "C", x = 2.3, y = 3.75)
plotText(label = "D", x = 4.3, y = 3.75)
plotText(label = "E", x = 6.3, y = 3.75)

plotGG(g_sharedness, x = 0.5, y=5.75, width = 2,height = 2)
plotText("F",  x = 0.5, y=5.75)

plotGG(g_score_prev_rep_teregs, x = 2.45, y=5.75, width = 3.45,height = 2)
plotText("G",  x = 2.5, y=5.75)

plotGG(g_nhits_prev_rep_teregs, x = 6, y=5.75, width = 2,height = 2)
plotText("H",  x = 6, y=5.75)

dev.off()

writexl::write_xlsx(list(`A (summarized to median for xlsx compatibility)`=group_by(g_var_exp$data,model,coef) |> summarise(med.var.explained=median(var.explained,na.rm=T)),
                         `B (subsampled for xlsx compatibility)`=sample_n(g_mf$unfiltered$data,replace = F,size = 999999),
                         C=g_mf$replicated$data,
                         D=g_repcor$male$data,
                         E=g_repcor$female$data,
                         `F`=g_sharedness$data,
                         G=g_score_prev_rep_teregs$data,
                         H=g_nhits_prev_rep_teregs$data),
                    path = ifelse(exists("snakemake"),
                                  snakemake@output$xlsx,
                                  "~/Downloads/test.xlsx"))
