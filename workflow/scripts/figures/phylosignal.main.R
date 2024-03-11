library(tidyverse)
library(patchwork)
library(plotgardener)
library(tidyverse)
library(clusterProfiler)
library(ggpubr)
library(phylosignal)
library(phylobase)
library(ggplotify)
library(ggtree)
library(tidytree)
library(ggtreeExtra)
library(ggnewscale)
library(paletteer)
library(ggpp)


# ------------------------------------------------------------------------------
# tf correlograms
# ------------------------------------------------------------------------------
ps_df <- read_rds("results/phylosignal/phylosignal_df.rds")
crlgs <- read_rds("results/phylosignal/goi_correlograms.rds")
cntrl_crlgs <- read_rds("results/phylosignal/control_correlograms.rds")
crlgs <- c(crlgs, cntrl_crlgs)

filtered_hits <- read_tsv("results/phylosignal/phylosignal_filtered_hits.tsv.gz")

h0 <- -1/(crlgs[[1]]$n - 1) # see the phylosignal paper/code for the expected null value calculation

crlg_df <- crlgs |>
  map(`$`,"res") |>
  map_df(as_tibble,.id="coef") |>
  set_names(c("coef","x","ci.upper","ci.lower","y"))

crlg_to_plot_df <- crlg_df |>
  mutate(TF = str_extract(coef,"(?<=score_).+"),
         sex = str_extract(coef,"female|male")) |>
  filter(coef!="ou") |>
  #filter(sex == "female" | coef %in% c("bm","random"))  |> 
  mutate(TF = if_else(is.na(TF),coef, TF)) |>
  group_by(TF) |>
  mutate(max_MoranI = head(y,1)) |>
  ungroup() |>
  mutate(TF = fct_reorder(TF, -max_MoranI)) |>
  mutate(TF = fct_relevel(TF, c("bm","random"),after=Inf)) |>
  arrange(TF) |>
  group_by(TF) |>
  mutate(id = LETTERS[cur_group_id()]) |>
  mutate(label=case_when(TF %in% c("bm","random")~TF,
                         T~id)) |>
  mutate(color = case_when(coef == "bm" ~ "darkgreen",
                           coef == "random"~ "red",
                           coef %in% filtered_hits$coef~"darkgreen",
                           T ~ "darkgray")) |>
  mutate(linetype = case_when(coef %in% c("bm","random") ~ "dashed", T~ "solid")) |>
  mutate(type = if_else(linetype == "dashed","control",coef))

plot_crlg <- function(df) {
  key <- dplyr::select(df, label, TF) |> filter(!TF %in% c("bm","random")) |> distinct()
  dfnpc <- tibble(x = 1, y = 1, tb = list(key))
  ggplot(df,aes(x,color=color,label=TF, y=y)) +
    geom_path(aes(group=TF, linetype=linetype)) +
    geom_path(aes(y=ci.upper), linetype="dotted") +
    geom_path(aes(y=ci.lower), linetype="dotted") +
    ggrepel::geom_text_repel(data= \(dat) {slice_min(group_by(dat,TF),x)},seed = 1,force_pull = 1000,force = 2.5,direction = "both",size=rel(3), 
                             position = position_nudge_keep(x = -0.03), max.iter = 11,color="black") +
    xlab("patristic distance") + ylab("coex. score autocorrelation") +
    geom_hline(yintercept = h0, color="gray") +
    theme_classic() +
    theme(text=element_text(size=5)) +
    scale_color_identity() +
    #geom_table_npc(mapping=aes(npcx=x, npcy=y,label=tb), data=dfnpc, table.theme = ttheme_gtminimal(base_size = 5), table.colnames = F,size=0.1) +
    scale_linetype_identity() +
    ylim(c(-0.1,0.125))
}


crlg_gs <- split(crlg_to_plot_df, crlg_to_plot_df$type) |>
  map(plot_crlg)


g_b <- crlg_gs$control


ord <- c("female_score_%s","male_score_%s") |>
  map(~sprintf(.x,c("pan","vvl","NfI","CG16779","Unr"))) |>
  unlist()

g_c <- Reduce(`+`, crlg_gs[! names(crlg_gs) %in% c("control")][ord]) + plot_layout(ncol=5, nrow=2)


# ------------------------------------------------------------------------------
# get tree
# ------------------------------------------------------------------------------

tr <- read_rds("results/te_sequence_similarity/te_sketch_tidytree.rds")

tr <- tr |> mutate(repClass= if_else(label %in% c("Stalker3T","TLD2_LTR"),"LTR",repClass))

g_a <- ggtree(tr, layout = "rectangular",right = T, ladderize = T, size=0.25) +
  #geom_tiplab2(size=unit(1.1,"pt"),hjust = -0.5) + 
  geom_tippoint(aes(color=repClass)) +
  theme(#plot.margin = unit(c(-140, -140, -170, -140), "pt"),
        line = element_line(linewidth=0.001),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.position = c(0.65,0.45),legend.justification = c("left","top"), plot.background = element_blank())

# ------------------------------------------------------------------------------
# make page
# ------------------------------------------------------------------------------

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_a, x = 0.5, y=0.5, width = 4,height = 2.65)
plotText("A", x = 0.5, y=0.5)

plotGG(g_b, x = 4.5, y=0.65, width = 3.5,height = 2.5)
plotText("B", x = 4.5, y=0.5)

plotGG(g_c, x = .25, y=3.25, width = 8,height = 3.5)
plotText("C", x = .5, y=3.25)
plotText("D", x = .5, y=5)


dev.off()
