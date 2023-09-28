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

cutoff <- 2

#i_hits <- ps_df |>   filter(score_type == "score") |> filter(metric == "I" & pval < 0.05) |>
#  count(TF) |>
#  filter(n==2)

crlgs <- read_rds("results/phylosignal/sig_correlograms.rds")
cntrl_crlgs <- read_rds("results/phylosignal/correlograms.rds")
crlgs <- c(crlgs, cntrl_crlgs)

h0 <- -1/(crlgs[[1]]$n - 1)

crlg_df <- crlgs |>
  map(`$`,"res") |>
  map_df(as_tibble,.id="coef")

colnames(crlg_df) <- c("coef","x","ci.upper","ci.lower","y")

get_col <- \(x) c("red", "orange", "darkblue")[match(x,c("bm","ou","random"))]

g_b <- crlg_df |>
  mutate(TF = str_extract(coef,"(?<=score_).+"),
        sex = str_extract(coef,"female|male")) |>
  filter(sex == "female") |>
  mutate(TF = if_else(is.na(TF),coef, TF)) |>
  ggplot(aes(x)) +
  geom_path(aes(y=y)) +
  geom_path(aes(y=ci.upper), linetype="dotted") +
  geom_path(aes(y=ci.lower), linetype="dotted") +
  xlab("patristic distance") + ylab("coex. score correlation") +
  geom_hline(yintercept = h0, color="gray") +
  facet_wrap(~TF, scales="free") +
  theme_classic() +
  theme(text=element_text(size=5))


crlg_to_plot_df <- crlg_df |>
  mutate(TF = str_extract(coef,"(?<=score_).+"),
         sex = str_extract(coef,"female|male")) |>
  filter(sex == "female" | coef %in% c("bm","random"))  |> 
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
                           T ~ "darkblue")) |>
  mutate(linetype = case_when(coef %in% c("bm","random") ~ "dashed", T~ "solid"))
  

key <- dplyr::select(crlg_to_plot_df, label, TF) |> filter(!TF %in% c("bm","random")) |> distinct()
dfnpc <- tibble(x = 1, y = 1, tb = list(key))
g_b <- ggplot(crlg_to_plot_df,aes(x,color=color,label=label, y=y)) +
  geom_path(aes(group=TF, linetype=linetype)) +
  ggrepel::geom_text_repel(data= \(dat) {slice_min(group_by(dat,TF),x)},seed = 1,force_pull = 1000,force = 2.5,direction = "both",size=rel(3), 
                           position = position_nudge_keep(x = -0.03), max.iter = 11,color="black") +
  xlab("patristic distance") + ylab("coex. score correlation") +
  geom_hline(yintercept = h0, color="gray") +
  theme_classic() +
  theme(text=element_text(size=5)) +
  scale_color_identity() +
  geom_table_npc(mapping=aes(npcx=x, npcy=y,label=tb), data=dfnpc, table.theme = ttheme_gtminimal(base_size = 5), table.colnames = F,size=0.1) +
  scale_linetype_identity()

# ------------------------------------------------------------------------------
# get tree
# ------------------------------------------------------------------------------

tr <- read_rds("results/te_sequence_similarity/te_sketch_tidytree.rds")


g_a <- ggtree(tr, layout = "circular",right = T, ladderize = T) +
  geom_tiplab2(size=unit(1.1,"pt"),hjust = -0.5) + 
  geom_tippoint(aes(color=repClass)) +
  theme(plot.margin = unit(c(-140, -140, -170, -140), "pt"),legend.title = element_blank(),
        legend.direction = "vertical",
        legend.position = c(0.65,0.4),legend.justification = c("left","top"), plot.background = element_blank())

# ------------------------------------------------------------------------------
# make page
# ------------------------------------------------------------------------------

dir.create("results/figures/")

pdf("results/figures/figure3.pdf",width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_a, x = -0.4, y=0.75, width = 6,height = 1.5)
plotText("A", x = 0.5, y=0.5)

plotGG(g_b, x = 4.5, y=0.75, width = 3.5,height = 2.5)
plotText("B", x = 4.5, y=0.5)

dev.off()
