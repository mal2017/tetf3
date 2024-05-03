Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)
library(phylosignal)
library(phylobase)
library(ggpp)

# ------------------------------------------------------------------------------
# tf correlograms
# ------------------------------------------------------------------------------
crlg_fl <- "results/phylosignal/goi_correlograms.rds"
crlg_fl <- snakemake@input$crlg
crlgs <- read_rds(crlg_fl)

filtered_fl <- "results/phylosignal/phylosignal_filtered_hits.tsv.gz"
filtered_fl <- snakemake@input$filt
filtered_hits <- read_tsv(filtered_fl)

# ------------------------------------------------------------------------------
# get expected value -  see phylosignal package for calculation
# ------------------------------------------------------------------------------
h0 <- -1/(crlgs[[1]]$n - 1) # see the phylosignal paper/code for the expected null value calculation

# ------------------------------------------------------------------------------
# get phylosignal's correlogram data as df for ggplot plotting
# ------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# annotate this with the best test for each factor
# and filter to only plot the TFs we've filtered
#-------------------------------------------------------------------------------

crlg_to_plot_df <- crlg_to_plot_df |>
  left_join(filtered_hits, by=c("coef","TF","sex")) |>
  mutate(label = sprintf("sex: %s\nbest metric: %s\np=%s\nadj. p=%s",sex,metric,round(pval,3),round(padj,3)))

# ------------------------------------------------------------------------------
# plotting func
# ------------------------------------------------------------------------------
plot_crlg <- function(df) {
  key <- dplyr::select(df, label, TF) |> filter(!TF %in% c("bm","random")) |> distinct()
  dfnpc <- tibble(x = 1, y = 1, tb = list(key))
  g <- ggplot(df,aes(x,color=color,label=TF, y=y)) +
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
  
  if (!unique(df$TF) %in% c("bm","random")) {
    g <- g + annotate("text",x=Inf,y=Inf,label=unique(df$label),hjust=1,vjust=1,size=2)
  }
  
  return(g)
}


crlg_gs <-  split(crlg_to_plot_df,crlg_to_plot_df$coef) |>
  map(plot_crlg)

write_rds(crlg_gs, snakemake@output$rds)