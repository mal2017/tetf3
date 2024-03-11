library(janitor)
library(modelr)
library(tidyverse)

tf <- 'pan'
tf <- snakemake@params$tf

# get info about pan coexpression
mods_fl <- "~/work/tetf3/upstream/final-models.collected-info.tsv.gz"
mods_fl <- snakemake@input$coex
mods <- read_tsv(mods_fl)

coex <- filter(mods, gene_symbol == tf)

coex <- coex |>
  dplyr::select(label=feature.y, sex=model, significant_x, coex.score = estimate.qnorm) |>
  group_by(label) |>
  summarise(is.coex = any(significant_x))

# get info about motif presence.
fimo_fl <- "results/motifs/fimo_on_tes/denovo/pan/fimo.tsv"
fimo_fl <- paste0(snakemake@input$fimo,"/fimo.tsv")
motif_instances <- read_tsv( fimo_fl,comment="#")

dat <- motif_instances |> 
  group_by(label=sequence_name, motif=motif_id) |>
  tally(sort = T) |>
  ungroup() |>
  pivot_wider(names_from = motif, values_from = n,values_fill = 0) |>
  left_join(coex, y=_, by="label") |>
  mutate(n_motifs = rowSums(across(dplyr::where(is.integer))))

dat2 <- dat |>
  pivot_longer(-c(label,is.coex),values_to = "n", names_to = "motif")

dat3 <- dat2 |> #filter(motif == "n_motifs") |> 
  mutate(n=replace_na(n,0)) |>
  nest(data=-c(motif)) #|>
#  mutate(logistic = map(data, ~{glm(is.coex ~ n, data=.x, family=binomial(link='logit'))})) |>
#  mutate(logistic.tidy = map(logistic, broom::tidy))

#dat4 <- dat3 |>
#  mutate(smoothed_data = map2(data, logistic, ~add_predictions(data=tibble(n=seq(0,max(.x$n),0.1)),model=.y, type="response")))


dat4 <- dat3 |>
  mutate(wilcox = map(data,~wilcox.test(n~is.coex,data=.x))) |>
  mutate(wilcox.tidy = map(wilcox,broom::tidy))

# make plots
dat5 <- dat4 |>
  mutate(motif = paste0("denovo::",motif)) |>
  unnest(wilcox.tidy) |>
  mutate(g_boxplot = pmap(list(data,p.value,motif),function(.x,.y,.z){
    ggplot(.x,aes(is.coex, n)) +
      geom_boxplot() +
      ylab(sprintf("N %s hits",.z)) +
      xlab("coexpressed") +
      annotate("text",x=-Inf,y=Inf,label=sprintf("Wilcoxon\np=%s",format.pval(.y,eps = 0.005)),hjust=0,vjust=1)
  }))

#https://stats.oarc.ucla.edu/r/dae/logit-regression/

write_rds(dat5, snakemake@output$rds)


