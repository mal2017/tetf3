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
fimo_fl <- "results/motifs/fimo_on_tes/pan/fimo.tsv"
fimo_fl <- paste0(snakemake@input$fimo,"/fimo.tsv")
motif_instances <- read_tsv( fimo_fl,comment="#")

dat <- motif_instances |> 
  group_by(label=sequence_name, motif=motif_alt_id) |>
  tally(sort = T) |>
  ungroup() |>
  pivot_wider(names_from = motif, values_from = n,values_fill = 0) |>
  left_join(coex, y=_, by="label") |>
  mutate(across(contains("STREME"),replace_na, 0)) |>
  #mutate(across(contains("STREME"),as.logical)) |>
  mutate(n_motifs = rowSums(across(dplyr::where(is.integer))))

dat2 <- dat |>
  janitor::clean_names() |>
  pivot_longer(-c(label,is_coex),values_to = "n", names_to = "motif")

dat3 <- dat2 |> filter(motif == "n_motifs") |> 
  nest(data=-c(motif)) |>
  mutate(logistic = map(data, ~{glm(is_coex ~ n, data=.x, family=binomial(link='logit'))})) |>
  mutate(logistic.tidy = map(logistic, broom::tidy),
         logistic.aug = map(logistic, broom::augment),
         logistic.glance = map(logistic, broom::glance))

dat4 <- dat3 |>
  mutate(test_data = map2(data, logistic, ~add_predictions(data=tibble(n=seq(0,max(.x$n),0.1)),model=.y, type="response")))

#https://stats.oarc.ucla.edu/r/dae/logit-regression/


write_rds(dat4, snakemake@output$rds)


