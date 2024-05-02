library(tidyverse)

# get te classes
te.classes.fl <- "~/work/tetf3/resources/Tidalbase_Dmel_TE_classifications_2015.txt"
te.classes.fl <- snakemake@input$te_classes

te.classes <- read_tsv(te.classes.fl) %>%
  dplyr::select(feature.y = Flybase_name, Class, repClass, repFamily) %>%
  distinct()

# fix some missing info
te.classes <- mutate(te.classes, repFamily = ifelse(feature.y == "TART-C","Jockey",repFamily)) %>%
  mutate(repFamily = ifelse(feature.y == "Tc3","TcMar-Tc1",repFamily)) %>%
  mutate(repFamily = ifelse(feature.y %in% c("TLD2_LTR","Stalker3T"),"Gypsy",repFamily)) %>%
  mutate(repClass = ifelse(feature.y == "TART-C","LINE",repClass)) %>%
  mutate(repClass = ifelse(feature.y == "Tc3","DNA",repClass)) |>
  mutate(repClass = ifelse(feature.y %in% c("TLD2_LTR","Stalker3T"),"LTR",repClass)) 

coex_dist_fl <- "results/te_sequence_similarity/coex_dist_df.rds"
coex_dist_fl <- snakemake@input$coex_dist

seq_dist_fl <- "results/te_sequence_similarity/te_sketch_dist.rds"
seq_dist_fl <- snakemake@input$seq_dist

coex_dist <- read_rds(coex_dist_fl)
seq_dist <- read_rds(seq_dist_fl)

# make into a tibble
seq_dist_df <- seq_dist |> 
  as.matrix() |>
  as_tibble(rownames="te1") |>
  pivot_longer(-te1, names_to = "te2",values_to = "seq_distance") |>
  filter(te1!=te2)


# note that this has both male and femal in it, so is 2x longer
coex_dist_df <- coex_dist |>
  dplyr::select(model,tab) |>
  unnest(tab) |>
  dplyr::rename(coex_dist="dissimilarity") |>
  filter(te1!=te2)


dat <- full_join(coex_dist_df, seq_dist_df, by=c("te1","te2")) |>
  left_join(dplyr::select(te.classes,feature.y,repFamily, repClass), by=c("te1"="feature.y")) |>
  left_join(dplyr::select(te.classes,feature.y,repFamily, repClass), by=c("te2"="feature.y"),suffix=c(".1",".2")) 

write_rds(dat,snakemake@output$rds)