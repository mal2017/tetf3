library(tidyverse)

mods_fl <- "upstream/final-models.collected-info.tsv.gz"
mods_fl <- snakemake@input$mods
mods <- read_tsv(mods_fl)

g_sharedness <- mods |> 
  dplyr::select(feature.x,feature.y,model,significant_x) |>
  pivot_wider(names_from = model, values_from = significant_x) |>
  filter(male | female) |>
  mutate(status = map2_chr(male,female, ~case_when(
    .x & .y ~ "shared",
    .x & !.y ~ "male-private",
    .y & !.x ~ "female-private",
    !(.y | .x) ~ "n.s.",
    T ~ "wtf"
  ))) |>
  ggplot(aes(status)) +
  geom_bar()


write_rds(g_sharedness,snakemake@output$gg)