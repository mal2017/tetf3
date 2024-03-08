library(tidyverse)

supp_fl <- "resources/pericent_enriched_pan_chips_by_inspection.json"
supp_fl <- snakemake@input$json
supporting <- jsonlite::read_json(supp_fl,simplifyVector = T)

repet <- "results/repetitiveness/chip_repetitiveness.rds"
repet <- snakemake@input$rds

repet <- read_rds(repet) |>
  filter(str_detect(target,"pan")) |>
  mutate(target = fct_reorder(target,ratio.te))

g <- repet |>
  mutate(visual.pericentromeric.enrichment=sample %in% supporting) |>
  ggplot(aes(visual.pericentromeric.enrichment,ratio.te)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()

write_rds(g,snakemake@output$gg)