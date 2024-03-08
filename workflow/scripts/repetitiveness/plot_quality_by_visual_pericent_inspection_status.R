library(tidyverse)

# ------------------------------------------------------------------------------
# quality vs repetitiveness
# ------------------------------------------------------------------------------
supp_fl <- "resources/pericent_enriched_pan_chips_by_inspection.json"
supp_fl <- snakemake@input$json
supporting <- jsonlite::read_json(supp_fl,simplifyVector = T)

repet <- "results/repetitiveness/chip_repetitiveness.rds"
repet <- snakemake@input$rep

qc_fl <- "results/repetitiveness/chip_qual_assessment.rds"
qc_fl <- snakemake@input$qc

qc_df0 <- read_rds(qc_fl)

qc_df <- filter(qc_df0, !str_detect(sample,"input") & str_detect(sample,"pan")) |> 
  mutate(experiment = str_extract(sample,"ENCSR.+(?=_rep)")) |>
  mutate(library = str_extract(sample,"pan.+rep\\d+")) |>
  mutate(visual.pericentromeric.enrichment=library %in% supporting) |>
  dplyr::select(library,visual.pericentromeric.enrichment,c("JS Distance"))

g <- qc_df |>
  ggplot(aes(visual.pericentromeric.enrichment,`JS Distance`)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2)

write_rds(g,snakemake@output$gg)