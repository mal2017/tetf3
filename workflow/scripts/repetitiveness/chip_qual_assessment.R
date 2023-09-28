library(tidyverse)


# used for matching inputs to IPs
ss_fl <- "upstream/sample_table_encode.csv"
ss_fl <- snakemake@input$ss
sample_sheet <- read_csv(ss_fl) # used for filtering later

# quality info from deeptools
chip_qual_fl <- "upstream/all_fingerprint.metrics.txt"
chip_qual_fl <- snakemake@input$fingerprint
chip_qual <- read_tsv(chip_qual_fl)


# ------------------------------------------------------------------------------
# get chip qc metrics for basic filtering
# ------------------------------------------------------------------------------
# see deeptools metrics for info

chip_qual_df <- mutate(chip_qual, sample =str_extract(Sample,"(?<=masked_bwa_mem2\\/).+(?=\\.srt.bam)")) |>
  dplyr::select(-Sample) |>
  relocate(sample)

input_df <- filter(chip_qual_df, str_detect(sample,"input|INPUT")) |>
  group_by(sample) |>
  slice_head(n=1) |>
  ungroup() |>
  dplyr::select(sample, `Synthetic JS Distance`, AUC, `X-intercept`, `Elbow Point`)

chip_qual_df <- filter(chip_qual_df, !str_detect(sample,"input|INPUT"))

# note encode assigns same input to multiple samples...
chip_qual_df <- sample_sheet |>
  dplyr::select(sample=sample_name, input, target) |>
  left_join(chip_qual_df,y=_) |>
  left_join(input_df, by=c(input="sample"), relationship="many-to-many", suffix=c(".IP",".input") )

chip_qual_df <- chip_qual_df |> 
  mutate(elbow_diff = `Elbow Point.IP`-`Elbow Point.input`) |>
  mutate(AUC_diff = AUC.IP-AUC.input)

chip_qual_df <- chip_qual_df |>
  mutate(signed.js = sign(elbow_diff) * `Synthetic JS Distance.IP`)


write_rds(chip_qual_df, snakemake@output$rds)