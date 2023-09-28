library(tidyverse)

# ------------------------------------------------------------------------------
# get data
# ------------------------------------------------------------------------------

# used for matching inputs to IPs
ss_fl <- "upstream/sample_table_encode.csv"
ss_fl <- snakemake@input$ss
sample_sheet <- read_csv(ss_fl) # used for filtering later

te_fa <- "upstream/Tidalbase_transposon_sequence.dmel.fixed.sensechecked.fa.gz"
te_fa <- snakemake@input$te_fa
tes <- names(Biostrings::readDNAStringSet(te_fa))

idxstats_fl <- "upstream/masked_idxstats.txt"
idxstats_fl <- snakemake@input$idxstats
idxstats <- read_tsv(idxstats_fl, col_names = c("seqname","seqlen","n_mapped","n_unmapped","file"))

qc_fl <- "results/repetitiveness/chip_qual_assessment.rds"
qc_fl <-  snakemake@input$qc
chip_qual_df <- read_rds(qc_fl)
# just filter for samples that have some detectable IP enrichment
# see deeptools metrics docs for details
acceptable_chips <- chip_qual_df |> filter(AUC_diff < 0 & signed.js > 0) |> pull(sample)

# ------------------------------------------------------------------------------
# te mapping ratio analysis
# ------------------------------------------------------------------------------

# sum read counts by subject seq type: te or genome
df <- idxstats |>
  group_by(file, mapping = if_else(seqname %in% tes,"TE","genome")) |>
  summarise(n_mapped=sum(n_mapped), .groups = "drop")

# clean up the df
df <- df |>
  mutate(sample =str_extract(file,"(?<=masked\\/).+(?=\\.idxstats.txt)")) |>
  relocate(sample) |>
  dplyr::select(-file)

# marry up ChIPs and IPs
df <- sample_sheet |>
  dplyr::select(sample=sample_name, input) |>
  filter(!is.na(input)) |>
  left_join(df, by= c("sample")) |>
  left_join(df, by=c(input="sample",mapping="mapping"), suffix = c(".IP",".input")) |>
  dplyr::select(sample, input, mapping, n_mapped.IP, n_mapped.input)

# get ratio of ratios
df2 <- df |>
  pivot_wider(names_from = c(mapping), values_from = c(n_mapped.IP, n_mapped.input)) |>
  mutate(ratio=(n_mapped.IP_TE/n_mapped.IP_genome)/(n_mapped.input_TE/n_mapped.input_genome)) |>
  dplyr::select(sample,ratio) |>
  left_join(sample_sheet, by=c(sample="sample_name")) |>
  dplyr::select(sample, target, ratio.te = ratio) |>
  left_join(df, by="sample")

fish_df <- df2 |>
  filter(sample %in% acceptable_chips) |>
  filter(str_detect(target,regex("H3K|pan|egg|setdb|var|gro|arm|nej",ignore_case=T))) |>
  nest(data=c(mapping,n_mapped.IP, n_mapped.input)) |>
  mutate(data=map(data,column_to_rownames,"mapping")) |>
  mutate(fish = map(data, fisher.test, alternative = "greater")) |>
  mutate(fish.tidy = map(fish,broom::tidy)) |>
  unnest(fish.tidy)

write_rds(fish_df, snakemake@output$rds)


