library(tidyverse)
library(gt)

regioner_df <- read_rds("results/csem_mosaics/regioner.rds")

x <- regioner_df |>
  filter(!((ChIP=="H3K9Me3") & te_set == "factor_specific")) |>
  mutate(masking = if_else(masking == "mask_heterochromatin","pericentromere",masking)) |>
  mutate(te_set = str_replace(te_set,"_specific","-coexpressed")) |>
  mutate(`observed overlaps`=map_int(regioneR_results,~{.x$numOverlaps$observed})) |>
  mutate(`mean(shuffled set overlaps)`=map_int(regioneR_results,~{.x$numOverlaps$permuted |> mean() |> round()})) |>
  dplyr::select(`peak set`=ChIP, 
                `excluded regions`=masking, 
                `putative fixed TE insertions` = te_set, 
                `observed overlaps`,
                `mean(shuffled set overlaps)`,
                `z-score`=z, 
                pval, 
                `alternative hypothesis`=hypothesis) |>
  filter(`peak set` %in% c("pan","H3K9Me3","vvl","NfI",'CG16779'))



tab <- x |>
  gt(row_group_as_column = T, groupname_col = "peak set") |>
  tab_spanner("regioneR (shuffling analysis)",
              columns = c("observed overlaps",
                          "mean(shuffled set overlaps)",
                          "z-score","pval",
                          "alternative hypothesis"))

dir.create("results/tables/")

gt::gtsave(tab, "results/tables/regioneR_results.table.docx")
