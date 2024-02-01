library(tidyverse)
library(gt)

regioner_df <- read_rds("results/csem_mosaics/regioner.rds")

x <- regioner_df |>
  filter(!((ChIP=="H3K9Me3") & te_set == "factor_specific")) |>
  dplyr::select(`peak set`=ChIP, 
                `masking approach`=masking, 
                `putative fixed TE insertions` = te_set, 
                `regioneR z-score`=z, 
                pval, 
                `alternative hypothesis`=hypothesis) |>
  filter(`peak set` %in% c("pan","H3K9Me3","vvl","NfI",'CG16779')) |>
  gt()


dir.create("results/tables/")

gt::gtsave(x, "results/tables/regioneR_results.table.docx")
