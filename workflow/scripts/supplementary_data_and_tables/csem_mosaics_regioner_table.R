library(tidyverse)
library(gt)

regioner_df <- read_rds("results/csem_mosaics/regioner.rds")

regioner_df |>
  dplyr::select(`peak set`=ChIP, 
                `masking approach`=masking, 
                `putative fixed TE insertions` = te_set, 
                `regioneR z-score`=z, 
                pval, 
                `alternative hypothesis`=hypothesis) |>
  grid.table()