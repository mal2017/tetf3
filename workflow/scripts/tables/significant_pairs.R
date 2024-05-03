Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(writexl)

x <- read_tsv("upstream/final-models.collected-info.tsv.gz")

x <- filter(x, significant_x)

x <- split(x,x$model)

write_xlsx(x, snakemake@output$xlsx)
