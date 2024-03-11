library(readxl)
library(writexl)
library(tidyverse)

fls <- Sys.glob("results/figures/data/*.xlsx")

names(fls) <- str_extract(fls,"(?<=data\\/)figure.+?(?=\\.)")

finalsheets <- list()

for (f in names(fls)) {
  sheetnames <- readxl::excel_sheets(fls[[f]])
  for (s in sheetnames) {
    newsheetname <- paste(f,s)
    finalsheets[[newsheetname]] <- read_xlsx(fls[[f]],sheet = s,)
  }

}

write_xlsx(finalsheets,snakemake@output$xlsx)
