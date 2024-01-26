library(tidyverse)

#seas <- Sys.glob("results/motifs/sea_csem_peaks/pan/*/sea.tsv")
seas <- snakemake@input$seas |> paste0("/sea.tsv")

seas <- seas %>% set_names(.,str_extract(.,"ENCS*.+rep\\d+"))

seas <- seas  %>%
  #head() %>%
  map(read_tsv, comment="#")

# remove empty
seas <- seas[map_lgl(seas,~{nrow(.x) >= 1})]

seas <- seas %>% 
  map(~mutate(.x, ID=as.character(ID))) %>% 
  bind_rows(.id="peak_set")

# export to snakemake output as tsv
seas %>% 
  write_tsv(snakemake@output[["tsv"]])



