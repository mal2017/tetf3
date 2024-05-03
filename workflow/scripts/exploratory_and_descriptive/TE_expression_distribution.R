Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)
library(SummarizedExperiment)

se_fl <- snakemake@input$se
se <- read_rds(se_fl)

# [str_detect(rownames(se),"FBgn", negate = T),]
res <- se %>%
  assay() %>%
  #.[,1] %>%
  as_tibble(rownames="feature") %>%
  pivot_longer(-feature) %>%
  mutate(feature.type = ifelse(str_detect(feature, "FBgn", negate = T),"TE","gene"))

write_rds(res,snakemake@output[["rds"]])

#ggplot(aes(value, fill=feature.type)) +
#  geom_histogram() +
#  scale_x_log10() +
#  scale_y_log10() +
#  facet_wrap(~feature.type, scales="free_y", )
