library(tidyverse)
library(jsonlite)

x <- read_json("results/stats/descriptive_lms.json",simplifyVector = F)

x |> enframe() |> unnest(value)
