Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(jsonlite)

x <- read_json("results/stats/descriptive_lms.json",simplifyVector = F)

x |> enframe() |> unnest(value)
