library(tidyverse)

pkgs <- read_lines(pipe("grep library workflow/scripts/*/*R")) |>
    str_extract_all("(?<=library\\().+?(?=\\))") |>
    flatten() |>
    unlist() |>
    str_remove_all("\"") |>
    unique()

walk(pkgs, library, character.only = TRUE)

write_lines(capture.output(sessionInfo()), "sessionInfo.txt")
