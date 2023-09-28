library(phylosignal)
library(parallel)
library(tidyverse)

df_fl <- "results/phylosignal/phylosignal_df.rds"
df_fl <- snakemake@input$df

ps_df <- read_rds(df_fl)

 cutoff <- 2

hits <- ps_df |>
   filter(score_type == "score") |>
   filter(metric!="Lambda") |>
      group_by(sex,TF, coef) |>
   summarise(n_tests_passed = sum(padj < 0.1),  .groups="drop") |>
   group_by(TF) |>
   filter(sum(n_tests_passed >= cutoff) == 2) |>
   ungroup()

#hits <- ps_df |>   filter(score_type == "score") |> filter(metric == "I" & pval < 0.01) |>
#  count(TF) |>
#  filter(n==2)

hits <- ps_df |> filter(TF %in% hits$TF) |> pull(coef) |> unique()

ps_fl <- "results/phylosignal/phylosignal.rds"
ps_fl <- snakemake@input$phylosignal

x <- readRDS(ps_fl)

tfs <- colnames(x$p4d@data)[! colnames(x$p4d@data) %in% c("zad_mean","bm","random","ou")]

names(tfs) <- tfs

# only hi-conf m/f hits for now. these take forever to generate.
tfs <- tfs[hits]


mc <- getOption("mc.cores", 4)

res <- mclapply(tfs, FUN=function(y) phyloCorrelogram(x$p4d, y))

saveRDS(res, snakemake@output$rds)