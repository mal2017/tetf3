library(phylosignal)
library(parallel)
library(tidyverse)

df_fl <- "results/phylosignal/phylosignal_df.rds"
df_fl <- snakemake@input$df

ps_df <- read_rds(df_fl)

# get the TFs that have at least cutoff significant results between both sexes
# and all phylosignal tests (except lambda, which seems to always be significant
# and seems to yield very oddly-behaved pvals)
cutoff <- 1

hits <- ps_df |>
  filter(score_type == "score") |>
  filter(metric!="Lambda") |> # lambda is basically always significant, so I don't trust it... remove here to reduce false pos
  group_by(sex,TF, coef) |>
  mutate(n_tests_sig = sum(padj < 0.1)) |>
  filter(n_tests_sig >=cutoff) |>
  slice_min(pval, with_ties = F) |>
  ungroup()

write_tsv(hits, snakemake@output$tsv)

# ------------------------------------------------------------------------------
# now generate exemplary correlograms for hits
# ------------------------------------------------------------------------------

to_plot <- filter(hits,TF %in% c("pan","vvl","CG16779","NfI","Unr") | 
                    (padj == min(padj) & n_tests_sig == max(n_tests_sig)))

# get phylosignal objects, previously calculated by another smk rule
ps_fl <- "results/phylosignal/phylosignal.rds"
ps_fl <- snakemake@input$phylosignal

x <- readRDS(ps_fl)

# list of TFs (coefs) to plot, only the 5 we kd'd for know
tfs <- colnames(x$p4d@data)[! colnames(x$p4d@data) %in% c("zad_mean","bm","random","ou")]
names(tfs) <- tfs
tfs <- tfs[to_plot$coef]

# generate correlograms
mc <- getOption("mc.cores", 4)
res <- mclapply(tfs, FUN=function(y) phyloCorrelogram(x$p4d, y))

saveRDS(res, snakemake@output$rds)


