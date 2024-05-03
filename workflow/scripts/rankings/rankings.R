Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)

# Read in data
main <- ifelse(exists("snakemake"), 
               snakemake@input$mods, 
               "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv() |>
  filter(valid)

# Read in data
reps <- ifelse(exists("snakemake"), 
               snakemake@input$reps, 
               "upstream/final-models.collected-info.reps.tsv.gz") %>%
  read_tsv() |>
  filter(valid)

dat <- bind_rows(main=main, reps=reps, .id="dataset")

# get variance explained by x (each gene expression term) as a single value
dat <- dat |>
  mutate(var_exp_x = sumsq_anova_x/total_variance)


# coex scores - abs and signed
dat <- dat %>% 
  dplyr::select(dataset,model,feature.x,gene_symbol,feature.y, estimate.qnorm, var_exp_x, significant_x) %>%
  mutate(across(estimate.qnorm,abs, .names="abs_{.col}")) |>
  mutate(signed_var_exp_x = sign(estimate.qnorm) * var_exp_x)


# i want examinine filtered vs overall results
filtering_funs <- c(sig = . %>% filter(significant_x),
                    nofilt = . %>% return)

rnks <- filtering_funs %>%
  map_df(exec, dat,.id="filtering_approach") %>%
  nest(data=c(-dataset, -model, -filtering_approach))

rnks2 <- rnks %>%
  #head(1) |>
  mutate(stats = map(data, ~summarise(group_by(.x, gene_id = feature.x,gene_symbol),across(where(is.numeric), .fns=c(mean=mean, max=max, sum=sum, min=min), .names="{.fn}_{.col}"),.groups = "drop"))) %>%
  dplyr::select(-data)

res <- rnks2 %>% unnest(stats)

res2 <- res %>%
  pivot_longer(-c(dataset, model, filtering_approach,  gene_id, gene_symbol), 
               names_to = "metric", values_to = "value") %>%
  nest(data=c(gene_id, gene_symbol, value))

# export ------------------------------------------------------------------------
# creates the tags that I use in the snakemake rule to name output files
# then looks in the smk object for them, then writes the appropriate data to them
res2 %>% 
  unite(of, -data) %>% mutate(of = str_replace(of,"\\.","_")) %>%
  deframe() %>%
  iwalk(.f=function(d, f) {
    write_tsv(d, snakemake@output[[f]])
  })
