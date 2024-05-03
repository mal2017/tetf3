Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)

# Read in data
main <- ifelse(exists("snakemake"), 
               snakemake@input$mods, 
               "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv()

# Read in data
reps <- ifelse(exists("snakemake"), 
               snakemake@input$reps, 
               "upstream/final-models.collected-info.reps.tsv.gz") %>%
  read_tsv()

dat <- bind_rows(main=main, reps=reps, .id="dataset")

# coex scores - abs and signed
coex <- dat %>% 
  dplyr::select(dataset,model,feature.x,gene_symbol,feature.y,contains("estimate.qnorm"), significant_x) %>%
  mutate(across(starts_with("estimate"),abs, .names="abs_{.col}")) %>%
  dplyr::select(dataset,model,feature.x,gene_symbol,feature.y, contains("estimate"), significant_x)

# i want examinine filtered vs overall results
filtering_funs <- c(sig = . %>% filter(significant_x),
                    nofilt = . %>% return)

coex <- filtering_funs %>%
  map_df(exec, coex,.id="filtering_approach") %>%
  nest(-dataset, -model, -filtering_approach)

res2 <- coex %>%
  mutate(stats = map(data, ~summarise(group_by(.x, gene_id = feature.x,gene_symbol),across(where(is.numeric), .fns=c(mean=mean, max=max, sum=sum), .names="{.fn}_{.col}", na.rm=T),.groups = "drop"))) %>%
  dplyr::select(-data)

res2 <- res2 %>%unnest(stats)

res2 <- res2 %>%
  pivot_longer(-c(dataset, model, filtering_approach,  gene_id, gene_symbol), 
               names_to = "metric", values_to = "value") %>%
  nest(gene_id, gene_symbol, value)

# export ------------------------------------------------------------------------
# creates the tags that I use in the snakemake rule to name output files
# then looks in the smk object for them, then writes the appropriate data to them
res2 %>% 
  unite(of, -data) %>% mutate(of = str_replace(of,"\\.","_")) %>%
  deframe() %>%
  iwalk(.f=function(d, f) {
    write_tsv(d, snakemake@output[[f]])
  })
