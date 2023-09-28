library(tidyverse)

rds_fl <- "results/pirna/pirna_genes_in_lms.rds"
rds_fl <- snakemake@input$rds

piRNA_genes_in_lms <- read_rds(rds_fl)

pirna_prop <- piRNA_genes_in_lms |>
  group_by(model,class,gene_symbol) |>
  summarise(coex = any(significant_x)) |>
  group_by(model, class) |>
  summarise(n_coex = sum(coex),n_not_coex=sum(!coex)) |>
  arrange(model,desc(class)) |>
  ungroup()

pirna_fish <- pirna_prop |> 
  nest(data=-model) |>
  mutate(data = map(data,column_to_rownames,"class")) |>
  mutate(fish = map(data, ~broom::tidy(fisher.test(.x)))) |>
  unnest(fish)

pirna_fish_label <- pirna_fish |>
  mutate(label = sprintf("\t\tFisher's Exact\n\t\t%s p=%s\n\t\tOdds ratio=%s", alternative,format.pval(p.value,1),round(estimate,2)))

write_rds(pirna_fish_label, snakemake@output$rds)
