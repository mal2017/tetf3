Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]
mods <- read_tsv(mods_path) #%>% filter(valid)

#pirna_path <- "results/resources/pirna_pathway.tsv"
pirna_path <- snakemake@input[["pirna"]]
pirna_tbl <- read_tsv(pirna_path)

# ------------------------------------------------------------------------------
# get clean table of piRNA genes in the lms
# ------------------------------------------------------------------------------
res <- mods %>%
  mutate(class = ifelse(feature.x %in% pirna_tbl$gene_ID,"ovary TE silencer","other")) %>%
  dplyr::select(model, gene_symbol, feature.y, class, estimate.qnorm, valid, significant_model, significant_x)


# ------------------------------------------------------------------------------
# fisher's asking if te silencers are enriched among the genes with >0 TE corrs
# ------------------------------------------------------------------------------
pirna_prop <- res |>
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


# ------------------------------------------------------------------------------
# plot
# ------------------------------------------------------------------------------

g_te_silencer_overrep_among_te_coex_genes <- pirna_fish_label |> 
  mutate(data = map(data, ~rownames_to_column(.x, "class"))) |>
  unnest(data) |>
  mutate(proportion = n_coex/(n_coex + n_not_coex)) |>
  ggplot(aes(str_wrap(class,10),proportion)) +
  geom_col() +
  geom_text(data=pirna_fish_label,aes(x=-Inf, y=1.1, label=label), hjust=0, size=rel(1.5)) +
  facet_wrap(~model) +
  coord_cartesian(ylim=c(0,1.3)) +
  xlab("gene type") +
  ylab("proportion >=1 significant TE correlation")

write_rds(g_te_silencer_overrep_among_te_coex_genes, snakemake@output[["gg"]])
saveRDS(res, snakemake@output[["rds"]])