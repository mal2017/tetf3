library(tidyverse)
library(clusterProfiler)
library(ggnewscale)

# mean abs
res <- ifelse(exists("snakemake"), snakemake@input$rds, "results/enrichment/sig_main_male_mean_abs_estimate_qnorm.gg_gsea.rds") %>% read_rds
res <- ifelse(exists("snakemake"), snakemake@input$rds, "results/enrichment/sig_main_female_mean_abs_estimate_qnorm.gg_gsea.rds") %>% read_rds

res <- ifelse(exists("snakemake"), snakemake@input$rds, "results/enrichment/nofilt_main_male_mean_estimate_qnorm.gg_gsea.rds") %>% read_rds
res <- ifelse(exists("snakemake"), snakemake@input$rds, "results/enrichment/nofilt_main_female_mean_estimate_qnorm.gg_gsea.rds") %>% read_rds

res <- ifelse(exists("snakemake"), snakemake@input$rds, "results/enrichment/nofilt_main_male_mean_estimate_qnorm.go_gsea.rds") %>% read_rds
res <- ifelse(exists("snakemake"), snakemake@input$rds, "results/enrichment/nofilt_main_female_mean_estimate_qnorm.go_gsea.rds") %>% read_rds

res <- ifelse(exists("snakemake"), snakemake@input$rds, "results/enrichment/sig_main_male_max_abs_estimate_qnorm.go_gsea.rds") %>% read_rds
res <- ifelse(exists("snakemake"), snakemake@input$rds, "results/enrichment/sig_main_female_max_abs_estimate_qnorm.go_gsea.rds") %>% read_rds


res <- ifelse(exists("snakemake"), snakemake@input$rds, "results/enrichment/nofilt_main_male_sum_estimate_qnorm.gg_gsea.rds") %>% read_rds
res <- ifelse(exists("snakemake"), snakemake@input$rds, "results/enrichment/nofilt_main_female_sum_var_exp_x.gg_gsea.rds") %>% read_rds

res %>% dotplot()

ranking_selector_funs <- list(main_male = . %>% filter(filtering_approach == "main" & dataset == "main" & model == "male" & metric == "mean_abs_estimate.qnorm"),
     main_female = . %>% filter(filtering_approach == "main" & dataset == "main" & model == "female" & metric == "mean_abs_estimate.qnorm"),
    indep_male = . %>% filter(filtering_approach == "replicate" & dataset == "reps" & model == "male" & metric == "mean_abs_estimate.qnorm"), 
     indep_female = . %>% filter(filtering_approach == "replicate" & dataset == "reps" & model == "female" & metric == "mean_abs_estimate.qnorm")) 

# sanity check
#pull(test,gsea)[[1]] %>% as_tibble() %>% filter(str_detect(ID,"C2H2")) %>% pull(core_enrichment) %>% str_count('/')
make_dotplot <- function(.x,.y) {
  dotplot(.x,showCategory=50, x="NES", size="str_count(core_enrichment,'/')") + 
  theme_linedraw() +
  labs(size="leading edge members") +
    ggtitle(.y)
}

ranking_selector_funs %>%
  map_df(exec,res, .id="label") %>%
  mutate(gg_dot = map2(gsea, label, make_dotplot)) %>%
  pull(gg_dot)




res %>%
  filter(metric == "mean_estimate.qnorm") %>%
  dplyr::select(model, filtering_approach, dataset, gsea.tidy) %>%
  unnest(gsea.tidy) %>%
  filter(ID == "ZAD_ZNF")


x <- read_tsv("upstream/final-models.collected-info.tsv.gz")

x %>%
  filter(feature.y=="1360" & gene_symbol == "pan") %>% glimpse()

tfs <- read_tsv("resources/Drosophila_melanogaster_TF.txt")


x %>% filter(significant_x) %>%
  dplyr::select(gene_symbol,feature.y) %>%
  distinct() %>%
  count(gene_symbol,sort = T) %>%
  filter(gene_symbol %in% tfs$Symbol) %>%
  print(n=Inf)
  group_by(feature.y, gene_symbol) %>%
  filter(n()==1)
