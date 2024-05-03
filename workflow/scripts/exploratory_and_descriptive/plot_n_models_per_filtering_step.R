Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]

mods <- vroom::vroom(mods_path)
  
funs <- list(
  tested = \(x) count(x, model),
  met.OLS.assumptions = \(x) {filter(x, valid) |> count(model)},
  `& reproducible` =  \(x) {filter(x, valid & reproducible) |> count(model)},
  `& significant model fit` = \(x)  {filter(x, valid & reproducible & significant_model) |> count(model)},
  `& signficant host gene expression term` = \(x)  {filter(x, valid & reproducible & significant_model & adj_p.value_anova_x < 0.1) |> count(model)},
  `& no fixed insertion near host gene or its coex. genes` = \(x)  {filter(x, valid & reproducible & significant_model & adj_p.value_anova_x < 0.1 & significant_x) |> count(model)}
)
  
res <- map_df(funs, .f =  function(x) {x(mods)}, .id="subset") 

res$subset <- fct_relevel(res$subset, names(funs))

g_filt <- res |>
  dplyr::rename(sex=model) |>
  mutate(subset = str_wrap(subset,width=20)) |>
  mutate(subset = fct_reorder(subset, n)) |> #pull(subset) %>% .[12]
  ggplot(aes(n, subset,fill=sex)) +
  geom_col(position = "dodge") +
  geom_text(data = \(x) filter(x, subset == subset[which.min(n)]), 
            aes(label=paste0("n=", n), x=300000),
            size=rel(2),
            position = position_dodge(width = 0.75)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) +
  ylab("filtering step") +
  xlab("N TE/gene pairs") +
  scale_fill_grey() +
  theme(legend.position = c(1,0), legend.justification = c("right","bottom"))



write_rds(g_filt,snakemake@output[["gg"]])
