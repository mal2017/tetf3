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

#res |>
#  ggplot(aes(subset, n, fill=model)) +
#  geom_col(position="dodge") +
#  theme(axis.text.x = element_text(angle=45, hjust=1)) 

write_rds(res,snakemake@output[["rds"]])
