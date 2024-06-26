Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)

#mods_fl <- "upstream/final-models.collected-info.tsv.gz"
mods_fl <- snakemake@input[["mods"]]

mods <- read_tsv(mods_fl)

dat <- mods %>%
  filter(significant_model) %>%
  mutate(sumsq_anova_remaining = explained_variance*total_variance) %>%
  dplyr::select(model,feature.x,feature.y,total_variance,contains("sumsq"), sumsq_anova_remaining) %>%
  mutate(across(contains("sumsq"),~{100*.x/total_variance})) %>%
  pivot_longer(contains("sumsq"), names_to = "coef", values_to = "var.explained") %>%
  mutate(coef = str_remove(coef,"sumsq_anova_"))


# rename for clear x axis labels
dat <- dat |> mutate(coef = case_when(coef == "overlap.coex.gene" ~ "overlap with highly coexpressed gene",
                               coef == "overlap" ~ "direct gene/TE overlap",
                               coef == "x" ~ "gene coexpression term",
                               coef == "wolbachia" ~ "wolbachia infection",
                               coef == "scaled.copies.y" ~ "TE copy number",
                               T ~ coef))

dat <- dat %>% 
  mutate(coef = fct_reorder(coef,var.explained,.fun = function(.x) {median(.x,na.rm=T)})) %>%
  mutate(var.explained = var.explained)

g <-  dat %>% ggplot(aes(coef,var.explained,fill=model)) +
  ggrastr::rasterise(geom_point(color="gray",position = position_jitterdodge(jitter.width=0.2),size=0.1,alpha=0.05),dpi=300) +
  geom_boxplot(outlier.shape=NA) +
  ylab("% explained variance")

saveRDS(g,snakemake@output[["rds"]])