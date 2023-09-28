# ----------------------------------------------------------------------------
#
#
# not used
#
# 
# ----------------------------------------------------------------------------

library(tidyverse)

pid.d <- read_rds("results/te_sequence_similarity/te_pid.rds")

coex.d <- read_rds("results/te_sequence_similarity/coex_dist_mat.rds")

mash.d <- read_rds("results/te_sequence_similarity/te_mash.rds") |> as.matrix()


# make in fam/out fam calls
te.classes <- ifelse(exists("snakemake"), snakemake@input$classes, "resources/Tidalbase_Dmel_TE_classifications_2015.txt") %>%
  read_tsv() %>%
  dplyr::select(feature.y = Flybase_name, Class, repClass, repFamily) %>%
  distinct()

te.classes <- mutate(te.classes, repFamily = ifelse(feature.y == "TART-C","Jockey",repFamily)) %>%
  mutate(repFamily = ifelse(feature.y == "Tc3","TcMar-Tc1",repFamily)) %>%
  mutate(repClass = ifelse(feature.y == "TART-C","LINE",repClass)) %>%
  mutate(repClass = ifelse(feature.y == "Tc3","DNA",repClass))

inout <- te.classes |>
  dplyr::select(x=feature.y, repFamily)

inout <- tidyr::crossing(inout, dplyr::rename(inout, y="x", repFamily.y="repFamily")) |>
  filter(x!=y) |>
  mutate(relationship = if_else(repFamily == repFamily.y,"in.family","out.family"))


inner_join(inout, pid.d) |>
  ggplot(aes(repFamily, pid1, color=relationship)) +
  geom_boxplot()


process_d <- \(x) as_tibble(x, rownames="feature1") |>
  pivot_longer(-feature1,names_to = "feature2",values_to = "d")


df <- full_join(process_d(coex.d),
          process_d(mash.d),by=c("feature1","feature2"), suffix = c(".coex",".mash")) |>
  inner_join(pid.d, by=c(feature1="x",feature2="y")) |>
  filter(feature1 != feature2) |>
  filter(!is.na(d.coex))

ggplot(df, aes(d.mash)) +
  geom_histogram()

# monte carlo #1
closest <- \(x) df$d.coex[which.min(abs(x-df$d.mash))]

sim <- runif(1e6) |>
  tibble(d.sim = _) |>
  mutate(coex.d.closest = map_dbl(d.sim, closest))

cor.test(formula = ~ d.sim + coex.d.closest, data = sim, method="pearson")

# monte carlo #2: uber-simple low cost monte carlo - simulates a perfectly uniform simulation
simple <- seq(0,1,0.001) |>
  tibble(d.sim = _) |>
  mutate(coex.d.closest = map_dbl(d.sim, closest))

ggplot(sim, aes(d.sim)) +
  geom_histogram()

ggplot(sim, aes(coex.d.closest)) +
  geom_histogram()

ggplot(df, aes(d.coex)) +
  geom_histogram()

ggplot(sim,aes(d.sim, coex.d.closest)) +
  geom_smooth(method="lm")

ggdensity::geom_hdr_lines()# +
#ggdensity::geom_hdr(method="histogram")


# ==============================================================================
# bootstrap
boot_res <- df |>
  bootstrap(n=1000) |>
  mutate(lmr = map(strap, ~{lm(d.coex~d.mash, data=.x)})) |>
  mutate(gl = map(lmr, broom::glance), tidy = map(lmr, ~filter(broom::tidy(.x),term=="d.mash"))) |> # , aug=map(lmr, broom::augment)
  unnest(c(tidy, gl), names_sep = "_")

ggplot(boot_res, aes(gl_p.value)) +
  geom_histogram() +
  scale_x_log10()

ggplot(boot_res, aes(tidy_estimate)) +
  geom_histogram()

ggplot(boot_res, aes(gl_r.squared)) +
  geom_histogram()


df |>
  ggplot(aes(perc_ident, d.mash)) +
    geom_point()
  
df |>
  ggplot(aes(d.mash, d.coex)) +
  geom_point() +
  ggdensity::geom_hdr_points() +
  ggpubr::stat_cor(method = "spearman")

df |>
  mutate(bin = cut_interval(d.mash,n=5)) |>
  ggplot(aes(bin, d.coex)) +
  geom_violin() +
  ggpubr::stat_compare_means()
  #geom_jitter(width = 0.2, alpha=0.1)
