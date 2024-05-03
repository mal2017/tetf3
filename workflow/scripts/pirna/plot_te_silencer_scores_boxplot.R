
Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)

teregs_fl <- "results/resources/pirna_pathway.tsv"
teregs_fl <- snakemake@input$teregs
teregs <- read_tsv(teregs_fl) |> pull(gene_ID)

# important note - this only uses "valid" models, meaning that some extreme coefficients from
# models that violate OLS assumptions are not used
f_nofilt_max_score <- read_tsv('results/rankings/nofilt_main_female_max_abs_estimate_qnorm.tsv.gz')
f_nofilt_mean_score <- read_tsv('results/rankings/nofilt_main_female_mean_abs_estimate_qnorm.tsv.gz')
f_nofilt_var_exp <- read_tsv('results/rankings/nofilt_main_female_max_var_exp_x.tsv.gz')

m_nofilt_max_score <- read_tsv('results/rankings/nofilt_main_male_max_abs_estimate_qnorm.tsv.gz')
m_nofilt_mean_score <- read_tsv('results/rankings/nofilt_main_male_mean_abs_estimate_qnorm.tsv.gz')
m_nofilt_var_exp <- read_tsv('results/rankings/nofilt_main_male_max_var_exp_x.tsv.gz')


f <- bind_rows(`max abs. score`=f_nofilt_max_score, 
              `mean abs. score`=f_nofilt_mean_score,
              `max var. explained` = f_nofilt_var_exp,
              .id = "score.type") |>
  mutate(prev.reported = gene_id %in% teregs) |>
  mutate(score.type = fct_relevel(score.type,"max var. explained"))
  

m <- bind_rows(`max abs. score`=m_nofilt_max_score, 
               `mean abs. score`=m_nofilt_mean_score,
               `max var. explained` = m_nofilt_var_exp,
               .id = "score.type") |>
  mutate(prev.reported = gene_id %in% teregs) |>
  mutate(score.type = fct_relevel(score.type,"max var. explained"))

plot_score_comparison <- \(x)  x |>
  #filter(score.type == "var. explained") |>
  ggplot(aes(prev.reported,value)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(label.y.npc = c("top","center"),size=2) +
  facet_wrap(~score.type, scales = "free_y") +
  ylab("coexpression metric") +
  xlab("known TE regulator")

g_f <- f |> plot_score_comparison()
g_m <- m |> plot_score_comparison()

write_rds(g_f, snakemake@output$female)
write_rds(g_m, snakemake@output$male)
