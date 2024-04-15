library(tidyverse)
library(corrr)
library(jsonlite)
library(ggrastr)

cols <- c("model","feature.x","gene_symbol","feature.y",
          "estimate.qnorm","significant_model","significant_x")

# Read in data
main <- ifelse(exists("snakemake"), 
               snakemake@input$mods, 
               "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv(col_select = all_of(cols))

reps <- ifelse(exists("snakemake"), 
               snakemake@input$indep, 
               "upstream/final-models.d2.collected-info.tsv.gz") %>%
  read_tsv(col_select = all_of(cols))

# combine main an indep results
dat <- full_join(main,reps,suffix=c(".main",".rep"),
                 by=c("model","feature.x","gene_symbol","feature.y"))

# correlation between each dataset's 'x' (gene expression coefs) for three filtering approaches
dat <- list(unfiltered = . %>% return,
     replicated = . %>% filter(significant_x.main & significant_x.rep),
     main_data = . %>% filter(significant_x.main)) %>%
  map_df(exec, dat,.id="result_set") %>%
  nest(-model,-result_set) %>%
  mutate(cor.test = map(data,~broom::tidy(cor.test(~ estimate.qnorm.main + estimate.qnorm.rep, data=.x,method = "pearson")))) %>%
  unnest(cor.test)

plot_scatter <- function(data, label) {
  ggplot(data,aes(estimate.qnorm.main, estimate.qnorm.rep)) +
    #ggdensity::geom_hdr_points() +
    #geom_hex() +
    rasterize(geom_point(size=0.1,alpha=0.2),dpi=300) +
    geom_smooth(method="lm") +
    xlab("main results") +
    ylab("independent dataset") +
    annotate("text", -Inf, Inf, label = str_wrap(label,width = 20), hjust = 0, vjust = 1, size=2)
}


res <- dat %>% 
  #filter(result_set!="unfiltered") %>%
  mutate(label = paste0(str_extract(method,"Pearson's|Spearman's",)," r=",round(estimate,digits = 3),"; ",alternative," p=",format.pval(p.value,digits=3))) %>%
  mutate(label=str_replace(label,"=<","<")) %>%
  dplyr::relocate(label) %>%
  mutate(gg= map2(data, label, plot_scatter))

dat %>%
  nest(-result_set,-model) %>%
  nest(-model) %>%
  jsonlite::write_json(snakemake@output$json, prettify=T)

write_rds(res, snakemake@output[["rds"]])
