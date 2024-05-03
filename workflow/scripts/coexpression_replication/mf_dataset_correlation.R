Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


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

# combine main an indep results
dat <- full_join(filter(main,model=="male"),filter(main,model=="female"),suffix=c(".male",".female"),
                 by=c("feature.x","gene_symbol","feature.y"))

# correlation between each dataset's 'x' (gene expression coefs) for three filtering approaches
dat <- list(unfiltered = . %>% return,
     replicated = . %>% filter(significant_x.male & significant_x.female),
     male_data = . %>% filter(significant_x.male),
     female_data = . %>% filter(significant_x.female)) %>%
  map_df(exec, dat,.id="result_set") %>%
  nest(-result_set) %>%
  mutate(cor.test = map(data,~broom::tidy(cor.test(~ estimate.qnorm.male + estimate.qnorm.female, data=.x,method = "pearson")))) %>%
  unnest(cor.test)

plot_scatter <- function(data, label) {
  ggplot(data,aes(estimate.qnorm.male, estimate.qnorm.female)) +
    #ggdensity::geom_hdr_points() +
    #geom_hex() +
    rasterize(geom_point(size=0.1,alpha=0.2),dpi=300) +
    geom_smooth(method="lm") +
    xlab("male results") +
    ylab("female results") +
    annotate("text", -Inf, Inf, label = str_wrap(label,width = 20), hjust = 0, vjust = 1, size=2)
}

res <- dat %>% 
  #filter(result_set!="unfiltered") %>%
  mutate(label = paste0(str_extract(method,"Pearson's|Spearman's",)," r=",round(estimate,digits = 3),"; ",alternative," p=",format.pval(p.value,digits=3))) %>%
  mutate(label=str_replace(label,"=<","<")) %>%
  dplyr::relocate(label) %>%
  mutate(gg= map2(data, label, plot_scatter))

dat %>%
  nest(-result_set) %>%
  jsonlite::write_json(snakemake@output$json, prettify=T)

write_rds(res, snakemake@output[["rds"]])
