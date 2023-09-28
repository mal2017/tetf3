# ----------------------------------------------------------------------------
#
#
# not used
#
# 
# ----------------------------------------------------------------------------

library(tidyverse)
library(ape)
library(furrr)

threads <- 2
threads <- snakemake@threads
plan(multisession, workers = threads)

# where x is the in group, y is the 'other'
get_dists <- function(x,y, d=coex.d) {
  d0 <- as.matrix(d)[x,y]
  
  if (identical(x,y)) {
    return(d0[upper.tri(d0,diag = F)])
  } else {
    return(as.vector(d0))
  }
}


# influenced by:
# Influenced by (this tutorial)[https://bioconductor.org/packages/devel/workflows/vignettes/fluentGenomics/inst/doc/fluentGenomics.html]
# get lms from snakemake
lms <- ifelse(exists("snakemake"),snakemake@input[["lms"]],"upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv()

te.names <- unique(lms$feature.y)

tfs <- ifelse(exists("snakemake"),snakemake@input[["tfs"]],"resources/Drosophila_melanogaster_TF.txt") %>%
  read_tsv()

te.classes <- ifelse(exists("snakemake"), snakemake@input$classes, "resources/Tidalbase_Dmel_TE_classifications_2015.txt") %>%
  read_tsv() %>%
  dplyr::select(feature.y = Flybase_name, Class, repClass, repFamily) %>%
  distinct()



# -------------------------- get dist ------------------------------------------
# see Chen et al. BMC Bioinfo 2023
# https://doi.org/10.1186/s12859-023-05161-y
coex.d <- lms %>% 
  filter(feature.x %in% tfs$Ensembl) %>%
  dplyr::select(feature.x, model, feature.y, estimate.qnorm) %>%
  pivot_wider(names_from = c(model,feature.x), values_from = estimate.qnorm) %>%
  column_to_rownames("feature.y") %>%
  t() %>%
  cor(method = "spearman") %>%
  {sqrt(1-abs(.))}


tes <- te.classes %>%
  filter(feature.y  %in% lms$feature.y) %>%
  group_by(superfamily=repFamily,class=repClass) %>%
  summarise(TEs = list(feature.y)) %>%
  mutate(other = map(TEs,.f= function(x) {setdiff(te.names, x)})) %>%
  filter(map_lgl(TEs, ~{length(.x) > 2}) & !is.na(superfamily)) %>%
  mutate(in.dist = map(TEs, .f= function(x) get_dists(x,x))) %>%
  mutate(out.dist = map2(TEs, other, .f= function(x,y) get_dists(x,y))) %>%
  mutate(wc = map2(in.dist, out.dist, ~broom::tidy(ks.test(.x,.y)))) %>%
  unnest(wc)

tes.test <- tes %>%
  dplyr::select(superfamily, class, p.value, in.dist, out.dist) %>%
  mutate(score = map2_dbl(in.dist,out.dist, .f= function(x,y) median(x)-median(y))) %>%
  arrange(score) %>% 
  mutate(rnk = row_number()) %>%
  mutate(label = case_when(p.value < 0.001 ~ "***",
                           p.value < 0.01 ~ "**",
                           p.value < 0.05 ~ "*",
                           T ~ ""))
  
theme_set(theme_classic() +
            theme(strip.text = element_text(size=rel(1.5)),
                  axis.text = element_text(size=rel(1.5)),
                  axis.title = element_text(size=rel(1.5)),
                  plot.caption = element_text(size=rel(1.5))))

g <- tes %>%
  #filter(p.value < 0.05) %>%
  dplyr::select(class,in.dist, out.dist) %>%
  pivot_longer(-c(class, superfamily), names_to = "group",values_to = "distance") %>%
  unnest(distance) %>%
  ggplot(aes(superfamily,distance)) +
  geom_boxplot(aes(fill=group)) +
  geom_text(data = tes.test, aes(label=label, y=1)) +
  scale_fill_grey(start=0.5, end=0.9,labels=c("intra-superfamily","other superfamily"),name="") +
  facet_grid(.~class, scales = "free",space = "free") +
  ylab("dissimilarity") +
  theme(legend.position = "bottom", legend.text = element_text(size=rel(1.5)))

ggsave("../JointGrpMar23/plots/te_superfamily_coex.png",dpi = 300, width=12, height=6)

coex.d %>% as.dist() %>% ape::bionjs() %>% ape::plot.phylo(type = "fan")

set.seed(2022)
boots0 <- lms %>%
  filter(significant_x) %>%
  filter(feature.x %in% tfs$Ensembl) %>%
  left_join(te.classes) %>%
  dplyr::select(grp=repFamily, feature.y) %>%
  distinct() %>%
  group_by(grp) %>%
  filter(n()>=2) %>% # doesn't make sense to get intra group distance for N=1
  nest(data=feature.y)

boots1 <- boots0 %>% 
   #head(n=2) %>% #for debuggin
   group_by(grp) %>%
   summarise(TEs=map(data,pull,feature.y),n_TEs=map_dbl(TEs,length)) %>%
   ungroup() %>%
   mutate(in.dist = future_map_dbl(TEs, get_dists,.options = furrr_options(seed = TRUE))) 

nperm <- 10000
nperm <- snakemake@params[["nperm"]]

bg <- sort(unique(boots1$n_TEs)) %>% set_names(.,as.character(.)) %>% map(get_boot_dists3,n=nperm)
expected <- bg %>% map(mean,na.rm=T) %>% enframe(name="n_TEs",value = "bg.dist") %>%
  unnest(bg.dist) %>%
  mutate(n_TEs = as.numeric(n_TEs))

get_pval <- function(d,n,np=nperm){
  sum(bg[[as.character(n)]] < d)/np
}

boots2 <- boots1 %>%
  mutate(pval = map2_dbl(in.dist,n_TEs,get_pval)) %>%
  mutate(padj = p.adjust(pval,method="BH")) %>%
  left_join(expected)

boots2 %>% filter(pval < 0.05) %>% filter(feature.x %in% tfs$Ensembl)

# ==============================================================================
# below here tests if TEs associated with a single gene
# tend to have more similar coexpression profiles overall than
# randomly selected TEs - this seems like a self fulfilling prophecy
# on the face of it, but points to non-random coexpression scores
# --------------- bootstrapping ------------------------------------------------
# random number consideration: https://www.r-bloggers.com/2020/09/future-1-19-1-making-sure-proper-random-numbers-are-produced-in-parallel-processing/
# set.seed(2022)
# boots0 <- lms %>% 
#   filter(significant_x) %>%
#   filter(feature.x %in% tfs$Ensembl) %>%
#   dplyr::select(feature.x,gene_symbol,feature.y) %>%
#   distinct() %>%
#   group_by(gene_symbol,feature.x) %>%
#   filter(n()>=2) %>% # doesn't make sense to get intra group distance for N=1
#   nest(data=feature.y)
# 
# boots1 <- boots0 %>% 
#   #head(n=2) %>% #for debuggin
#   group_by(gene_symbol,feature.x) %>%
#   summarise(TEs=map(data,pull,feature.y),n_TEs=map_dbl(TEs,length)) %>%
#   ungroup() %>%
#   mutate(in.dist = future_map_dbl(TEs, get_dists,.options = furrr_options(seed = TRUE))) 
# 
# nperm <- 1000
# nperm <- snakemake@params[["nperm"]]
# 
# bg <- 2:length(te.names) %>% set_names(.,as.character(.)) %>% map(get_boot_dists3,n=nperm)
# expected <- bg %>% map(mean,na.rm=T) %>% enframe(name="n_TEs",value = "bg.dist") %>% 
#   unnest(bg.dist) %>%
#   mutate(n_TEs = as.numeric(n_TEs))
# 
# get_pval <- function(d,n,np=nperm){
#   sum(bg[[as.character(n)]] < d)/np
# }
# 
# 
# boots2 <- boots1 %>%
#   mutate(pval = map2_dbl(in.dist,n_TEs,get_pval)) %>%
#   mutate(padj = p.adjust(pval,method="BH"))
# 
# boots2 %>% filter(pval < 0.05) %>% filter(feature.x %in% tfs$Ensembl)

boots2 %>%
  left_join(expected) %>%
  dplyr::select(gene_symbol,in.dist,bg.dist) %>%
  pivot_longer(c(in.dist, bg.dist), names_to = "metric", values_to = "value") %>%
  ggplot(aes(value, fill=metric)) +
  geom_density() 


bg$`2` %>% tibble(dist = .) %>%
  ggplot(aes(dist)) +
  geom_histogram()


write_rds(kmer.dist,snakemake@output[["kmer_dist"]])
write_rds(boots2, snakemake@output[["each_gene"]])
write_rds(bg, snakemake@output[["null_model"]])