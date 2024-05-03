Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


# read in results of supercell analysis
supercell <- readRDS("results/calderon22/fca_reanalysis_supercell.rds")

GE<-supercell$GE
SC<-supercell$SC

rm(supercell);gc()

library(plotgardener)
library(tidyverse)
library(SuperCell)
library(patchwork)

tfs <- read_tsv("resources/Drosophila_melanogaster_TF.txt")

tes <- jsonlite::read_json("upstream/te_element_lookup.json") %>%
  names()

# read in TE TF correlations we identifed from the supercell results
tf_te_correlations <- read_rds("results/calderon22/fca_reanalysis_correlations.rds") |>
  dplyr::select(lineage,res=res.spqn) |>
  unnest(res)

# ------------------------------------------------------------------------------
# plot exemplary and control results from supercell analysis
# ------------------------------------------------------------------------------

plot_indiv_relationship <- function(x, y, caption="") {
  require(SuperCell)
  supercell_GeneGenePlot(GE,
                         gene_x = x, 
                         gene_y = y,
                         supercell_size = SC$supercell_size,
                         color.use = "black") -> res
  
  g <- res$p$data |> 
    as_tibble() |>
    ggplot(aes(x, y)) +
    #ggdensity::geom_hdr_points(method="mvnorm",aes(size=size)) +
    ggrastr::rasterise(geom_point(size=0.5,alpha=0.5),dpi=300) +
    xlab(x) +
    ylab(y) +
    labs(title=sprintf("%s vs %s",x,y),subtitle = sprintf("raw weighted corr:%s",round(res$w.cor[[1]],digits = 3)),caption = caption) +
    geom_smooth(method="lm",se=F,color="red",linetype="dashed",linewidth=1) +
    scale_size_continuous(name="barcodes/supercell") +
    scale_color_viridis_d(name="density")
  
  return(g)
}

# https://elifesciences.org/articles/68573
# summarizes refs for deposition of piwi in embryo - espec piwi
# which is also deposited in embryonic soma
# see franz et al 2017 (probing canonicity of wnt pathway) for selection
# of control genes
g_pan_te_expemplary1 <- plot_indiv_relationship("pan","invader2")
g_pan_te_exemplary <- plot_indiv_relationship("pan","gypsy10")
g_pan_poscon <- plot_indiv_relationship("pan","Toll-7")
g_pan_negcon <- plot_indiv_relationship("pan","CG4115")

g_pan_exemplary <- g_pan_poscon + g_pan_negcon + g_pan_te_expemplary1  + g_pan_te_exemplary + 
  plot_layout(nrow=1,guides = "collect") & 
  theme(legend.position = "bottom")


# https://www.mdpi.com/2077-0383/8/4/560
#g_unr_myc <- plot_indiv_relationship("Unr","Myc")


# http://flybase.org/reports/FBgn0263352.htm#interactions
# mihailovic 2012, RNA, also JIL-1, Nup153
g_unr_1 <- plot_indiv_relationship("Unr","brat")
g_unr_2 <- plot_indiv_relationship("Unr","pasha")
g_unr_3 <- plot_indiv_relationship("Unr","jockey")
g_unr_4 <- plot_indiv_relationship("Unr","gypsy11")

g_unr_exemplary <- g_unr_1 + g_unr_2 + g_unr_3 + g_unr_4 +
  plot_layout(nrow=1,guides = "collect") & 
  theme(legend.position = "bottom")



# asking if pan is highly correlated with TEs in general
g_pan_highly_corr_with_tes <- 
  tf_te_correlations %>%
  filter(feature %in% tfs$Symbol & y %in% tes) |>
  dplyr::select(feature,y,coef,p,padj) |>
  distinct() |>
  mutate(feature2 = if_else(feature%in%c("pan","Unr","vvl","NfI","CG16779"),feature,"other")) |>
  mutate(feature2 = fct_reorder(feature2,coef)) |>
  mutate(feature2 = fct_relevel(feature2,"other")) |>
  ggplot(aes(feature2,coef)) +    
  geom_boxplot() +
  ggpubr::stat_compare_means(ref="other",size=2) +
  xlab("") + ylab("weighted correlation")



# plotting page 1 --------------------------------------------------------------------
theme_set(theme_classic() + theme(text=element_text(size=unit(7,"pt"))))

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 8.5, height = 11)

pageCreate(width =8.5, height = 11, default.units = "inches", showGuides = interactive())
figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

plotGG(g_pan_exemplary,x=0.4,y=.5,width=7.6,height=2.5)
plotText(label = "A", x = 0.5, y = .5)
plotText(label = "B", x = 2.4, y = .5)
plotText(label = "C", x = 4.15, y = .5)
plotText(label = "D", x = 6.15, y = .5)


plotGG(g_unr_exemplary,x=0.4,y=3.05,width=7.6,height=2.5)
plotText(label = "E", x = 0.5, y = 3.05)
plotText(label = "F", x = 2.4, y = 3.05)
plotText(label = "G", x = 4.15, y = 3.05)
plotText(label = "H", x = 6.15, y = 3.05)

dev.off()

writexl::write_xlsx(list(`A`=g_pan_poscon$data,
                         `B`=g_pan_negcon$data,
                         C=g_pan_te_exemplary$data,
                         D=g_pan_te_expemplary1$data,
                         E=g_unr_1$data,
                         `F`=g_unr_2$data,
                         G=g_unr_3$data,
                         H=g_unr_4$data),
                    path = ifelse(exists("snakemake"),
                                  snakemake@output$xlsx,
                                  "~/Downloads/test.xlsx"))
