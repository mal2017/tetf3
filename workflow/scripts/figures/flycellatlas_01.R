Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


# read in results of supercell analysis
supercell <- readRDS("results/calderon22/fca_reanalysis_supercell.rds")

GE<-supercell$GE
SC<-supercell$SC

rm(supercell);gc()

library(plotgardener)
library(tidyverse)
library(DiagrammeRsvg)
library(rsvg)
library(DiagrammeR)
library(SuperCell)

tfs <- read_tsv("resources/Drosophila_melanogaster_TF.txt")

tes <- jsonlite::read_json("upstream/te_element_lookup.json") %>%
  names()

# read in TE TF correlations we identifed from the supercell results
feature_correlations <- read_rds("results/calderon22/fca_reanalysis_correlations.rds") |>
  dplyr::select(lineage,res=res.spqn) |>
  unnest(res)

# ------------------------------------------------------------------------------
# plot supercell sizes to help explain that analysis
# ------------------------------------------------------------------------------

g_supercell_size <- SC$supercell_size |>
  enframe(value = "n_cells") |>
  ggplot(aes(n_cells)) +
  geom_histogram() +
  xlab("barcodes per supercell") + ylab("N") +
  scale_x_log10()

# ------------------------------------------------------------------------------
# plot results of supercell correlation analysis
# ------------------------------------------------------------------------------

# asking if pan is highly correlated with TEs in general
g_pan_highly_corr_with_tes <- 
feature_correlations %>%
  filter(feature %in% tfs$Symbol & y %in% tes) |>
  dplyr::select(y,feature,coef,p,padj) |>
    distinct() |>
  mutate(feature2 = if_else(feature%in%c("pan","Unr","vvl","NfI","CG16779"),feature,"other")) |>
  mutate(feature2 = fct_reorder(feature2,coef)) |>
  mutate(feature2 = fct_relevel(feature2,"other")) |>
  ggplot(aes(feature2,coef)) +    
  geom_boxplot() +
  ggpubr::stat_compare_means(ref="other",size=1) +
  xlab("") + ylab("weighted correlation")


g_a <- grViz("
digraph boxes_and_circles {

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]

  # several 'node' statements
  node [shape = box,
        fontname = Arial]
  data [label='adult head scRNA-seq data\n(FlyCellAtlas)'];
  merge [label='quantify TE/gene expression with alevin-fry'];
  metacell [label='create metacells with SuperCell'];
  corr [label='calculate weighted TE/gene correlations'];

  data -> merge;
  merge -> metacell;
  metacell -> corr;
}
")

cartoon_temp <- tempfile()
#grVizToPNG(g_a, filename = "project_overview.png")
export_svg(g_a) |>
  charToRaw() |>
  rsvg_svg(cartoon_temp,width = 2000,height = 2000)

g_a_cartoon <- magick::image_read_svg(cartoon_temp) |> magick::image_ggplot(interpolate = T)


# plotting page 1 --------------------------------------------------------------------
theme_set(theme_classic() + theme(text=element_text(size=unit(7,"pt"))))

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 8.5, height = 11)

pageCreate(width =8.5, height = 11, default.units = "inches", showGuides = interactive())
figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

pa <- plotGG(plot = g_a_cartoon, x = 1, y=0.3, width = 3, height=2.25)
plotText(label = "A", x = 1, y = 0.5)

pb <- plotGG(plot = g_supercell_size, x = 4.5, y=0.5, width = 3, height=2.25)
plotText(label = "B", x = 4.5, y = 0.5)

pc <- plotGG(plot = g_pan_highly_corr_with_tes, x = 1,  y=3, width = 3.25, height=2.25)
plotText(label = "C", x = 1, y = 3)


dev.off()

writexl::write_xlsx(list(B=g_supercell_size$data,
                         C=g_pan_highly_corr_with_tes$data),
                    path = ifelse(exists("snakemake"),
                                  snakemake@output$xlsx,
                                  "~/Downloads/test.xlsx"))
