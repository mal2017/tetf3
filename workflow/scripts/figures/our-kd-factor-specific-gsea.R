Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)
library(patchwork)
library(plotgardener)

# ------------------------------------------------------------------------------
# barplot
# ------------------------------------------------------------------------------
g_rws <- read_rds("results/signatures/ourKD_gsea_randomwalks.gg_df.rds") |>
  filter(!signature_name %in% c("all_tes","TE.regulators","siRNA"))

g_rws <- g_rws |>
  dplyr::select(comparison,gg) |>
  deframe() |>
  map( ~{ .x & theme(axis.title = element_text(size=5), axis.text = element_text(size=5), plot.title = element_text(size=7, hjust=0.5))})

# ------------------------------------------------------------------------------
# factor-specific TE leading edge bar chart
# ------------------------------------------------------------------------------
g_a <- read_rds("results/signatures/ourKD_gsea_barplots.gg_list.rds")$factor.specific

g_b <- g_rws$knockdown2_Unr_female_head_Mef2.R_control_female_head_Mef2.R
g_c <- g_rws$knockdown2_vvl_female_head_Mef2.R_control_female_head_Mef2.R

# ------------------------------------------------------------------------------
# scrna tfoi boxplot
# ------------------------------------------------------------------------------
tfs <- read_tsv("resources/Drosophila_melanogaster_TF.txt")

tes <- jsonlite::read_json("upstream/te_element_lookup.json") %>%
  names()

# read in TE TF correlations we identifed from the supercell results
feature_correlations <- read_rds("results/calderon22/fca_reanalysis_correlations.rds") |>
  dplyr::select(lineage,res=res.spqn) |>
  unnest(res)

# asking if pan is highly correlated with TEs in general
g_d <- 
  feature_correlations %>%
  filter(feature %in% tfs$Symbol & y %in% tes) |>
  dplyr::select(y,feature,coef,p,padj) |>
  distinct() |>
  mutate(feature2 = if_else(feature%in%c("pan","Unr","vvl","NfI","CG16779"),feature,"other")) |>
  mutate(feature2 = fct_reorder(feature2,coef)) |>
  mutate(feature2 = fct_relevel(feature2,"other")) |>
  ggplot(aes(feature2,coef)) +    
  geom_boxplot(outlier.size = 0.2) +
  ggpubr::stat_compare_means(ref="other",size=2,label = "p.format") +
  xlab("") + ylab("weighted correlation")



# ------------------------------------------------------------------------------
# create page

theme_set(theme_classic() + 
            theme(text = element_text(size=7)) +
            theme(plot.title = element_text(hjust = 0.5))
)

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())
figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

plotGG(g_a, x=0.5, y=0.5, width=3.75, height=2.5)
plotText("A", x = 0.5, y=0.5)

plotGG(g_b, x = 4.25, y=0.7, width = 3.25,height = 2.5)
plotText("B", x = 4.5, y=0.5)


plotGG(g_c, x = .5, y=3.5, width = 3.25,height = 2.5)
plotText("C", x = .5, y=3.5)

plotGG(g_d, x = 4.25, y=3.5, width = 3.25,height = 2.5)
plotText("D", x = 4.25, y=3.5)

dev.off()


writexl::write_xlsx(list(A=g_a$data,
                         B=g_b$data,
                         C=g_c$data),
                    path = ifelse(exists("snakemake"),
                                  snakemake@output$xlsx,
                                  "~/Downloads/test.xlsx"))