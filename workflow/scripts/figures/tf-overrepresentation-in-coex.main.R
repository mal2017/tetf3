library(patchwork)
library(plotgardener)
library(tidyverse)
library(clusterProfiler)
library(ggpubr)
library(ggdensity)


# ------------------------------------------------------------------------------
# generate gsea random walk plots describing enriched gene groups

male_gsea <- read_rds("results/enrichment/sig_main_male_max_abs_estimate_qnorm.gg_gsea.rds")
female_gsea <- read_rds("results/enrichment/sig_main_female_max_abs_estimate_qnorm.gg_gsea.rds")

make_dotplot <- function(.x) {
  as_tibble(.x) |>
    filter(p.adjust < 0.1) |>
    mutate(group = if_else(str_detect(ID,"AnimalTFDB|ZAD_ZNF"),ID,"Flybase gene group")) |>
    mutate(lab=fct_reorder(str_wrap(str_trunc(ID,50),20, whitespace_only = F), -abs(pvalue))) |>
    mutate(`leading edge members`=str_count(core_enrichment,'/')) |>
    ggplot(aes(lab, -log10(pvalue),fill=group)) +
    geom_col() +
    scale_fill_grey(start = 0.3,end = 0.6) +
    ylab("log10(p)") + xlab("") +
    theme(legend.spacing = unit(0.01,"in")) +
    labs(size="leading edge genes") +
    theme(legend.key.size = unit(.05,"inches"), axis.text = element_text(size=unit(5,"pt"))) +
    geom_text(aes(label=sprintf("n=%s",`leading edge members`), y= 0.1*max(-log10(pvalue))), size=rel(2),color="white",hjust=0.5) +
    theme(axis.text.x = element_text(angle=45, hjust=1))
}

g_a <- make_dotplot(male_gsea)

g_b <- make_dotplot(female_gsea)



g <- g_a + g_b + plot_layout(guides="collect",widths = c((6/4)*3,3))

# ------------------------------------------------------------------------------
# make page
# ------------------------------------------------------------------------------


theme_set(theme_classic() + 
            theme(text = element_text(size=unit(7,"pt"))) +
            theme(plot.title = element_text(hjust = 0.5))
)

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())
figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

figtitle = ifelse(exists("snakemake"),snakemake@params$figtitle,"Figure X")
plotText(figtitle,x=0,y=0,just = c("left","top"))

plotGG(g, x = 1, y=0.5, width = 7,height = 2.75)
plotText("A", x = 1.25, y=0.5)
plotText("B", x = 4.25, y=0.5)

dev.off()

writexl::write_xlsx(list(A=g_a$data,
                         B=g_b$data),
                    path = ifelse(exists("snakemake"),
                                  snakemake@output$xlsx,
                                  "~/Downloads/test.xlsx"))
