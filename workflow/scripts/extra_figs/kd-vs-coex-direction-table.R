library(tidyverse)
library(DescTools)
library(plotgardener)

kd <- read_rds("results/signatures/ourKD_gsea.rds") |>
  mutate(lab = str_replace_all(str_extract(comparison,"(?<=knockdown2_).+(?=_control)"), "_", " / ")) |>
  mutate(signature_type = if_else(signature_name %in% c("all_tes","TE.regulators"),signature_name,"Factor.specific")) |>
  dplyr::select(lab, kd, signature_type, padj,direction=NES) |>
  mutate(direction=if_else(padj < 0.1,sign(direction),0)) |>
  dplyr::select(-padj) |>
  mutate(signature_type = case_match(signature_type,"Factor.specific"~"Factor-specific TEs\nexpression change",
                                     "TE.regulators"~"TE silencers\nexpression change",
                                     "all_tes"~"All TEs\nexpression change")) |>
  pivot_wider(names_from = signature_type, values_from = "direction")

lms <- read_tsv("upstream/final-models.collected-info.tsv.gz") |>
  filter(gene_symbol %in% kd$kd)

toplot <- lms |>
  filter(significant_x) |>
  group_by(gene_symbol) |>
  summarize(`TE coexpression` = Mode(sign(estimate.qnorm)),.groups = "drop") |>
  left_join(kd, by=c(gene_symbol ="kd")) |>
  dplyr::relocate(lab) |>
  pivot_longer(-c(lab,gene_symbol)) |>
  mutate(name=fct_relevel(as_factor(name),"TE coexpression")) |>
  mutate(name=fct_relevel(name,"TE silencers\nexpression change",after=1)) |>
  mutate(value=case_match(value,1~"positive",
                          -1~"negative",
                          0~"n.s."))

gg <- toplot |>
ggplot(aes(name,lab)) +
    geom_point(aes(shape=value,fill=value),size=rel(7)) +
  scale_shape_manual(values = c("positive"=24,"negative"=25,"n.s."=1)) +
  scale_fill_manual(values = c("positive"="darkgreen","negative"="red","n.s."="darkgray")) +
  xlab("") +
  ylab("knockdown experiment") +
  theme(legend.title = element_blank())



# ------------------------------------------------------------------------------
# create page
# ------------------------------------------------------------------------------

theme_set(theme_classic() + 
            theme(text = element_text(size=7))
)

dir.create("results/figures2/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(gg, x = 1, y=0.5, width = 6.5,height = 4)

dev.off()
