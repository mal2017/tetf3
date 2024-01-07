library(tidyverse)
library(patchwork)
library(plotgardener)

lkup_path <- "results/resources/gene_symbol_lookup.tsv.gz"
res_path <- "results/deg/s2rplus.res.tsv.gz"

res0 <- read_tsv(res_path)
lkup <- read_tsv(lkup_path)

# ------------------------------------------------------------------------------
# barplot
# ------------------------------------------------------------------------------

res <- res0 %>% left_join(lkup, by=c(feature = "gene_ID")) |>
  filter(!str_detect(feature,"FBgn")) |>
  arrange(logFC) |>
  mutate(rnk = dense_rank(logFC)) |>
  mutate(knockdown = if_else(comparison %in% c("pan","NfI","CG1679","vvl","Unr"),comparison,"other"))

g_b <- res %>%
  group_by(comparison, knockdown) |>
  summarize(logFC = mean(logFC), .groups = "drop") |> 
  mutate(rnk = dense_rank(logFC)) |>
  ggplot(aes(rnk,logFC)) +
  geom_col(width = 1.01) +
  geom_col(data =  . %>% filter(knockdown != "other"), width=1.01, color="red") +
  geom_text(data = . %>% filter(knockdown != "other"), 
            aes(label=comparison),nudge_x = 10, nudge_y = -0.1,
            color="red", fontface="italic") +
  annotate(label=paste0("n=",length(unique(res$comparison))),x=Inf,y=-Inf, geom="text", hjust=1, vjust=0) +
  xlab("rank") + ylab("mean TE log fold change (limma)")

# ------------------------------------------------------------------------------
# pan FC distrib KS test
# ------------------------------------------------------------------------------

set.seed(2023)
resampled_res <- sample_n(res,size = 1e6, replace = T) |>
  mutate(knockdown="resampled") |>
  bind_rows(filter(res,comparison=="pan"))


ks_res <- ks.test(logFC ~ knockdown, data=resampled_res) |> 
  broom::tidy() |>
  mutate(label = paste0(method," p=",format.pval(p.value, digits=3)))

ks_res$label |> str_wrap(width = 30)

g_c <- resampled_res |>
  ggplot(aes(logFC,color=knockdown)) +
  stat_ecdf() +
  annotate(geom="text",label=str_wrap(ks_res$label,width = 30),x = max(resampled_res$logFC), y=0,hjust=1,vjust=0, size=rel(2)) +
  ylab("probability (CDF)") + xlab("log fold change (limma)") +
  scale_color_manual(values=c("pan"="red","resampled"="black"))

# ------------------------------------------------------------------------------
# All TE leading edge bar chart
# ------------------------------------------------------------------------------
x <- read_rds("results/signatures/ourKD_gsea.rds")

x <- x |> mutate(lab = str_replace_all(str_extract(comparison,"(?<=knockdown2_).+(?=_control)"), "_", " / "))

plot_barchart <- function(dat, cmp) {
  dat |>
    filter(signature_name == cmp) |>
    mutate(lab = fct_reorder(lab, pvalue)) |>
    ggplot(aes(-log10(pvalue), lab)) +
    geom_col() +
    geom_vline(xintercept = -log10(0.05),color="red",linetype="dashed") +
    ylab("RNAi / sex / sample / driver")
}

g_a <- plot_barchart(x, "all_tes")



# ------------------------------------------------------------------------------
# create page

theme_set(theme_classic() + 
            theme(text = element_text(size=7)) +
            theme(plot.title = element_text(hjust = 0.5))
)

dir.create("results/figures2/")

pdf(snakemake@output$pdf,width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_a, x=0.5, y=0.5, width=3.75, height=2.5)
plotText("A", x = 0.5, y=0.5)


plotGG(g_b, x = 4.25, y=0.7, width = 3.85,height = 2)
plotText("B", x = 4.5, y=0.5)


plotGG(g_c, x = .5, y=3.5, width = 3.7,height = 2)
plotText("C", x = .5, y=3.5)

dev.off()


