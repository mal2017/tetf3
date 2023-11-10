library(tidyverse)
library(patchwork)
library(plotgardener)

lkup_path <- "results/resources/gene_symbol_lookup.tsv.gz"
res_path <- "results/deg/s2rplus.res.tsv.gz"

res0 <- read_tsv(res_path)
lkup <- read_tsv(lkup_path)

res <- res0 %>% left_join(lkup, by=c(feature = "gene_ID")) |>
  filter(!str_detect(feature,"FBgn")) |>
  arrange(logFC) |>
  mutate(rnk = dense_rank(logFC)) |>
  mutate(knockdown = if_else(comparison=="pan",comparison,"other"))

g_a <- res %>%
  group_by(comparison) |>
  summarize(logFC = mean(logFC)) |> 
  mutate(rnk = dense_rank(logFC)) |>
  ggplot(aes(rnk,logFC)) +
  geom_col(width = 1.01) +
  geom_col(data =  . %>% filter(comparison == "pan"), width=1.01, color="red") +
  geom_text(data = . %>% filter(comparison == "pan"), 
            aes(label=comparison),nudge_x = 10, nudge_y = -0.1,
            color="red", fontface="italic") +
  annotate(label=paste0("n=",length(unique(res$comparison))),x=Inf,y=-Inf, geom="text", hjust=1, vjust=0) +
  xlab("rank") + ylab("mean TE log fold change (limma)")


set.seed(2023)
resampled_res <- sample_n(res,size = 1e6, replace = T) |>
  mutate(knockdown="resampled") |>
  bind_rows(filter(res,comparison=="pan"))


ks_res <- ks.test(logFC ~ knockdown, data=resampled_res) |> 
  broom::tidy() |>
  mutate(label = paste0(method," p=",format.pval(p.value, digits=3)))

ks_res$label |> str_wrap(width = 30)

g_b <- resampled_res |>
  ggplot(aes(logFC,color=knockdown)) +
  stat_ecdf() +
  annotate(geom="text",label=str_wrap(ks_res$label,width = 30),x = max(resampled_res$logFC), y=0,hjust=1,vjust=0, size=rel(2)) +
  ylab("probability (CDF)") + xlab("log fold change (limma)") +
  scale_color_manual(values=c("pan"="red","resampled"="black"))



# ------------------------------------------------------------------------------
# create page

theme_set(theme_classic() + 
            theme(text = element_text(size=7)) +
            theme(plot.title = element_text(hjust = 0.5))
)

dir.create("results/figures/")

pdf("results/figures/s2rplus_supplement-01.pdf",width = 8.5, height = 11)

pageCreate(height = 11, showGuides=interactive())

plotGG(g_a, x = 0.5, y=0.5, width = 3.85,height = 2)
plotText("A", x = 0.5, y=0.5)


plotGG(g_b, x = 4.5, y=0.5, width = 3.7,height = 2)
plotText("B", x = 4.5, y=0.5)

dev.off()


