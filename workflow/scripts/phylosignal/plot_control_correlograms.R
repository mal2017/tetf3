library(phylosignal)

ps_fl <- "results/phylosignal/phylosignal.rds"
ps_fl <- snakemake@input$phylosignal
x <- readRDS(ps_fl)

bm <- phyloCorrelogram(x$p4d, "bm")
ou <- phyloCorrelogram(x$p4d, "ou")
random <- phyloCorrelogram(x$p4d, "random")

res <- list(bm = bm, random = random, ou=ou)

saveRDS(res, snakemake@output$rds)