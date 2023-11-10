
library(tidyverse)
library(universalmotif)
library(memes)
library(tesseract)


# ------------------------------------------------------------------------------
# known motifs
# ------------------------------------------------------------------------------
known_pan_meme <- "results/motifs/known_motifs/all_known.meme"
known_pan_meme <- snakemake@input$known_meme
known_pan <- read_meme(known_pan_meme)

names(known_pan) <- known_pan %>% map_chr( `@`, name)
names(known_pan) <- paste0("known::", names(known_pan))
known_pan <- imap(known_pan, ~{.x@name <- .y; .x})
# ------------------------------------------------------------------------------
# de novo motifs
# ------------------------------------------------------------------------------
motifs_dir <- ifelse(exists("snakemake"), snakemake@input[["denovo"]],
                 "results/motifs/streme_per_tf/pan/")

# get the motifs worth comparing to known (eval/pval < 0.05 (lower for homer, bc it reports more seemngly significant motifs))
if (snakemake@params$motif_program == "streme") {
  motifs <- paste0(motifs_dir, "/streme.xml") |> memes::importStremeXML() |> pull(motif)
  names(motifs) <- motifs %>% map_chr( `@`, name)
} else if (snakemake@params$motif_program == "meme") {
  motifs <- paste0(motifs_dir, "/meme.txt") |> memes::importMeme() |> filter(eval<0.05) |> pull(motif)
  names(motifs) <- motifs %>% map_chr( `@`, name)
} else if (snakemake@params$motif_program == "homer") {
  motifs <-paste0(motifs_dir, "/homerMotifs.all.motifs")
  motifs <- read_homer(motifs)
  names(motifs) <- map_chr(motifs,`@`,consensus)
  motifs <- motifs[map_lgl(motifs, ~{.x@pval <= 0.001})]
}

names(motifs) <- paste0("denovo::", names(motifs))
motifs <- imap(motifs, ~{.x@name <- .y; .x})
# ------------------------------------------------------------------------------
# comparison
# ------------------------------------------------------------------------------
all_motifs <- c(motifs, known_pan)

# default args, except for score.strat - the sampling distribution
# of pearson's rho is skewed, so averaging after FZT makes more sense
# I used the same approach to average the rho values for gene x gene correlation
# among salmon replicates
# and normalize scores, which favors more complete alignments

# settings
METHOD = "PCC"
SCORE.STRAT= "fzt"
RELATIVE_ENTROPY = F
NORMALIZE_SCORES=T
USETYPE="PPM"

# get pval df

p_df <- compare_motifs(motifs = all_motifs, 
                       compare.to = 1:length(all_motifs),
                       method = METHOD,nthreads = 4,
                       relative_entropy = RELATIVE_ENTROPY,
                       max.p=1, max.e = Inf, # all comparisons returned
                       normalise.scores = NORMALIZE_SCORES,
                       use.type = USETYPE,
                       score.strat = SCORE.STRAT)

p_df <- p_df |>
  as_tibble() |>
  filter(subject!=target) |>
  filter(str_detect(target,"denovo") | str_detect(subject,"denovo")) |>
  filter(!(str_detect(target,"denovo") & str_detect(subject,"denovo"))) |>
  mutate(denovo=map2_chr(target,subject, ~if_else(str_detect(.x,"denovo"),.x,.y))) |>
  mutate(known=map2_chr(target,subject, ~if_else(str_detect(.x,"denovo"),.y,.x))) |>
  mutate(padj = p.adjust(Pval, method="BH")) |>
  arrange(padj) |>
  mutate(motifs = map2(known,denovo, ~{c(all_motifs[[.x]],all_motifs[[.y]])})) %>%
  mutate(gg = map(motifs, ~view_motifs(.x,method = METHOD, score.strat = SCORE.STRAT, text.size = 7,  normalise.scores = NORMALIZE_SCORES, use.type = USETYPE)))

# get target factor names
p_df <- known_pan |>
  map_chr(`@`,"altname") |>
  enframe(name = "known",value = "name") |>
  right_join(p_df, by="known") |>
  dplyr::relocate(denovo,known)

# get sim/dist mat
mat <- compare_motifs(motifs = all_motifs,
                      method = METHOD,
                      normalise.scores = NORMALIZE_SCORES,
                      relative_entropy = RELATIVE_ENTROPY,
                      use.type = USETYPE,
                      score.strat = SCORE.STRAT)

saveRDS(p_df, snakemake@output[["motif_comparison"]])
saveRDS(mat, snakemake@output[["motif_similarity"]])

