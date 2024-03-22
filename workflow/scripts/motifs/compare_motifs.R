library(tidyverse)
library(universalmotif)
library(memes)

# ------------------------------------------------------------------------------
# known motifs
# ------------------------------------------------------------------------------
known_motifs_meme <- "results/motifs/known_motifs/all_known.meme"
known_motifs_meme <- snakemake@input$known_meme
known_motifs <- read_meme(known_motifs_meme)

names(known_motifs) <- known_motifs %>% map_chr( `@`, name)
names(known_motifs) <- paste0("known::", names(known_motifs))
known_motifs <- imap(known_motifs, ~{.x@name <- .y; .x})
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
  
  # allow graceful exit if there's an empty meme file (NfI does this)
  if (file.info(paste0(motifs_dir,"/meme.txt"))$size == 0) {
    saveRDS(tibble(), snakemake@output[["motif_comparison"]])
    saveRDS(matrix(), snakemake@output[["motif_similarity"]])
    quit()
  }
  motifs <- paste0(motifs_dir, "/meme.txt") |> memes::importMeme() |> pull(motif)
  names(motifs) <- motifs %>% map_chr( `@`, name)
} else if (snakemake@params$motif_program == "homer") {
  motifs <-paste0(motifs_dir, "/homerMotifs.all.motifs")
  motifs <- read_homer(motifs)
  names(motifs) <- map_chr(motifs,`@`,consensus)
}

names(motifs) <- paste0("denovo::", names(motifs))
motifs <- imap(motifs, ~{.x@name <- .y; .x})
# ------------------------------------------------------------------------------
# comparison
# ------------------------------------------------------------------------------
all_motifs <- c(known_motifs,motifs)

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

# this is a little more computation than necessary, but
# we compare all to all to make it easy to get each denovo motif as a target
p_df0 <- compare_motifs(motifs = all_motifs, 
                       compare.to = 1:length(all_motifs),
                       method = METHOD,nthreads = 2,
                       relative_entropy = RELATIVE_ENTROPY,
                       max.p=1, max.e = Inf, # all comparisons returned
                       normalise.scores = NORMALIZE_SCORES,
                       use.type = USETYPE,
                       score.strat = SCORE.STRAT)

# rename and select cols to make clear what is denovo and which motif each denovo
# is compared to
p_df <- p_df0 |>
  as_tibble() |>
  filter(subject!=target) |>
  filter(str_detect(target,"denovo")) |>
  filter(!str_detect(subject,"denovo")) |>
  dplyr::rename(denovo="target",known="subject") |>
  dplyr::select(denovo,known,universalmotif_pval=Pval)

# annotate with the Evalue of the discovery program
denovo_discovery_signif <- motifs |> 
  enframe(name = "denovo",value = "motif_obj") |>
  mutate(discovery_eval=map_dbl(motif_obj, ~{ifelse(is_empty(.x@eval),NA,.x@eval)})) |>
  mutate(discovery_pval=map_dbl(motif_obj, ~{ifelse(is_empty(.x@pval),NA,.x@pval)})) |>
  dplyr::select(-motif_obj)

p_df <- left_join(p_df, denovo_discovery_signif,by="denovo")

# get target factor names
p_df <- known_motifs |>
  map_chr(`@`,"altname") |>
  enframe(name = "known",value = "name") |>
  left_join(p_df, y=_, by="known") |>
  dplyr::relocate(denovo,known)

# get adjusted p, grouping on tests related to the same denovo motif
p_df <- p_df |>
  group_by(denovo) |>
  mutate(universalmotif_padj = p.adjust(universalmotif_pval, method="BH")) |>
  ungroup()

# get motif alignments using identical params as used for hypothesis testing above
p_df <- mutate(p_df, 
               known_motif_obj = map(known,~{all_motifs[[.x]]}),
               denovo_motif_obj = map(denovo,~{all_motifs[[.x]]})) %>%
  mutate(gg = pmap(list(known_motif_obj,denovo_motif_obj,universalmotif_padj),.f =  function(x,y,z) {
    if (z >= 0.1) {
      return(NA)
    } else {
      return(view_motifs(c(x,y),
                  method = METHOD, 
                  score.strat = SCORE.STRAT, 
                  sort.positions = T,
                  text.size = 7,  
                  normalise.scores = NORMALIZE_SCORES, use.type = USETYPE))
    }
    
    }))



saveRDS(p_df, snakemake@output[["motif_comparison"]])

