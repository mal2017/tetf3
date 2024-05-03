Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)
library(universalmotif)
library(memes)

# ------------------------------------------------------------------------------
# best meme motif(s) by match to known - the ones we want to highlight
# ------------------------------------------------------------------------------
comparison_fl <- "results/motifs/comparison/pan_denovo_comparison.meme.gg_df.rds"
comparison_fl <- snakemake@input$comparison
comparison <- read_rds(comparison_fl)

best_hits <- comparison |>
  filter(name == "pan" | str_detect(class,"Arch")) |>
  filter(discovery_eval < 3) |>
  group_by(denovo) |>
  slice_min(universalmotif_pval,with_ties = F) |>
  ungroup() |>
  filter(str_detect(name,"Archbold|^pan") | !str_detect(known,"known::MA"))

denovo_to_highlight <- best_hits |> pull(denovo) |>
  str_remove("denovo::") |>
  unique()
  
known_to_highlight <- best_hits |> 
  dplyr::select(denovo,known) |>
  mutate(across(everything(),~{str_remove(.x,".+::")})) |>
  deframe() 

# ------------------------------------------------------------------------------
# known motifs
# ------------------------------------------------------------------------------
known_pan_meme <- "results/motifs/known_motifs/all_known.meme"
known_pan_meme <- snakemake@input$known_meme
known_pan0 <- read_meme(known_pan_meme)
names(known_pan0) <- map_chr(known_pan0,`@`,name)
known_pan <- known_pan0[known_to_highlight]
known_pan <- unlist(known_pan)
# ------------------------------------------------------------------------------
# de novo motifs
# ------------------------------------------------------------------------------
# - meme motifs to highlight
meme_motifs_dir <- ifelse(exists("snakemake"), snakemake@input[["meme"]],
                         "results/motifs/meme_per_tf/pan")

meme <- paste0(meme_motifs_dir, "/meme.txt") |> memes::importMeme() |>  
  pull(motif) |> unlist()

names(meme) <- map_chr(meme,`@`,name)
meme <- meme[denovo_to_highlight]

# - all streme motifs
streme_motifs_dir <- ifelse(exists("snakemake"), snakemake@input[["streme"]],
                 "results/motifs/streme_per_tf/pan/")

streme <- paste0(streme_motifs_dir, "/streme.xml") |> memes::importStremeXML() |> pull(motif) |> unlist()


# - all homer motifs
homer_motifs_dir <- ifelse(exists("snakemake"), snakemake@input[["homer"]],
                            "results/motifs/homer_per_tf/pan/")
homer <- paste0(homer_motifs_dir, "/homerMotifs.all.motifs") |> 
  read_homer() |> unlist()

motifs <- c(meme,known_pan,streme,homer) |> unlist()
names(motifs) <- map_chr(motifs,`@`,name)
# ------------------------------------------------------------------------------
# comparison
# ------------------------------------------------------------------------------

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

p_df <- compare_motifs(motifs = motifs, 
                       compare.to = 1:length(denovo_to_highlight),
                       method = METHOD,nthreads = 4,
                       relative_entropy = RELATIVE_ENTROPY,
                       max.p=1, max.e = Inf, # all comparisons returned
                       normalise.scores = NORMALIZE_SCORES,
                       use.type = USETYPE,
                       score.strat = SCORE.STRAT) 

p_df <- p_df |> as_tibble() |> filter(!target %in%p_df$subject)

p_df <- p_df |>
  filter(Pval < 0.05) |>
  arrange(Pval)

p_df2 <- p_df |>
  mutate(motif_origin = case_when(str_detect(target,"Degenerate")~"Archbold (degenerate)",
                                  str_detect(target,'^\\d')~"HOMER",
                                  str_detect(target,"^m\\d+_")~"STREME",
                                  str_detect(target,"^MA0")~"jaspar",
                                  T~"other")) |>
  filter(motif_origin %in% c("HOMER","STREME")) |>
  group_by(motif_origin,subject) |>
  slice_min(Pval)

p_df3 <- p_df2 |>
  group_by(subject,motif_origin) |>
  slice_min(Pval) |>
  ungroup() |>
  arrange(subject,Pval) |>
  group_by(subject) |>
  mutate(rnkscore=mean(logPval)) |>
  filter(n() >= 2) |>
  ungroup() |>
  arrange(rnkscore) |>
  dplyr::select(-rnkscore)


plotting_grps <- p_df3 |>
  dplyr::select(subject,target) |>
  nest(data=-subject) |>
  mutate(data=map(data,~pull(.x,target))) |>
  deframe() |>
  imap(~{c(.y,.x)})


gg_l <- plotting_grps |>
  imap(~{
    view_motifs(c(motifs[known_to_highlight[.y]],motifs[.x]),
                method = METHOD, 
                score.strat = SCORE.STRAT, 
                text.size = 12,  
                normalise.scores = NORMALIZE_SCORES, 
                use.type = USETYPE)
  }
  ) 

saveRDS(p_df3, snakemake@output[["motif_comparison"]])
saveRDS(motifs, snakemake@output[["motifs_um"]])
saveRDS(gg_l, snakemake@output[["gg"]])

