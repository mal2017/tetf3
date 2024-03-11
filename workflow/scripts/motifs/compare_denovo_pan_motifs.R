library(tidyverse)
library(universalmotif)
library(memes)

# ------------------------------------------------------------------------------
# known motifs
# ------------------------------------------------------------------------------
known_pan_meme <- "results/motifs/known_motifs/all_known.meme"
known_pan_meme <- snakemake@input$known_meme
known_pan <- read_meme(known_pan_meme)

names(known_pan) <- known_pan %>% map_chr( `@`, altname)
known_pan <- known_pan[names(known_pan)=="pan"][[2]]
known_pan@name <- known_pan@name |> paste0("JASPAR::",y=_)
# ------------------------------------------------------------------------------
# de novo motifs
# ------------------------------------------------------------------------------
streme_motifs_dir <- ifelse(exists("snakemake"), snakemake@input[["streme"]],
                 "results/motifs/streme_per_tf/pan/")

streme <- paste0(streme_motifs_dir, "/streme.xml") |> memes::importStremeXML() |> pull(motif)

streme <- streme$m4_AAAAATGGCRCATRG
streme@name <- streme@name |> paste0("STREME::",y=_)

meme_motifs_dir <- ifelse(exists("snakemake"), snakemake@input[["meme"]],
                            "results/motifs/meme_per_tf/pan")

meme <- paste0(meme_motifs_dir, "/meme.txt") |> memes::importMeme() |> filter(eval<0.05) |> pull(motif)
meme <- meme[[1]]
meme@name <- meme@name |> paste0("MEME::",y=_)

homer_motifs_dir <- ifelse(exists("snakemake"), snakemake@input[["homer"]],
                            "results/motifs/homer_per_tf/pan/")
homer <- paste0(homer_motifs_dir, "/homerMotifs.all.motifs") |> 
  read_homer()

names(homer) <- map_chr(homer,`@`,consensus)
homer <- homer[map_lgl(homer,~{.x@pval <= 0.001})]
homer <- homer[["GCTAKTTWGMTG"]]
homer@name <- homer@name |> paste0("HOMER::",y=_)

motifs <- c(known_pan,meme,streme,homer)
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
                       compare.to = 1,
                       method = METHOD,nthreads = 4,
                       relative_entropy = RELATIVE_ENTROPY,
                       max.p=1, max.e = Inf, # all comparisons returned
                       normalise.scores = NORMALIZE_SCORES,
                       use.type = USETYPE,
                       score.strat = SCORE.STRAT) 

p_df <- p_df |>
  as_tibble() |>
  mutate(padj = p.adjust(Pval, method="BH")) |>
  arrange(padj)

gg <- view_motifs(motifs,
                  method = METHOD, 
                  score.strat = SCORE.STRAT, 
                  text.size = 12,  
                  normalise.scores = NORMALIZE_SCORES, 
                  use.type = USETYPE)

saveRDS(p_df, snakemake@output[["motif_comparison"]])
saveRDS(motifs, snakemake@output[["motifs_um"]])
saveRDS(gg, snakemake@output[["gg"]])

