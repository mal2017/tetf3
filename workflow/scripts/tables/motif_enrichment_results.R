Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))

library(tidyverse)
library(writexl)
library(universalmotif)

# ------------------------------------------------------------------------------
# de novo results
# ------------------------------------------------------------------------------
# fun to import list of  and conver to df
import2 <- function(x) {
  
  program <- str_extract(x,"\\w+(?=_per_tf)")
  fac <- str_extract(x,"(?<=per_tf\\/).+?(?=\\/)")
  
  func <- list("meme"=universalmotif::read_meme,
       "homer"=universalmotif::read_homer,
       "streme"=universalmotif::read_meme)[[program]]

  m <- func(x)
  
  if (length(m)>0) {
    map_df(m, universalmotif::as.data.frame) |>
      dplyr::select(name,altname,consensus,eval,pval) |>
      mutate(program = program,TF=fac) |>
      dplyr::relocate(TF,program)
  }
}

# sheet 1-3 de novo
pages <- list(`MEME results`=Sys.glob("results/motifs/meme_per_tf/*/meme.txt") |> map_df(import2) |> dplyr::select(-pval),
           `STREME results`=Sys.glob("results/motifs/streme_per_tf/*/streme.txt") |> map_df(import2) |> dplyr::select(-pval),
           `HOMER results`=Sys.glob("results/motifs/homer_per_tf/*/homerMotifs.all.motifs") |> map_df(import2) |> dplyr::select(-eval)
  )


pages$`STREME results` <- Sys.glob("results/motifs/streme_per_tf_empirical_fdr/*_empirical_fdr.tsv") |>
  map_df(~{
    read_tsv(.x) |>
      dplyr::select(TF=te_group,name=motif_ID,empirical_fdr = fdr) |>
      drop_na()
  }) |>
  left_join(pages$`STREME results`,y=_, by=c("TF","name"))

# fimo meme hits
pages$`de novo pan motif hits on TE seqs` <- read_tsv("results/motifs/fimo_on_tes/denovo/pan/fimo.tsv")

# ------------------------------------------------------------------------------
# known matches to de novo results
# ------------------------------------------------------------------------------

sig_comparisons <- list(meme = "results/motifs/comparison/pan_denovo_comparison.meme.rds",
     homer="results/motifs/comparison/pan_denovo_comparison.homer.rds",
     streme = "results/motifs/comparison/pan_denovo_comparison.streme.rds") |>
  map_df(read_rds,.id="denovo_program") |>
  dplyr::select(-c(known_motif_obj,denovo_motif_obj,gg)) |>
  filter(universalmotif_padj < 0.1)

pages$`significant similarities to known motifs` <- sig_comparisons

# ------------------------------------------------------------------------------
# known motifs in csem peaks
# ------------------------------------------------------------------------------

# known motifs in csem peaks
pages$`known motifs in CSEM peaks` <- read_tsv("results/motifs/csem_peak_sea.known.pan.tsv.gz")

write_xlsx(pages, snakemake@output$xlsx)
