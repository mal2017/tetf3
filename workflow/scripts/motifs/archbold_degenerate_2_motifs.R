library(tidyverse)
library(universalmotif)
library(Biostrings)
library(tesseract)
library(cluster)

# ------------------------------------------------------------------------------
# - archbold motifs
# ------------------------------------------------------------------------------
eng <- tesseract("eng")
text <- tesseract::ocr("https://doi.org/10.1371/journal.pgen.1004591.s007", engine = eng)

hmg_motifs_df <- read_delim(text, skip = 8)

hmg_motifs_df <- hmg_motifs_df  %>%
  unite(dm3_location, c("Chr.","Location"),sep = ":") |>
  dplyr::select(dm3_location, Sequence) |>
  mutate(nm = paste0("Archbold14-degenerate::",Sequence))

seqs <- hmg_motifs_df |>
  dplyr::select(nm, Sequence) |>
  deframe() |>
  DNAStringSet()


# -----------------------------------------------------------------------------------
# get reasonable k for clustering degenerate sequences
# ------------------------------------------------------------------------------------
# idea from https://www.bioconductor.org/packages/devel/bioc/vignettes/universalmotif/inst/doc/SequenceSearches.pdf
# but reworked and with kmeans run at the end
# also see https://bioconductor.org/books/3.17/OSCA.advanced/clustering-redux.html
# for info on how why I reassign "closest" below
get_clustering <- function(s, klet_size=3, klusters=2) {
  pmat <- get_bkg(s,k=klet_size, merge.res = F) |>
    as_tibble() |>
    dplyr::select(sequence, klet, "probability") |>
    pivot_wider(names_from = "klet", values_from = "probability") |>
    column_to_rownames("sequence")
  
  dmat <- dist(pmat)
  cl <- kmeans(pmat, centers = klusters)
  sil <- cluster::silhouette(cl$cluster,dist=dmat) |>
    as_tibble() |>
    mutate(across(c(cluster,neighbor),as.factor)) |>
    mutate(closest = as.factor(ifelse(sil_width > 0,cluster, neighbor)))
  
  tibble(cl = list(cl), 
         d=list(dmat), 
         p=list(pmat),
         sil = list(as_tibble(sil)))
}


# makes a tibble that contains klet freq mat, dist mat, kmeens objs + qc for
# a number of combinations of klet size and centers (kmeans k param).
# also some diagnostic plots
set.seed(2023)
sils <- expand_grid(klet_size=2:4, centers = 2:6) |>
  mutate(clustering = map2(klet_size, centers, ~{get_clustering(seqs,.x,.y)})) |>
  unnest(clustering) |> 
  mutate(mean_width = map_dbl(sil,~{mean(.x$sil_width)})) |>
  mutate(gg = map(sil,.f=~{
    ggplot(.x, aes(cluster,sil_width, color=closest)) +
      geom_jitter(width=0.2) +
      coord_cartesian(ylim=c(-1,1)) +
      xlab("cluster") +
      ylab("silhouette width")
  })) |>
  mutate(ggt = pmap(list(klet_size,centers,mean_width), 
                    .f=function(x,y,z) {sprintf("klet size:%s / centers: %s / mean sil. width: %s",x,y,format(z,digits = 3))})) |>
           mutate(gg = map2(gg,ggt, ~{.x + labs(title = .y)})) |>
  dplyr::select(-ggt) |> 
  mutate(n_mislabeled = map_int(sil,~{sum(as.character(.x$cluster)!=as.character(.x$closest))}))

# identify best clustering, prioritizing no 
# candidate mislabeled and then picking by mean sil width
best_sil <- filter(sils, n_mislabeled == 0) |>
  slice_max(mean_width, n=1, with_ties = T) |>
  slice_max(centers, n=1, with_ties = F)

motifs <- split(names(best_sil$cl[[1]]$cluster), best_sil$cl[[1]]$cluster) |>
  map(~{seqs[.x]}) |>
  map(create_motif) |>
  map(trim_motifs, min.ic=1)

names(motifs) <- names(motifs) |> str_glue("Archbold.Degenerate.{ix}",ix=_)

motifs <- imap(motifs, ~{.x@name <- .y; .x@altname <- .y; .x})

write_meme(motifs, snakemake@output$meme)
write_rds(sils, snakemake@output$sils)
Biostrings::writeXStringSet(seqs,snakemake@output$seqs)


#library(patchwork)
#map(motifs, universalmotif::view_motifs) |> Reduce(`+`,x=_)



