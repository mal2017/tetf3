library(tidyverse)
library(Category)

zads <- read_tsv("results/resources/zad_genes.tsv")
lkup0 <- read_tsv("resources/Drosophila_melanogaster_TF.txt")
fbgn_lkup <- read_tsv("https://ftp.flybase.net/releases/FB2023_06/precomputed_files/genes/fbgn_annotation_ID_fb_2023_06.tsv.gz",skip = 3)

# ------------------------------------------------------------------------------
# generate lookup list for gene symbols -> ids
# ------------------------------------------------------------------------------
# some point upstream I made the names syntactic, so I reconstruct these so I can look up with this
lkup <- fbgn_lkup |> 
  mutate(syntactic_symbol = str_replace_all(`##gene_symbol`,"[\\(\\)-]",".")) |> 
  dplyr::select(syntactic_symbol, Ensembl=`primary_FBgn#`) |>
  distinct() |>
  deframe()

# some IDs are included as CG#s
annotID_lkup <- dplyr::select(fbgn_lkup, syntactic_symbol=annotation_ID ,  Ensembl=`primary_FBgn#`) |> deframe()

# if any of these can rescue unknowns, add them here
for (g in names(annotID_lkup)) {
  if (!g %in% names(lkup)) {
    lkup[[g]] <- annotID_lkup[[g]]
  }
}

# weird cases caused by AnimalTFDB not using current gene symbles and/or just using FBgn instead
edge_cases <- list(h="FBgn0263118", FBgn0263118='FBgn0288888')
lkup <- c(lkup,edge_cases)

# fun convers phylosignal coefs to a named list
clean_and_rename <- \(x) pull(x,coef) |> str_extract("(?<=score_).+") |> unique() |> set_names() |> map(~lkup[[.x]])

# ------------------------------------------------------------------------------
# import hits and universe
# ------------------------------------------------------------------------------

# hit renaming for use with downstream tools
hits <- read_tsv("results/phylosignal/phylosignal_filtered_hits.tsv.gz") |> clean_and_rename() 
universe <- read_rds("results/phylosignal/phylosignal_df.rds") |> clean_and_rename()
universe <- universe[!is.na(names(universe))] # NAs in names because the pos and negative controls for phylosignal have no gene name -  get rid of this here

stopifnot(!any(is.na(hits)))
stopifnot(!any(is.na(universe)))

# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# nothing w.r.t zads
clusterProfiler::enricher(gene = as.character(unlist(hits)),
         universe = as.character(unlist(universe)),
         TERM2GENE = relocate(mutate(tibble(gene_id=zads$ensembl_gene_id[zads$ensembl_gene_id %in% unlist(universe)]),term="ZAD-ZNF"),term),pvalueCutoff = 1)

# not used
convert <- function(x,type="ENTREZID") {
  bitr(x,fromType = "FLYBASE",toType = type,OrgDb =  org.Dm.eg.db::org.Dm.eg.db)
}

# nothing in go
clusterProfiler::enrichGO(unlist(hits),OrgDb = org.Dm.eg.db::org.Dm.eg.db, keyType = "FLYBASE",universe = unlist(universe),ont = "ALL",pvalueCutoff = 1)



peps <- Biostrings::readAAStringSet("results/peptides/longest.fasta")
names(peps) <- names(peps) |> str_extract(".+(?=\\.FBpp)")

poi <- peps[names(peps) %in% unlist(hits)]
