library(tidyverse)
library(DESeq2)
library(Biostrings)
library(sva)
library(plyranges)

# ------------------------------------------------------------------------------
# setup
# ------------------------------------------------------------------------------

se_path <- "upstream/kd.se.gene.0.rds"
se_path <- snakemake@input[["se"]]

gene_universe <- read_tsv("resources/fbgn_fbtr_fbpp_expanded_fb_2021_04.tsv.gz", skip = 4)

allowed_genes <- gene_universe %>% filter(gene_type %in% c("protein_coding_gene","ncRNA_gene"))

se <- read_rds(se_path)

colData(se) <- colData(se) |>
  as_tibble(rownames = "name")|>
  mutate(across(where(is_character),.f=~str_replace(.x,"NFI","NfI"))) |>
  mutate(batch = str_extract(path_r1,"novo\\d*")) |>
  mutate(knockdown2 = paste(knockdown,tissue,driver,sep="_")) |>
  column_to_rownames("name") |>
  DataFrame()

se <- se [rownames(se) %in% allowed_genes$gene_ID | 
               !str_detect(rownames(se),"FBgn"),]


# ------------------------------------------------------------------------------
# pre filtering 
# ------------------------------------------------------------------------------
keep <- rowSums(assay(se) >= 10) >= 5 &
  rowSums(assay(se) >= 1) >= 2

se <- se[keep,]

# ------------------------------------------------------------------------------
# separate heads and the gonads to minimize the amount of batch correction
# necessary
# ------------------------------------------------------------------------------
se.head <- se[,se$driver %in% c("Mef2.R")]

se.gonad <- se[, (se$knockdown %in% c("CG16779","control") & se$driver %in% c("aTub")) | 
           (se$knockdown %in% c("control","pan") & se$driver %in% c("aTub","tj"))]

# ------------------------------------------------------------------------------
# batch correction for heads
# ------------------------------------------------------------------------------
mat <-assay(se.head) |> round()
mode(mat) <- "integer"

# remove batch
adjusted <- ComBat_seq(counts = mat, 
                       batch = se.head$batch, 
                       full_mod = T)

mode(adjusted) <- "integer"

# se is raw, se2 becomes the Combat adjusted data
se.head.adjusted <- se.head; assay(se.head.adjusted) <- adjusted

# ------------------------------------------------------------------------------
# set up dds's and run DESeq
# ------------------------------------------------------------------------------
dds_l <- list(gonad = se.gonad,
              head.raw = se.head,
              head = se.head.adjusted) |>
  map(~DESeqDataSet(.x,~knockdown2)) |>
  map(DESeq)

# ------------------------------------------------------------------------------
# differential expression on full set at once, raw and corrected
# using all possible contrasts with appropriate control
# ------------------------------------------------------------------------------
# function for extracting results
get_res <- function(x) {
  contrasts <- expand_grid(num = as.character(x$knockdown2),ctrl = as.character(x$knockdown2)) |>
    filter(str_detect(ctrl,"control") & !str_detect(num,"control")) |>
    filter(str_extract(num,"(?<=_).+") == str_extract(ctrl,"(?<=_).+")) |>
    mutate(contrast = map2(num,ctrl,~c("knockdown2",.x,.y))) |>
    distinct() |>
    pull(contrast) %>%
    set_names(.,map_chr(.,~paste(.x,collapse="_")))
  map(contrasts,~results(x,contrast = .x, format = "GRanges")) %>%
    map(~{mutate(.x,feature=names(.))})
}

# extract results (as GRanges)
res_l <- dds_l |> map(get_res)

# make a tibble of these results for easy access w/out using GRanges
combined_res_df <- res_l |> 
  map_df(~{map_df(.x,as_tibble,.id="comparison")},.id="group") |>
  mutate(comparison=str_extract(comparison,".+(?=_control)"))

#combined_res_df |>
#  filter(!str_detect(feature,"FBgn") & padj < 0.1) |>
#  dplyr::count(group,comparison,direction=sign(log2FoldChange)) |>
#  pivot_wider(names_from = direction, values_from = n)

# ------------------------------------------------------------------------------
# reorganize so this slots in with current downstream workflow (espec GSEA)
# ------------------------------------------------------------------------------
res2_l <- list(adjusted = c(res_l$gonad, res_l$head),
               raw = c(res_l$gonad, res_l$head.raw))

# ------------------------------------------------------------------------------
# export
# ------------------------------------------------------------------------------
saveRDS(res2_l,snakemake@output[["grs"]])
saveRDS(dds_l,snakemake@output[["dds"]])
saveRDS(combined_res_df,snakemake@output[["df"]])

