library(tidyverse)
library(DESeq2)
library(Biostrings)
library(sva)
library(plyranges)

#se_path <- "upstream/kd.se.gene.0.rds"
se_path <- snakemake@input[["se"]]

gene_universe <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2022_04.tsv.gz", skip = 4)

allowed_genes <- gene_universe %>% filter(gene_type %in% c("protein_coding_gene"))

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()

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

se <- se[,se$driver %in% c("aTub","tj","Mef2.R")]

# ----------- count filtering --------------------------------------------------
keep <- rowSums(assay(se,"counts") >= 10) >= 2

se <- se[keep,]

# annotation -------------------------------------------------------------------

elementMetadata(se) <- elementMetadata(se) %>%
  as_tibble() %>%
  left_join(lkup,by=c(gene_id="gene_ID")) %>%
  mutate(gene_symbol = if_else(is.na(gene_symbol),gene_id,gene_symbol)) %>%
  DataFrame()

# adjust for batch ---------- --------------------------------------------------
assay(se) <- assay(se) %>% round()
mat <-assay(se)

mode(mat) <- "integer"

# remove batch
# adding extra covar via covar_mod paramter seems to skyrocket counts of certain genes, 
# ie to as much as 1e12
# see https://github.com/zhangyuqing/ComBat-seq/issues/20
# https://github.com/zhangyuqing/ComBat-seq/issues/15
# this seems seems to be related to using covar_mod, so we won't do that
# https://www.biostars.org/p/456699/
adjusted <- ComBat_seq(counts = mat, 
                       batch = se$batch,
                       group = se$knockdown,full_mod = T)

mode(adjusted) <- "integer"

# se is raw, se2 becomes the Combat adjusted data
se2 <- se; assay(se2) <- adjusted

# differential expression on full set at once, raw and corrected----------------

ses <- list(raw = se, adjusted = se2)

# set up the contrasts we want to use
contrasts <- list(CG16779.head = c("knockdown2","CG16779_female_head_Mef2.R","control_female_head_Mef2.R"),
                  Unr = c("knockdown2","Unr_female_head_Mef2.R","control_female_head_Mef2.R"),
                  NfI = c("knockdown2","NfI_female_head_Mef2.R","control_female_head_Mef2.R"),
                  vvl = c("knockdown2","vvl_female_head_Mef2.R","control_female_head_Mef2.R"),
                  pan.testis = c("knockdown2","pan_male_gonad_aTub","control_male_gonad_aTub"),
                  pan.tj.ovary = c("knockdown2","pan_female_gonad_tj","control_female_gonad_tj"),
                  CG16779.testis = c("knockdown2","CG16779_male_gonad_aTub","control_male_gonad_aTub"),
                  CG16779.tj.ovary = c("knockdown2","CG16779_female_gonad_tj","control_female_gonad_tj"),
                  pan.head = c("knockdown2","pan_female_head_Mef2.R","control_female_head_Mef2.R"))

contrasts <- contrasts %>% set_names(.,map_chr(.,paste,collapse="_"))

# function for creating dds
run_deseq <- . %>%
  estimateSizeFactors(type="poscounts") %>%
  DESeq()

# create dds's
dds_list <- map(ses,~DESeqDataSet(.,~knockdown2)) %>%
  map(run_deseq)

# function for extracting results
get_res <- function(dds) {
  map(contrasts,~lfcShrink(dds,contrast = .x,type = "normal", format = "GRanges")) %>%
    map(~{mutate(.x,feature=names(.))})
}

# extract results (as GRanges)
res_list <- map(dds_list,get_res) 

# make a tibble of these results for easy access w/out using GRanges
combined_res_df <- res_list |> 
  map_df(~{map_df(.x,as_tibble,.id="comparison")},.id="adjustment") |>
  mutate(comparison=str_extract(comparison,".+(?=_control)"))

saveRDS(res_list,snakemake@output[["grs"]])
saveRDS(dds_list,snakemake@output[["dds"]])
saveRDS(combined_res_df,snakemake@output[["df"]])

