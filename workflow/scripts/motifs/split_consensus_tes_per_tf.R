Sys.setenv(R_PROFILE=".Rprofile")
source(Sys.getenv("R_PROFILE"))


library(tidyverse)
library(rtracklayer)

tfs <- ifelse(exists("snakemake"),snakemake@input[["tfs"]],
               "resources/Drosophila_melanogaster_TF.txt") %>%
  read_tsv()

lms0 <- ifelse(exists("snakemake"),snakemake@input[["mods"]],
    "upstream/final-models.collected-info.tsv.gz") %>% 
    read_tsv() %>%
    filter(significant_x) #%>%
    #filter(feature.x %in% tfs$Ensembl)

tes_2_use_0 <- lms0 %>% 
  dplyr::select(gene_symbol,feature.y) %>%
  distinct() %>%
  group_by(gene_symbol) %>%
  summarise(n = n()) |>
  filter(n >=10) |>
  left_join(lms0) |>
  pull(feature.y) |>
  unique()

lms <- lms0 |>
    filter(gene_symbol %in% c("pan","CG16779","vvl","Unr","NfI"))

fa <- ifelse(exists("snakemake"),snakemake@input[["tes"]],
        "resources/Tidalbase_transposon_sequence.fasta") %>%
        import(format="fasta")

# restrict to TEs that are coexpressed with factors
# coexpressed with many TEs - this is presumably
# enriched for non-sequence-specific TE/gene interactions and
# therefore a reasonable background set  for comparing putatively
# TF-regulated groups of TEs to groups that may not be TF-regulated
fa <- fa[names(fa) %in% tes_2_use_0]

seqs <- lms %>% 
    dplyr::select(gene_symbol,feature.y) %>%
    distinct() %>%
    split(.,.$gene_symbol) %>%
    map(pull,feature.y) %>%
    map(.f = ~{list(coex = fa[.x],other = fa[!names(fa) %in% .x])})

seqs <- seqs[c("pan","CG16779","vvl","Unr","NfI")]

odir0 <- ifelse(exists("snakemake"),snakemake@output[["odir"]],
               "results/analysis/motifs/consensus_tes_per_tf/")

for (s in names(seqs)) {
    message(s)
    sq <- seqs[[s]]
    
    for (k in names(sq)) {
        odir <- sprintf("%s/%s",odir0,s)
        dir.create(odir,recursive = T)
        export(sq[[k]],paste0(odir,"/",k,".fasta"))
    }
}
