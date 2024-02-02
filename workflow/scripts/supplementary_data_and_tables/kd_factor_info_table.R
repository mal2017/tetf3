library(tidyverse)
library(gt)

# Bgee's expression presence absence calls
bgee_calls <- read_tsv("resources/Bgee_expression_calls.Drosophila_melanogaster_expr_advanced_all_conditions.tsv.gz")

expr_df <- bgee_calls |>
  #filter(`Gene name` %in% c("CG16779")) |>
  filter(`Gene name` %in% c("pan","NfI","CG16779","Unr","vvl")) |>
  filter(Expression=="present") |>
  filter(Strain %in% c("wild-type")) |>
  filter(`Developmental stage name` %in% c("fully formed stage","life cycle")) |>
  filter(str_starts(`Anatomical entity ID`,"UBERON")) |> # simplify to uberon terms
  filter((`Anatomical entity name` %in% c("testis","ovary","female reproductive system")) | (`Anatomical entity name` %in% c("insect adult head","head") & Sex == "any")) |> # simplify to broad body parts
  dplyr::select(`Gene symbol`=`Gene name`,`Gene ID`,tissue=`Anatomical entity name`,Expression,FDR) |>
  mutate(Expression = map2_chr(Expression,FDR,~sprintf("%s (%s)",.x,format.pval(.y)))) |>
  dplyr::select(-FDR) |>
  pivot_wider(names_from = "tissue", values_from = Expression)

# tf family
atfdb <- read_tsv("resources/Drosophila_melanogaster_TF.txt") |>
  dplyr::select(`Gene symbol`=Symbol,Family) |>
  distinct()

# factor description
desc_df <- read_tsv("https://ftp.flybase.net/releases/FB2023_06/precomputed_files/genes/gene_snapshots_fb_2023_06.tsv.gz",
         skip=4,col_names = c("Gene ID",'Gene Symbol','Gene Name',"datestamp","snapshot")) |>
  filter(`Gene ID` %in% expr_df$`Gene ID`)

desc_df <- desc_df |>
  mutate(snapshot=str_extract(snapshot,"encodes a.+")) |>
  mutate(snapshot=replace_na(snapshot,"")) |>
  dplyr::select(`Gene ID`,`Gene symbol`=`Gene Symbol`,`Flybase snapshot`=snapshot)

tab <- left_join(expr_df, atfdb) |>
  left_join(desc_df) |>
  arrange(desc(`Flybase snapshot`)) |>
  mutate(`Flybase snapshot` = str_extract(`Flybase snapshot`,".+?(?=\\.)")) |>
  dplyr::relocate(c(Family,`Flybase snapshot`),.after = `Gene ID`) |>
  gt() |>
  tab_spanner("Bgee expression call (FDR)",columns = c("testis","female reproductive system","insect adult head")) |>
  tab_style(locations = cells_body(columns="Flybase snapshot"),style = cell_text(style="italic",size = "xx-small")) |>
  tab_style(locations = cells_body(columns=c("testis","female reproductive system","insect adult head")),style = cell_text(size = "x-small")) |>
  cols_width(`Flybase snapshot` ~ px(150))


gtsave(tab,"results/tables/kd_factor_info.table.docx")
