# For finding TF-coexpressed TE signatures in knockdown data

rule ourKD_gsea:
    input:
        mods = config.get("MERGED_MODELS"),
        res = rules.this_study_kd_deseq2.output.grs,
        pirna = rules.make_pirna_gene_list.output.tsv,
    output:
        rds = "results/signatures/ourKD_gsea.rds",
    script:
        "../scripts/signatures/ourKD_gsea.R"


rule s2rplus_coex_te_gsea_by_de:
    input:
        deg = rules.s2rplus_limma.output.tsv,
        coex = config.get("MERGED_MODELS")
    output:
        rds = "results/signatures/s2rplus_te_gsea.rds",
    script:
        "../scripts/signatures/tfrnai_gsea_de.R"

rule signatures:
    input:
        rules.ourKD_gsea.output,
        rules.s2rplus_coex_te_gsea_by_de.output