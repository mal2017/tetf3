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

rule plot_ourKD_gsea_randomwalks:
    input:
        rds = rules.ourKD_gsea.output.rds,
    output:
        gg_df = "results/signatures/ourKD_gsea_randomwalks.gg_df.rds",
    script:
        "../scripts/signatures/plot_ourKD_gsea_randomwalks.R"

rule plot_ourKD_gsea_barplots:
    input:
        rds = rules.ourKD_gsea.output.rds,
    output:
        gg_list = "results/signatures/ourKD_gsea_barplots.gg_list.rds",
    script:
        "../scripts/signatures/plot_ourKD_gsea_barplots.R"