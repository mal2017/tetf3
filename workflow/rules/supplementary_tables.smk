rule table_te_regulator_chip_prox:
    input:
        "results/resources/pirna_pathway.tsv",
        "results/resources/gene_symbol_lookup.tsv.gz",
        "results/pirna/encode_peaks_dist_to_pirna.gr.rds",
    output:
        docx= "results/tables/table_te_regulator_chip_prox.docx",
    script:
        "../scripts/tables/known_TE_regulators_in_chips_table.v2.R"

rule table_te_regulator_in_kds:
    input:
        "results/resources/pirna_pathway.tsv",
        "results/resources/gene_symbol_lookup.tsv.gz",
        rules.this_study_kd_deseq2.output.grs,
    output:
        xlsx= "results/tables/table_known_te_regulators_in_knockdown.xlsx",
    script:
        "../scripts/tables/known_TE_regulators_in_knockdowns_table.R"

rule motif_results_table:
    input:
        expand("results/motifs/streme_per_tf_empirical_fdr/{tf}_empirical_fdr.tsv",tf=TFSOI),
        expand("results/motifs/fimo_on_tes/denovo/{tf}",tf="pan"),
        "results/motifs/csem_peak_sea.known.pan.tsv.gz",
    output:
        xlsx= "results/tables/motif_results.xlsx",
    script:
        "../scripts/tables/motif_enrichment_results.R"


rule tables:
    input:
        "results/tables/table_te_regulator_chip_prox.docx",
        "results/tables/table_known_te_regulators_in_knockdown.xlsx",
        "results/tables/descriptive_lms.xlsx",
        "results/tables/motif_results.xlsx",