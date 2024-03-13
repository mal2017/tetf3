rule table_te_regulator_chip_prox:
    input:
        "results/resources/pirna_pathway.tsv",
        "results/resources/gene_symbol_lookup.tsv.gz",
        "results/pirna/encode_peaks_dist_to_pirna.gr.rds",
    output:
        docx= "results/tables/table_te_regulator_chip_prox.docx",
    script:
        "../scripts/tables/known_TE_regulators_in_chips_table.v2.R"

rule tables:
    input:
        "results/tables/table_te_regulator_chip_prox.docx",