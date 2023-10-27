rule pirna_enrichment_in_kd:
    input:
        pirna = rules.make_pirna_gene_list.output.tsv,
        res = rules.this_study_kd_deseq2.output.grs,
    output:
        de_pirna_tbl = "results/pirna_genes_de_in_kd_tbl.rds",
        de_pirna_fisher = "results/pirna_genes_de_in_kd_fisher.rds",
    script:
        "../scripts/pirna/pirna_in_kds.R"