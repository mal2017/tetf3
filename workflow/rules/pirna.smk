rule pirna_enrichment_in_kd:
    input:
        pirna = rules.make_pirna_gene_list.output.tsv,
        res = rules.this_study_kd_deseq2.output.grs,
    output:
        de_pirna_tbl = "results/pirna/pirna_genes_de_in_kd_tbl.rds",
        de_pirna_fisher = "results/pirna/pirna_genes_de_in_kd_fisher.rds",
    script:
        "../scripts/pirna/pirna_in_kds.R"

rule remap_peaks_dist_to_pirna:
    input:
        txdb = rules.make_txdb.output.txdb,
        pirna = rules.make_pirna_gene_list.output.tsv,
        remap = rules.annotate_fixed_insertions.output.remap,
    output:
        rds = "results/pirna/remap_peaks_dist_to_pirna.gr.rds",
    script:
        "../scripts/pirna/remap_peaks_near_pirna_genes.R"


rule collect_pirna_pathway_analysis:
    input:
        rules.pirna_enrichment_in_kd.output.de_pirna_fisher,
        rules.remap_peaks_dist_to_pirna.output.rds,
