
rule plot_pirna_genes_in_our_kd:
    input:
        res = rules.this_study_kd_deseq2.output.grs,
        pirna = rules.make_pirna_gene_list.output.tsv
    output:
        rds = "results/other_plots/pirna_genes_in_our_kd.rds",
        png = "results/other_plots/pirna_genes_in_our_kd.png",
    script:
        "../scripts/other_plots/plot_piRNA_genes_in_our_kd.R"

rule plot_pirna_genes_in_our_kd_all:
    input:
        res = rules.this_study_kd_deseq2.output.grs,
        pirna = rules.make_pirna_gene_list.output.tsv
    output:
        rds = "results/other_plots/pirna_genes_in_our_kd_all.rds",
        png = "results/other_plots/pirna_genes_in_our_kd_all.png",
    script:
        "../scripts/other_plots/plot_piRNA_genes_in_our_kd_all.R"