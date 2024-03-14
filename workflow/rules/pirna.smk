rule plot_te_silencer_overrepresentation:
    input:
        mods = config.get("MERGED_MODELS"),
        pirna = rules.make_pirna_gene_list.output.tsv
    output:
        rds = "results/pirna/te_silencers_in_lms.rds",
        gg = "results/pirna/te_silencer_overrepresentation.gg.rds"
    script:
        "../scripts/pirna/plot_te_silencer_overrepresentation.R"        

rule plot_te_silencer_n_hits_boxplot:
    input:
        mods = rules.plot_te_silencer_overrepresentation.output.rds,
    output:
        female = "results/pirna/te_silencer_n_hits_boxplot.females.gg.rds",
        male = "results/pirna/te_silencer_n_hits_boxplot.males.gg.rds"
    script:
        "../scripts/pirna/plot_te_silencer_nhits_boxplot.R"

rule plot_te_silencer_scores_boxplot:
    input:
        rules.rankings.output,
        teregs = rules.make_pirna_gene_list.output.tsv,
    output:
        female = "results/pirna/te_silencer_scores_boxplot.females.gg.rds",
        male = "results/pirna/te_silencer_scores_boxplot.males.gg.rds"
    script:
        "../scripts/pirna/plot_te_silencer_scores_boxplot.R"

rule encode_peaks_dist_to_pirna:
    input:
        txdb = rules.make_txdb.output.txdb,
        pirna = rules.make_pirna_gene_list.output.tsv,
    output:
        rds = "results/pirna/encode_peaks_dist_to_pirna.gr.rds",
    script:
        "../scripts/pirna/encode_peaks_dist_to_piRNA_genes.R"


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
