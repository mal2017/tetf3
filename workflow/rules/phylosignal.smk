rule phylosignal_te_tree:
    input:
        zads = rules.get_zad_genes.output.tsv,
        tfs  = config.get("TFS"),
        tree = "results/te_sequence_similarity/te_sketch_tidytree.rds",
        mods = config.get("MERGED_MODELS"),
    output:
        phylosignal = "results/phylosignal/phylosignal.rds",
        tab = "results/phylosignal/phylosignal_df.rds",
    script:
        "../scripts/phylosignal/phylosignal.R"

rule plot_control_correlograms:
    input:
        phylosignal = rules.phylosignal_te_tree.output.phylosignal,
    output:
        rds = "results/phylosignal/control_correlograms.rds",
    threads:
        4
    script:
        "../scripts/phylosignal/plot_control_correlograms.R"

rule filter_phylo_and_plot_correlograms:
    input:
        phylosignal = rules.phylosignal_te_tree.output.phylosignal,
        df = rules.phylosignal_te_tree.output.tab,
    threads:
        4   
    output:
        rds = "results/phylosignal/goi_correlograms.rds",
        tsv = "results/phylosignal/phylosignal_filtered_hits.tsv.gz",
    script:
        "../scripts/phylosignal/filter_phylo_and_plot_correlograms.R"

rule plot_main_fig_correlograms:
    input:
        crlg = rules.filter_phylo_and_plot_correlograms.output.rds,
        filt = rules.filter_phylo_and_plot_correlograms.output.tsv,
    output:
        rds = "results/phylosignal/main_fig_correlograms.gg.rds",
    script:
        "../scripts/phylosignal/plot_main_fig_correlograms.R"
