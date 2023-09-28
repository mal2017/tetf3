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
        rds = "results/phylosignal/correlograms.rds",
    script:
        "../scripts/phylosignal/plot_control_correlograms.R"

rule plot_sig_correlograms:
    input:
        phylosignal = rules.phylosignal_te_tree.output.phylosignal,
        df = rules.phylosignal_te_tree.output.tab,
    output:
        rds = "results/phylosignal/sig_correlograms.rds",
    script:
        "../scripts/phylosignal/plot_sig_correlograms.R"