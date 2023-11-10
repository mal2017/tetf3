# ---------------------------------------------------------------------
# looking for association of motif presence and co-expression signal
# ---------------------------------------------------------------------

rule plot_motif_and_coex_on_tree:
    input:
        tidytree = "results/te_sequence_similarity/te_sketch_tidytree.rds",
        coex = config.get("MERGED_MODELS"),
        fimo = rules.fimo_denovo_motifs_tes.output.odir, 
    output:
        rds = "results/integrative/motif_and_coex_on_tree.{tf}.plot.rds",
    params:
        tf = "{tf}",
    script:
        "../scripts/integrative/plot_motif_and_coex_on_tree2.R"


rule n_denovo_vs_sig_coef:
    input:
        coex = config.get("MERGED_MODELS"),
        fimo = rules.fimo_denovo_motifs_tes.output.odir,
    output:
        rds = "results/integrative/n_denovo_vs_sig_coef.{tf}.rds",
    params:
        tf= "{tf}"
    script:
        "../scripts/integrative/n_denovo_vs_sig_coef.R"


    