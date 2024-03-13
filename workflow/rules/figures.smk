# ---------------------------------------------------------------------------------------------------
# figures related to results section 1
# ---------------------------------------------------------------------------------------------------
rule figure1:
    input:
        rules.plot_n_features_scatter.output,
        rules.plot_te_silencer_overrepresentation.output,
        rules.plot_te_silencer_scores_boxplot.output,
        rules.plot_te_silencer_n_hits_boxplot.output,
        rules.plot_n_models_per_filtering_step.output.gg,
    output:
        pdf="results/figures/panels/figure1_main.overview.pdf",
        xlsx="results/figures/data/figure1_main.overview.xlsx"
    params:
        figtitle="Figure 1"
    script:
        "../scripts/figures/overview.main.R"

rule figure1_supp_01:
    input:
        rules.plot_mf_shared_coex_hits.output,
        rules.independent_dataset_correlation.output.rds,
        rules.mf_dataset_correlation.output.rds,
        rules.plot_variance_explained.output.rds,
        rules.plot_te_silencer_scores_boxplot.output,
        rules.plot_te_silencer_n_hits_boxplot.output,
    output:
        pdf="results/figures/panels/figure1_supp_01.extra-overview.pdf",
        xlsx="results/figures/data/figure1_supp_01.extra-overview.xlsx"
    params:
        figtitle="Supplement 01 to Figure 1"
    script:
        "../scripts/figures/extra-overiew.R"

# ---------------------------------------------------------------------------------------------------
# figures related to results section 2
# ---------------------------------------------------------------------------------------------------

rule figure2:
    input:
        male_gg_gsea = "results/enrichment/sig_main_male_max_abs_estimate_qnorm.gg_gsea.rds",
        female_gg_gsea = "results/enrichment/sig_main_female_max_abs_estimate_qnorm.gg_gsea.rds",
    output:
        pdf="results/figures/panels/figure2_main.tf-overrepresentation.pdf",
        xlsx="results/figures/data/figure2_main.tf-overrepresentation.xlsx"
    params:
        figtitle="Figure 2"
    script:
        "../scripts/figures/tf-overrepresentation-in-coex.main.R"

# ---------------------------------------------------------------------------------------------------
# figures related to results section 3
# ---------------------------------------------------------------------------------------------------

rule figure3:
    input:
        rules.plot_ourKD_gsea_randomwalks.output,
        rules.plot_ourKD_gsea_barplots.output,
    output:
        pdf = "results/figures/panels/figure3_main.knockdowns.pdf"
    params:
        figtitle="Figure 3"
    script:
        "../scripts/figures/knockdowns.main.R"

rule figure3_supp_01:
    """
    factor specific kd enrichment
    """
    input:
        rules.plot_ourKD_gsea_randomwalks.output,
        rules.plot_ourKD_gsea_barplots.output,
    output:
        pdf="results/figures/panels/figure3_supp_01.factor_specific_gsea.pdf"
    params:
        figtitle="Supplement 01 to Figure 3"
    script:
        "../scripts/figures/our-kd-factor-specific-gsea.R"

# ---------------------------------------------------------------------------------------------------
# figures related to results section 4
# ---------------------------------------------------------------------------------------------------

rule figure4:
    input:
        "results/motifs/comparison/pan_denovo_comparison.meme.gg_df.rds",
        "upstream/csem_mosaics/bigwigs/"
    output:
       pdf="results/figures/panels/figure4_main.pan-motifs.pdf"
    params:
        figtitle="Figure 4"
    script:
      "../scripts/figures/pan-motifs.main.R"


rule figure4_supp_01:
    input:
        rules.csem_peaks_regioner.output,
    output:
        pdf="results/figures/panels/figure4_supp_01.csem-regioner.pdf",
    params:
        figtitle="Supplement 01 to Figure 4"
    script:
        "../scripts/figures/csem-regioner.R"


rule figure4_supp_02:
    input:
        n_denovo_vs_sig_coef = "results/motifs/n_denovo_vs_sig_coef.pan.rds",
        streme = "results/motifs/comparison/pan_denovo_comparison.streme.gg_df.rds",
        homer = "results/motifs/comparison/pan_denovo_comparison.homer.gg_df.rds",
        streme_empirical_fdr = "results/motifs/streme_per_tf_empirical_fdr/pan_empirical_fdr.tsv",
        streme_motif_comparison = "results/motifs/comparison/pan_denovo_comparison.streme.rds",
        streme_motif_similarity = "results/motifs/comparison/pan_denovo_similarity.streme.rds",
        homer_motif_comparison = "results/motifs/comparison/pan_denovo_comparison.homer.rds",
        homer_motif_similarity = "results/motifs/comparison/pan_denovo_similarity.homer.rds",
        denovo_motif_comparison = "results/motifs/comparison/pan_within_denovo_comparison.rds",
        denovo_motifs_um = "results/motifs/comparison/pan_within_denovo.universal_motif.rds",
        denovo_comparison_gg = "results/motifs/comparison/pan_within_denovo.gg.rds",
    output:
        pdf="results/figures/panels/figure4_supp_02.denovo-motif-pan-supp.pdf"
    params:
        figtitle="Supplement 02 to Figure 4"
    script:
        "../scripts/figures/denovo-motif-pan-supp.R"

rule figure4_supp_03:
    input:
        dds = "results/ripseq/unr_ripseq.dds.rds",
        tsv = "results/ripseq/unr_ripseq.tsv.gz",
        au = "results/ripseq/unr_ripseq_features_au_content.tsv.gz",
        au_in_region = "results/ripseq/unr_ripseq_features_au_content_in_region.tsv.gz",
        attta_in_region = "results/ripseq/unr_ripseq_features_attta_sites_in_region.tsv.gz",
        gsea = "results/ripseq/unr_bound_tx_in_kd.gsea.rds",
        p4d = "results/ripseq/unr_ripseq_phylosignal.p4d.rds",
    params:
        relpos = config.get("UNR_RIPSEQ_TX_RELATIVE_POSITION"),
        figtitle="Supplement 03 to Figure 4"
    output:
        pdf="results/figures/panels/figure4_supp_03.ripseq-ares.pdf",
    script:
        "../scripts/figures/unr_ripseq.R"


rule figure4_supp_04:
    input:
        rules.fca_heads_reanalysis.output,
    output:
        pdf="results/figures/panels/figure4_supp_04.fca.pdf",
    params:
        figtitle="Supplement 04 to Figure 4"
    script:
        "../scripts/figures/flycellatlas.R"

# ---------------------------------------------------------------------------------------------------
# figures related to results section 5
# ---------------------------------------------------------------------------------------------------


rule figure5:
    input:
        crlg = rules.filter_phylo_and_plot_correlograms.output.rds,
        control_crlg = rules.plot_control_correlograms.output.rds,
        tre = "results/te_sequence_similarity/te_sketch_tidytree.rds",
    output:
        pdf="results/figures/panels/figure5_main.phylosignal.pdf"
    params:
        figtitle="Figure 5"
    script:
        "../scripts/figures/phylosignal.main.R"

rule figure5_supp_01:
    input:
        "results/te_sequence_similarity/te_sketch_tidytree.rds",
    output:
        pdf="results/figures//panels/figure5_supp_01.dashing2-tree.pdf"
    params:
        figtitle="Supplement 01 to Figure 5"
    script:
        "../scripts/figures/dashing2-te-phylo-supplement-01.R"

rule figure5_supp_02:
    input:
        dds = "results/ripseq/unr_ripseq.dds.rds",
        tsv = "results/ripseq/unr_ripseq.tsv.gz",
        au = "results/ripseq/unr_ripseq_features_au_content.tsv.gz",
        au_in_region = "results/ripseq/unr_ripseq_features_au_content_in_region.tsv.gz",
        attta_in_region = "results/ripseq/unr_ripseq_features_attta_sites_in_region.tsv.gz",
        gsea = "results/ripseq/unr_bound_tx_in_kd.gsea.rds",
        p4d = "results/ripseq/unr_ripseq_phylosignal.p4d.rds",
    params:
        relpos = config.get("UNR_RIPSEQ_TX_RELATIVE_POSITION"),
        figtitle="Supplement 03 to Figure 5"
    output:
        pdf="results/figures/panels/figure5_supp_02.ripseq-phylosignal.pdf"
    script:
        "../scripts/figures/unr_ripseq_phylosignal.R"

rule figure5_supp_03:
    input:
        rules.plot_de_volcanos.output,
    output:
        pdf="results/figures/panels/figure5_supp_03.pirna-volcano.pdf"
    params:
        figtitle="Supplement 03 to Figure 5"
    script:
        "../scripts/figures/pirna-and-tes-in-kd-volcano.R"

# ---------------------------------------------------------------------------------------------------
# figures related to minor details in methods
# ---------------------------------------------------------------------------------------------------

rule methods_figs:
    input:
        rules.plot_check_kds_by_chip_prox.output,
    output:
        pdf="results/figures/panels/methods_supp_01.pdf"
    params:
        figtitle="Supplement 01 to Methods"
    script:
        "../scripts/figures/methods_supp_01.R"


# ---------------------------------------------------------------------------------------------------
# bring it together
# ---------------------------------------------------------------------------------------------------

rule figures:
    input:
        rules.figure1.output.pdf,
        rules.figure1_supp_01.output.pdf,
        rules.figure2.output.pdf,
        rules.figure3.output.pdf,
        rules.figure3_supp_01.output.pdf,
        rules.figure4.output.pdf,
        rules.figure4_supp_01.output.pdf,
        rules.figure4_supp_02.output.pdf,
        rules.figure4_supp_03.output.pdf,
        rules.figure4_supp_04.output.pdf,
        rules.figure5.output.pdf,
        rules.figure5_supp_01.output.pdf,
        rules.figure5_supp_02.output.pdf,
        rules.figure5_supp_03.output.pdf,
        rules.methods_figs.output.pdf,
    output:
        "results/figures/figures.pdf"
    shell:
        """
        gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile={output} {input}
        """

rule figuredata:
    input:
        rules.figure1.output.xlsx,
        rules.figure1_supp_01.output.xlsx,
        rules.figure2.output.xlsx,
    output:
        xlsx="results/figures/figuredata.xlsx"
    script:
        "../scripts/figures/combine_figuredata.R"