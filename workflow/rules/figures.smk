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
        "../scripts/figures/extra-overview.R"


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
        pdf = "results/figures/panels/figure3_main.knockdowns.pdf",
        xlsx = "results/figures/data/figure3_main.knockdowns.xlsx"
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
        rules.fca_heads_reanalysis.output,
    output:
        pdf="results/figures/panels/figure3_supp_01.factor_specific_gsea.pdf",
        xlsx="results/figures/data/figure3_supp_01.factor_specific_gsea.xlsx"
    params:
        figtitle="Supplement 01 to Figure 3"
    script:
        "../scripts/figures/our-kd-factor-specific-gsea.R"

rule figure3_supp_02:
    input:
        rules.fca_heads_reanalysis.output,
        rules.pan_unr_summed_te_counts.output,
    output:
        pdf="results/figures/panels/figure3_supp_02.fca.pdf",
        xlsx="results/figures/data/figure3_supp_02.fca.xlsx"
    params:
        figtitle="Supplement 02 to Figure 3"
    script:
        "../scripts/figures/flycellatlas_01.R"


# ---------------------------------------------------------------------------------------------------
# figures related to results section 4
# ---------------------------------------------------------------------------------------------------

rule figure4:
    input:
        "results/motifs/comparison/pan_denovo_comparison.meme.gg_df.rds",
        "results/ripseq/unr_ripseq_maplot.gg.rds",
        "results/ripseq/unr_ripseq_au_richness.boxplot.gg_list.rds",
        "results/ripseq/unr_ripseq_attta_richness.boxplot.gg_list.rds",
    output:
       pdf="results/figures/panels/figure4_main.pan-motifs.pdf",
       xlsx="results/figures/data/figure4_main.pan-motifs.xlsx"
    params:
        figtitle="Figure 4"
    script:
      "../scripts/figures/pan-motifs.main.R"


rule figure4_supp_01:
    input:
        rules.csem_peaks_regioner.output,
    output:
        pdf="results/figures/panels/figure4_supp_01.csem-regioner.pdf",
        xlsx="results/figures/data/figure4_supp_01.csem-regioner.xlsx",
    params:
        figtitle="Supplement 01 to Figure 4",
    script:
        "../scripts/figures/csem-regioner.R"


rule figure4_supp_02:
    input:
        n_denovo_vs_sig_coef = "results/motifs/n_denovo_vs_sig_coef.pan.rds",
        #streme_empirical_fdr = "results/motifs/streme_per_tf_empirical_fdr/pan_empirical_fdr.tsv",
        streme_motif_comparison = "results/motifs/comparison/pan_denovo_comparison.streme.rds",
        homer_motif_comparison = "results/motifs/comparison/pan_denovo_comparison.homer.rds",
        denovo_motif_comparison = "results/motifs/comparison/pan_within_denovo_comparison.rds",
        denovo_motifs_um = "results/motifs/comparison/pan_within_denovo.universal_motif.rds",
        denovo_comparison_gg = "results/motifs/comparison/pan_within_denovo.gg.rds",
        sea_csem = 'results/motifs/upstream_csem_known_pan_sea.gg.rds',
    output:
        pdf="results/figures/panels/figure4_supp_02.denovo-motif-pan-supp.pdf",
        xlsx="results/figures/data/figure4_supp_02.denovo-motif-pan-supp.xlsx",
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
        au_richness = rules.plot_au_richness.output.lineplot,
        atta_richness = rules.plot_attta_richness.output.boxplot,
    params:
        relpos = config.get("UNR_RIPSEQ_TX_RELATIVE_POSITION"),
        figtitle="Supplement 03 to Figure 4"
    output:
        pdf="results/figures/panels/figure4_supp_03.ripseq-ares.pdf",
        xlsx="results/figures/data/figure4_supp_03.ripseq-ares.xlsx",
    script:
        "../scripts/figures/unr_ripseq.R"


rule figure4_supp_04:
    input:
        rules.fca_heads_reanalysis.output,
    output:
        pdf="results/figures/panels/figure4_supp_04.fca.pdf",
        xlsx="results/figures/data/figure4_supp_04.fca.xlsx"
    params:
        figtitle="Supplement 04 to Figure 4"
    script:
        "../scripts/figures/flycellatlas_02.R"

# ---------------------------------------------------------------------------------------------------
# figures related to results section 5
# ---------------------------------------------------------------------------------------------------


rule figure5:
    input:
        #crlg = rules.filter_phylo_and_plot_correlograms.output.rds,
        control_crlg = rules.plot_control_correlograms.output.rds,
        tre = "results/te_sequence_similarity/te_sketch_tidytree.rds",
        crlgs = rules.plot_main_fig_correlograms.output.rds,
        ripseq_crlg = "results/ripseq/unr_ripseq_phylosignal.crlg.gg.rds",
    output:
        pdf="results/figures/panels/figure5_main.phylosignal.pdf",
        xlsx="results/figures/data/figure5_main.phylosignal.xlsx",
    params:
        figtitle="Figure 5"
    script:
        "../scripts/figures/phylosignal.main.R"

rule figure5_supp_01:
    input:
        "results/te_sequence_similarity/te_sketch_tidytree.rds",
    output:
        pdf="results/figures/panels/figure5_supp_01.dashing2-tree.pdf",
        xlsx="results/figures/data/figure5_supp_01.dashing2-tree.xlsx",
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
        tree = "results/ripseq/unr_ripseq_phylosignal.tree.gg.rds"
    params:
        relpos = config.get("UNR_RIPSEQ_TX_RELATIVE_POSITION"),
        figtitle="Supplement 02 to Figure 5"
    output:
        pdf="results/figures/panels/figure5_supp_02.ripseq-phylosignal.pdf",
        xlsx="results/figures/data/figure5_supp_02.ripseq-phylosignal.xlsx",
    script:
        "../scripts/figures/unr_ripseq_phylosignal.R"

# rule figure5_supp_03:
#     input:
#         rules.plot_de_volcanos.output,
#     output:
#         pdf="results/figures/panels/figure5_supp_03.pirna-volcano.pdf",
#         xlsx="results/figures/data/figure5_supp_03.pirna-volcano.xlsx",
#     params:
#         figtitle="Supplement 03 to Figure 5"
#     script:
#         "../scripts/figures/pirna-and-tes-in-kd-volcano.R"

rule figure5_supp_03:
    input:
        rds = rules.plot_ourKD_gsea_randomwalks.output.gg_df,
    output:
        pdf="results/figures/panels/figure5_supp_04.knockdowns_silencer_enrichment.pdf",
        xlsx="results/figures/data/figure5_supp_04.knockdowns_silencer_enrichment.xlsx",
    params:
        figtitle="Supplement 03 to Figure 5"
    script:
        "../scripts/figures/knockdowns.silencer_enrichment.R"

# ---------------------------------------------------------------------------------------------------
# figures related to minor details in methods
# ---------------------------------------------------------------------------------------------------

# rule methods_figs:
#     input:
#         rules.plot_check_kds_by_chip_prox.output,
#     output:
#         pdf="results/figures/panels/methods_supp_01.pdf"
#     params:
#         figtitle="Supplement 01 to Methods"
#     script:
#         "../scripts/figures/methods_supp_01.R"


# ---------------------------------------------------------------------------------------------------
# bring it together
# ---------------------------------------------------------------------------------------------------

rule collect_figures:
    input:
        rules.figure1.output,
        rules.figure1_supp_01.output,
        rules.figure2.output,
        rules.figure3.output,
        rules.figure3_supp_01.output,
        rules.figure3_supp_02.output,
        rules.figure4.output,
        rules.figure4_supp_01.output,
        rules.figure4_supp_02.output,
        rules.figure4_supp_03.output,
        rules.figure4_supp_04.output,
        rules.figure5.output,
        rules.figure5_supp_01.output,
        rules.figure5_supp_02.output,
        rules.figure5_supp_03.output,


# rule figures:
#     input:
#         rules.figure1.output.pdf,
#         rules.figure1_supp_01.output.pdf,
#         rules.figure2.output.pdf,
#         rules.figure3.output.pdf,
#         rules.figure3_supp_01.output.pdf,
#         rules.figure3_supp_02.output.pdf,
#         rules.figure4.output.pdf,
#         rules.figure4_supp_01.output.pdf,
#         rules.figure4_supp_02.output.pdf,
#         rules.figure4_supp_03.output.pdf,
#         rules.figure4_supp_04.output.pdf,
#         rules.figure5.output.pdf,
#         rules.figure5_supp_01.output.pdf,
#         rules.figure5_supp_02.output.pdf,
#         rules.figure5_supp_03.output.pdf,
#     output:
#         "results/figures/figures.pdf"
#     shell:
#         """
#         for f in {input}; do convert -trim -density 600 -quality 100 "$f" "${{f%pdf}}jpg"; done
#         gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile={output} {input}
#         """

# rule supp_figures:
#     input:
#         rules.figure1_supp_01.output.pdf,
#         rules.figure3_supp_01.output.pdf,
#         rules.figure3_supp_02.output.pdf,
#         rules.figure4_supp_01.output.pdf,
#         rules.figure4_supp_02.output.pdf,
#         rules.figure4_supp_03.output.pdf,
#         rules.figure4_supp_04.output.pdf,
#         rules.figure5_supp_01.output.pdf,
#         rules.figure5_supp_02.output.pdf,
#         rules.figure5_supp_03.output.pdf,
#     output:
#         "results/figures/supp_figures.pdf"
#     shell:
#         """
#         gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile={output} {input}
#         """

# rule main_figures:
#     input:
#         rules.figure1.output.pdf,
#         rules.figure2.output.pdf,
#         rules.figure3.output.pdf,
#         rules.figure4.output.pdf,
#         rules.figure5.output.pdf,
#     output:
#         "results/figures/main_figures.pdf"
#     shell:
#         """
#         gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile={output} {input}
#         """


rule figuredata:
    input:
        rules.figure1.output.xlsx,
        rules.figure1_supp_01.output.xlsx,
        rules.figure2.output.xlsx,
        rules.figure3.output.xlsx,
        rules.figure3_supp_01.output.xlsx,
        rules.figure3_supp_02.output.xlsx,
        rules.figure4.output.xlsx,
        rules.figure4_supp_01.output.xlsx,
        rules.figure4_supp_02.output.xlsx,
        rules.figure4_supp_03.output.xlsx,
        rules.figure4_supp_04.output.xlsx,
        rules.figure5.output.xlsx,
        rules.figure5_supp_01.output.xlsx,
        rules.figure5_supp_02.output.xlsx,
        rules.figure5_supp_03.output.xlsx,
    output:
        xlsx="results/figures/figuredata.xlsx"
    script:
        "../scripts/figures/combine_figuredata.R"