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
        pdf="results/figures/figure1_main.overview.pdf"
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
        pdf="results/figures/figure1_supp_01.extra-overview.pdf"
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
        pdf="results/figures/figure2_main.tf-overrepresentation.pdf"
    script:
        "../scripts/figures/tf-overrepresentation-in-coex.main.R"

# # ---------------------------------------------------------------------------------------------------
# # figures related to results section 3
# # ---------------------------------------------------------------------------------------------------

rule figure3:
    input:
        rules.plot_ourKD_gsea_randomwalks.output,
        rules.plot_ourKD_gsea_barplots.output,
    output:
        pdf = "results/figures/figure3_main.knockdowns.pdf"
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
        pdf="results/figures/figure3_supp_01.factor_specific_gsea.pdf"
    script:
        "../scripts/figures/our-kd-factor-specific-gsea.R"

rule figure3_supp_02:
    input:
        rules.plot_de_volcanos.output,
    output:
        pdf="results/figures/figure3_supp_02.pirna-volcano.pdf"
    script:
        "../scripts/figures/pirna-and-tes-in-kd-volcano.R"


# # ---------------------------------------------------------------------------------------------------
# # figures related to results section 4
# # ---------------------------------------------------------------------------------------------------

rule figure4:
    input:
        "results/motifs/comparison/pan_denovo_comparison.meme.gg_df.rds",
        "upstream/csem_mosaics/bigwigs/"
    output:
       pdf="results/figures/figure4_main.pan-motifs.pdf"
    script:
      "../scripts/figures/pan-motifs.main.R"

rule figure4_supp_01:
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
        pdf="results/figures/figure4_supp_01.denovo-motif-pan-supp.pdf"
    script:
        "../scripts/figures/denovo-motif-pan-supp.R"


rule figure4_supp_02:
    """
    pan tracks and qc
    """
    input:
        "upstream/csem_mosaics/bigwigs/",
        rules.plot_quality_by_visual_pericent_inspection_status.output,
        rules.plot_quality_by_visual_pericent_inspection_status.output,
    output:
        pdf="results/figures/figure4_supp_02.csem-tracks-and-qc-pan-profile.pdf",
    script:
        "../scripts/figures/csem-tracks-and-qc-pan-profile.R"


# rule figure4_supp_01:
#     input:
#         n_denovo_vs_sig_coef = "results/integrative/n_denovo_vs_sig_coef.pan.rds",
#         pan_tree = "results/integrative/motif_and_coex_on_tree.pan.plot.rds",
#     output:
#         pdf="results/figures/figure4_supp_01.pan-meme-denovo.pdf"
#     script:
#         "../scripts/figures/pan-meme-denovo.R"

# rule figure4_supp_02:
#     input:
#         motif_comparison = "results/motifs/comparison/pan_denovo_comparison.homer.rds",
#         motif_similarity = "results/motifs/comparison/pan_denovo_similarity.homer.rds",
#     output:
#         pdf="results/figures/figure4_supp_02.pan-homer-denovo.pdf"
#     script:
#         "../scripts/figures/pan-homer-denovo.R"


# rule figure4_supp_03:
#     input:
#         denovo_empirical_fdr = "results/motifs/streme_per_tf_empirical_fdr/pan_empirical_fdr.tsv",
#         motif_comparison = "results/motifs/comparison/pan_denovo_comparison.streme.rds",
#         motif_similarity = "results/motifs/comparison/pan_denovo_similarity.streme.rds",
#     output:
#         pdf="results/figures/figure4_supp_03.pan-streme-denovo.pdf"
#     script:
#         "../scripts/figures/pan-streme-denovo.R"



# rule figure4_supp_07:
#     input:
#         rules.csem_peaks_regioner.output,
#     output:
#         pdf="results/figures/figure4_supp_07.csem-regioner.pdf",
#     script:
#         "../scripts/figures/csem-regioner.R"



# rule figure5_and_supp:
#     input:
#         dds = "results/ripseq/unr_ripseq.dds.rds",
#         tsv = "results/ripseq/unr_ripseq.tsv.gz",
#         au = "results/ripseq/unr_ripseq_features_au_content.tsv.gz",
#         au_in_region = "results/ripseq/unr_ripseq_features_au_content_in_region.tsv.gz",
#         attta_in_region = "results/ripseq/unr_ripseq_features_attta_sites_in_region.tsv.gz",
#         gsea = "results/ripseq/unr_bound_tx_in_kd.gsea.rds",
#         p4d = "results/ripseq/unr_ripseq_phylosignal.p4d.rds",
#     params:
#         relpos = config.get("UNR_RIPSEQ_TX_RELATIVE_POSITION")
#     output:
#         pdf="results/figures/figure5.ripseq-ares.pdf",
#         pdf2="results/figures/figure5_supp_01.ripseq-phylosignal.pdf"
#     script:
#         "../scripts/figures/unr_ripseq.R"

# rule figure5_supp_02:
#     input:
#         rules.fca_heads_reanalysis.output,
#     output:
#         pdf="results/figures/figure5_supp_02.fca.pdf",
#     script:
#         "../scripts/figures/flycellatlas.R"

# rule figure5_supp_03:
#     input:
#         rules.calderon22_reanalysis.output,
#     output:
#         pdf="results/figures/figure5_supp_03.calderon22.pdf",
#     script:
#         "../scripts/figures/calderon22.R"

# # # ---------------------------------------------------------------------------------------------------
# # # figures related to results section 5
# # # ---------------------------------------------------------------------------------------------------


# rule figure6:
#     input:
#         crlg = rules.filter_phylo_and_plot_correlograms.output.rds,
#         control_crlg = rules.plot_control_correlograms.output.rds,
#         tre = "results/te_sequence_similarity/te_sketch_tidytree.rds",
#     output:
#         pdf="results/figures/figure6_main.phylosignal.pdf"
#     script:
#         "../scripts/figures/phylosignal.main.R"

# rule figure6_supp_01:
#     input:
#         "results/te_sequence_similarity/te_sketch_tidytree.rds",
#     output:
#         pdf="results/figures/figure6_supp_01.dashing2-tree.pdf"
#     script:
#         "../scripts/figures/dashing2-te-phylo-supplement-01.R"



# ---------------------------------------------------------------------------------------------------
# bring it together
# ---------------------------------------------------------------------------------------------------

rule figures:
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

