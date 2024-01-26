# ---------------------------------------------------------------------------------------------------
# figures related to results section 1
# ---------------------------------------------------------------------------------------------------
rule figure1:
    input:
        ncoex_scatter_and_hist = rules.ncoex_scatter_and_hist.output.rds,
        pirna_coex_w_te_prop = rules.pirna_coex_w_te_prop.output.rds,
    output:
        pdf="results/figures2/figure1_main.overview.pdf"
    script:
        "../scripts/figures2/overview.main.R"

rule figure1_supp_01:
    input:
        rules.independent_dataset_correlation.output.rds,
        rules.mf_dataset_correlation.output.rds,
        "results/exploratory_and_descriptive/g_variance_explained.rds",
        config.get("MERGED_MODELS"),
        coex_dist = rules.coex_metric_distributions.output.rds,
        n_models = rules.n_models_per_filtering_step.output.rds,
    output:
        pdf="results/figures2/figure1_supp_01.extra-overview.pdf"
    script:
        "../scripts/figures2/extra-overiew.R"

# ---------------------------------------------------------------------------------------------------
# figures related to results section 2
# ---------------------------------------------------------------------------------------------------

rule figure2:
    input:
        male_gg_gsea = "results/enrichment/sig_main_male_max_abs_estimate_qnorm.gg_gsea.rds",
        female_gg_gsea = "results/enrichment/sig_main_female_max_abs_estimate_qnorm.gg_gsea.rds",
        coex_vs_seq = rules.coex_vs_seq_similarity.output.rds,
    output:
        pdf="results/figures2/figure2_main.tf-overrepresentation.pdf"
    script:
        "../scripts/figures2/tf-overrepresentation-in-coex.main.R"

# # ---------------------------------------------------------------------------------------------------
# # figures related to results section 3
# # ---------------------------------------------------------------------------------------------------

rule figure3:
    input:
        "results/deg/ourKD.de.grs.rds",
        "results/signatures/ourKD_gsea.rds",
    output:
        pdf = "results/figures2/figure3_main.knockdowns.pdf"
    script:
        "../scripts/figures2/knockdowns.main.R"

rule figure3_supp_01:
    input:
        "results/deg/s2rplus.res.tsv.gz",
        rules.pirna_enrichment_in_kd.output.de_pirna_fisher,
        rules.remap_peaks_dist_to_pirna.output.rds,
    output:
        pdf="results/figures2/figure3_supp_01.s2rplus.pdf"
    script:
        "../scripts/figures2/s2rplus-all-te-and-our-kd-factor-specific-kd-enrichment.R"


rule figure3_supp_02:
    input:
        "results/deg/ourKD.de.grs.rds",
        "results/resources/pirna_pathway.tsv"
    output:
        pdf="results/figures2/figure3_supp_02.pirna-volcano.pdf"
    script:
        "../scripts/figures2/pirna-and-tes-in-kd-volcano.R"

rule figure3_supp_03:
    input:
        "results/signatures/ourKD_gsea.rds",
        "upstream/final-models.collected-info.tsv.gz"
    output:
        pdf="results/figures2/figure3_supp_03.kd-vs-coex-direction-table.pdf"
    script:
        "../scripts/figures2/kd-vs-coex-direction-table.R"



# # ---------------------------------------------------------------------------------------------------
# # figures related to results section 4
# # ---------------------------------------------------------------------------------------------------

rule figure4:
    input:
        motif_comparison = "results/motifs/comparison/pan_denovo_comparison.meme.rds",
        motif_similarity = "results/motifs/comparison/pan_denovo_similarity.meme.rds",
    output:
       pdf="results/figures2/figure4_main.pan-motifs.pdf"
    script:
      "../scripts/figures2/pan-motifs.main.R"

rule figure4_supp_01:
    input:
        n_denovo_vs_sig_coef = "results/integrative/n_denovo_vs_sig_coef.pan.rds",
        pan_tree = "results/integrative/motif_and_coex_on_tree.pan.plot.rds",
    output:
        pdf="results/figures2/figure4_supp_01.pan-meme-denovo.pdf"
    script:
        "../scripts/figures2/pan-meme-denovo.R"

rule figure4_supp_02:
    input:
        motif_comparison = "results/motifs/comparison/pan_denovo_comparison.homer.rds",
        motif_similarity = "results/motifs/comparison/pan_denovo_similarity.homer.rds",
    output:
        pdf="results/figures2/figure4_supp_02.pan-homer-denovo.pdf"
    script:
        "../scripts/figures2/pan-homer-denovo.R"


rule figure4_supp_03:
    input:
        denovo_empirical_fdr = "results/motifs/streme_per_tf_empirical_fdr/pan_empirical_fdr.tsv",
        motif_comparison = "results/motifs/comparison/pan_denovo_comparison.streme.rds",
        motif_similarity = "results/motifs/comparison/pan_denovo_similarity.streme.rds",
    output:
        pdf="results/figures2/figure4_supp_03.pan-streme-denovo.pdf"
    script:
        "../scripts/figures2/pan-streme-denovo.R"


rule figure4_supp_04:
    """
    pan tracks and qc
    """
    output:
        pdf="results/figures2/figure4_supp_04.csem-tracks-and-qc-pan-profile.pdf",
    script:
        "../scripts/figures2/csem-tracks-and-qc-pan-profile.R"

rule figure4_supp_05:
    """
    validate csem peaks
    """
    input:
        "results/motifs/csem_peak_sea.pan.tsv.gz"
    output:
        pdf="results/figures2/figure4_supp_05.known-pan-motifs-in-csem-peaks.pdf",
    script:
        "../scripts/figures2/known-pan-motifs-in-csem-peaks.R"


rule figure4_supp_06:
    """
    h3k9 profiles and repetitiveness
    """
    input:
        repetitiveness = rules.chip_repetitiveness.output.rds,
    output:
        pdf ="results/figures2/figure4_supp_06.h3k9me3-profile-repetitiveness.pdf",
    script:
        "../scripts/figures2/h3k9me3-profile-repetitiveness.R"


rule figure4_supp_07:
    input:
        rules.csem_peaks_regioner.output,
    output:
        pdf="results/figures2/figure4_supp_07.csem-regioner.pdf",
    script:
        "../scripts/figures2/csem-regioner.R"

rule figure4_supp_08:
    input:
        "resources/putatively_bound_insertions.rds",
        "workflow/scripts/utils/plotting.R"
    output:
        pdf="results/figures2/figure4_supp_08.exemplary-bound-locus-1.pdf",
    script:
        "../scripts/figures2/exemplary-bound-locus-1.R"

rule figure4_supp_09:
    input:
        rules.calderon22_reanalysis.output,
    output:
        pdf="results/figures2/figure5_supp_09.calderon22.pdf",
    script:
        "../scripts/figures2/calderon22.R"

# # ---------------------------------------------------------------------------------------------------
# # figures related to results section 5
# # ---------------------------------------------------------------------------------------------------


rule figure5:
    input:
        crlg = rules.filter_phylo_and_plot_correlograms.output.rds,
        control_crlg = rules.plot_control_correlograms.output.rds,
        tre = "results/te_sequence_similarity/te_sketch_tidytree.rds",
    output:
        pdf="results/figures2/figure5_main.phylosignal.pdf"
    script:
        "../scripts/figures2/phylosignal.main.R"

rule figure5_supp_01:
    input:
        "results/te_sequence_similarity/te_sketch_tidytree.rds",
    output:
        pdf="results/figures2/figure5_supp_01.dashing2-tree.pdf"
    script:
        "../scripts/figures2/dashing2-te-phylo-supplement-01.R"


rule figure5_supp_02:
    input:
        
    output:
        pdf="results/figures2/figure5_supp_02.te-regulators-in-kd-and-remap.pdf"
    script:
        "../scripts/figures2/te-regulators-in-kd-and-remap-peak-prox.R"


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
        rules.figure3_supp_03.output,

        rules.figure4.output,
        rules.figure4_supp_01.output,
        rules.figure4_supp_02.output,
        rules.figure4_supp_03.output,
        rules.figure4_supp_04.output,
        rules.figure4_supp_05.output,
        rules.figure4_supp_06.output,
        rules.figure4_supp_07.output,
        rules.figure4_supp_08.output,
        rules.figure4_supp_09.output,

        rules.figure5.output,
        rules.figure5_supp_01.output,
        rules.figure5_supp_02.output,

