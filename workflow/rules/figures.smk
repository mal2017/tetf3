# ---------------------------------------------------------------------------------------------------
# figures related to results section 1
# ---------------------------------------------------------------------------------------------------
rule figure1:
    input:
        coex_dist = rules.coex_metric_distributions.output.rds,
        n_models = rules.n_models_per_filtering_step.output.rds,
        ncoex_scatter_and_hist = rules.ncoex_scatter_and_hist.output.rds,
        pirna_coex_w_te_prop = rules.pirna_coex_w_te_prop.output.rds,
    output:
        "results/figures/figure1.pdf"
    script:
        "../scripts/figures/figure1.R"

rule basic_exploratory_supplement_01:
    input:
        "results/coexpression_replication/intermediate/replicate_dataset_correlation.rds",
        "results/coexpression_replication/intermediate/mf_dataset_correlation.rds",
        "results/exploratory_and_descriptive/g_variance_explained.rds",
        config.get("MERGED_MODELS")
    output:
        pdf="results/figures/basic_exploratory_supplement_01.pdf"
    script:
        "../scripts/figures/basic-exploratory-supplement-01.R"

# ---------------------------------------------------------------------------------------------------
# figures related to results section 2
# ---------------------------------------------------------------------------------------------------

rule figure2:
    input:
        male_gg_gsea = "results/enrichment/sig_main_male_max_abs_estimate_qnorm.gg_gsea.rds",
        female_gg_gsea = "results/enrichment/sig_main_female_max_abs_estimate_qnorm.gg_gsea.rds",
        coex_vs_seq = rules.coex_vs_seq_similarity.output.rds,
    output:
        "results/figures/figure2.pdf"
    script:
        "../scripts/figures/figure2.R"

# ---------------------------------------------------------------------------------------------------
# figures related to results section 3
# ---------------------------------------------------------------------------------------------------

rule figure3:
    input:
        crlg = rules.plot_sig_correlograms.output.rds,
        control_crlg = rules.plot_control_correlograms.output.rds,
    output:
        "results/figures/figure3.pdf"
    script:
        "../scripts/figures/figure3.R"

# ---------------------------------------------------------------------------------------------------
# figures related to results section 4
# ---------------------------------------------------------------------------------------------------

rule figure4:
    input:
        "results/deg/ourKD.de.grs.rds",
        "results/signatures/ourKD_gsea.rds",
    output:
        "results/figures/figure4.pdf"
    script:
        "../scripts/figures/figure4.R"

rule s2rplus_supplement_01:
    input:
        "results/deg/s2rplus.res.tsv.gz"
    output:
        pdf="results/figures/s2rplus_supplement-01.pdf"
    script:
        "../scripts/figures/s2rplus_supplement_01.R"

rule calderon22_reanalysis_supplement:
    input:
        rules.calderon22_reanalysis.output,
    output:
        pdf="results/figures/calderon22_reanalysis_supplement-01.pdf",
        pdf2="results/figures/calderon22_reanalysis_supplement-02.pdf" 
    script:
        "../scripts/figures/calderon22_reanalysis_supplement.R"

# ---------------------------------------------------------------------------------------------------
# figures related to results section 5
# ---------------------------------------------------------------------------------------------------

rule figure5:
    input:
        motif_comparison = "results/motifs/comparison/pan_denovo_comparison.meme.rds",
        motif_similarity = "results/motifs/comparison/pan_denovo_similarity.meme.rds",
    output:
       "results/figures/figure5.pdf"
    script:
      "../scripts/figures/figure5_v2.R"

rule archbold_merging_supplement_01:
    input:
        sils = rules.archbold2motifs.output.sils,
        meme = rules.archbold2motifs.output.meme,
    output:
        pdf="results/figures/archbold_merging_supplement-01.pdf"
    script:
        "../scripts/figures/archbold_merging_supplement-01.R"

rule denovo_motifs_on_tes_supplement:
    input:
        n_denovo_vs_sig_coef = "results/integrative/n_denovo_vs_sig_coef.pan.rds",
        pan_tree = "results/integrative/motif_and_coex_on_tree.pan.plot.rds",
    output:
        pdf="results/figures/denovo-motifs-on-tes-supplement-01.pdf"
    script:
        "../scripts/figures/denovo-motifs-on-tes-supplement-01.R"


rule denovo_motifs_on_tes_supplement_homer:
    input:
        motif_comparison = "results/motifs/comparison/pan_denovo_comparison.homer.rds",
        motif_similarity = "results/motifs/comparison/pan_denovo_similarity.homer.rds",
    output:
        pdf="results/figures/denovo-motifs-on-tes-homer-supplement-01.pdf"
    script:
        "../scripts/figures/denovo-motifs-on-tes-homer-supplement-01.R"


rule denovo_motifs_on_tes_supplement_streme:
    input:
        denovo_empirical_fdr = "results/motifs/streme_per_tf_empirical_fdr/pan_empirical_fdr.tsv",
        motif_comparison = "results/motifs/comparison/pan_denovo_comparison.streme.rds",
        motif_similarity = "results/motifs/comparison/pan_denovo_similarity.streme.rds",
    output:
        pdf="results/figures/denovo-motifs-on-tes-streme-supplement-02.pdf"
    script:
        "../scripts/figures/denovo-motifs-on-tes-streme-supplement-02.R"


rule csem_tracks_and_quality_supplement_01:
    """
    pan tracks and qc
    """
    output:
        pdf="results/figures/csem-tracks-and-quality-supplement-pan-profile.pdf",
    script:
        "../scripts/figures/csem-tracks-and-quality-supplement-01.R"


rule csem_tracks_h3k9me3_profile_repetitiveness_supplement:
    """
    h3k9 profiles and repetitiveness
    """
    input:
        repetitiveness = rules.chip_repetitiveness.output.rds,
    output:
        pdf ="results/figures/csem-tracks-h3k9me3-profile-repetitiveness-supplement.pdf",
    script:
        "../scripts/figures/csem-tracks-and-quality-supplement-02.R"


rule csem_regioner_supplement_01:
    input:
        rules.csem_peaks_regioner.output,
    output:
        pdf="results/figures/csem-regioner-supplement-01.pdf",
    script:
        "../scripts/figures/csem-regioner-supplement-01.R"

rule exemplary_bound_loci_supplement:
    input:
        "resources/putatively_bound_insertions.rds",
        "workflow/scripts/utils/plotting.R"
    output:
        pdf1="results/figures/exemplary_bound_loci/exemplary_bound_locus_1.pdf",
    script:
        "../scripts/figures/exemplary_bound_loci_supplement-01.R"



# ---------------------------------------------------------------------------------------------------
# bring it together
# ---------------------------------------------------------------------------------------------------

rule figures:
    input:
        rules.figure1.output,
        rules.basic_exploratory_supplement_01.output,

        rules.figure2.output,
        rules.figure3.output,

        rules.figure4.output,
        rules.s2rplus_supplement_01.output,
        rules.calderon22_reanalysis_supplement.output,

        rules.figure5.output,
        rules.archbold_merging_supplement_01.output,
        rules.denovo_motifs_on_tes_supplement.output,
        rules.denovo_motifs_on_tes_supplement_homer.output,
        rules.denovo_motifs_on_tes_supplement_streme.output,
        #rules.csem_tracks_and_quality_supplement_01.output,
        #rules.csem_tracks_and_quality_supplement_02.output,
        rules.csem_tracks_h3k9me3_profile_repetitiveness_supplement.output,
        rules.csem_regioner_supplement_01.output,
        rules.exemplary_bound_loci_supplement.output,