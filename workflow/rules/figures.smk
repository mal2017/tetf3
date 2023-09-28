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
        pan_tree = "results/integrative/motif_and_coex_on_tree.pan.plot.rds",
        repetitiveness = rules.chip_repetitiveness.output.rds,
        motif_comparison = "results/motifs/comparison/pan_denovo_comparison.rds",
        motif_similarity = "results/motifs/comparison/pan_denovo_similarity.rds",
        n_denovo_vs_sig_coef = "results/integrative/n_denovo_vs_sig_coef.pan.rds"
    output:
       "results/figures/figure5.pdf"
    script:
      "../scripts/figures/figure5.R"

# ---------------------------------------------------------------------------------------------------
# bring it together
# ---------------------------------------------------------------------------------------------------

rule figures:
    input:
        rules.figure1.output,
        rules.figure2.output,
        rules.figure3.output,
        rules.figure4.output,
        rules.figure5.output,
        rules.basic_exploratory_supplement_01.output,
        rules.calderon22_reanalysis_supplement.output,