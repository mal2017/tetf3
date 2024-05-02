rule stats_descriptive_lms:
    input:
        mods = config.get("MERGED_MODELS"),
    output:
        json = "results/stats/descriptive_lms.json",
        xlsx = "results/tables/descriptive_lms.xlsx"
    script:
        "../scripts/tables/descriptive_lms.R"

rule lms_hit_enrichment:
    input:
        gsea_tbl_f = "results/enrichment/sig_main_female_max_abs_estimate_qnorm.gg_gsea.rds",
        gsea_tbl_m = "results/enrichment/sig_main_male_max_abs_estimate_qnorm.gg_gsea.rds"
    output:
        xlsx = "results/tables/lms_hit_enrichment.xlsx"
    script:
        "../scripts/tables/lms_hit_enrichment.R"

rule significant_pairs:
    input:
        config.get("MERGED_MODELS")
    output:
        xlsx = "results/tables/significant_pairs.xlsx",
    script:
        "../scripts/tables/significant_pairs.R"

rule metacell_table:
    input:
        supercell = "upstream/fca_supercells.rds",
        corr = "results/calderon22/fca_reanalysis_correlations.rds",
    output:
        xlsx= "results/tables/metacell_results.xlsx",
        sce = "results/tables/metacell_sce.rds",
    script:
        "../scripts/tables/metacell_table.R"


rule kds_te_enrichment:
    input:
        gsea_tbl = rules.ourKD_gsea.output.rds,
    output:
        xlsx = "results/tables/kds_te_enrichment.xlsx"
    script:
        "../scripts/tables/kds_te_enrichment.R"

rule motif_results_table:
    input:
        expand("results/motifs/streme_per_tf_empirical_fdr/{tf}_empirical_fdr.tsv",tf=TFSOI),
        expand("results/motifs/fimo_on_tes/denovo/{tf}",tf="pan"),
        "results/motifs/csem_peak_sea.known.pan.tsv.gz",
        "results/motifs/comparison/pan_denovo_comparison.meme.rds",
        "results/motifs/comparison/pan_denovo_comparison.homer.rds",
        "results/motifs/comparison/pan_denovo_comparison.streme.rds",
    output:
        xlsx= "results/tables/motif_results.xlsx",
    script:
        "../scripts/tables/motif_enrichment_results.R"

rule csem_mosaics_regioner_table:
    input:
        "results/csem_mosaics/regioner.rds",
    output:
       xlsx= "results/tables/csem_mosaics_regioner_table.xlsx",
    script:
        "../scripts/tables/csem_mosaics_regioner_table.R"

rule phylosignal_table:
    input:
        unr_phylosignal = "results/ripseq/unr_ripseq_phylosignal.tbl.rds",
        overall_phylosignal = "results/phylosignal/phylosignal_df.rds",
    output:
        xlsx= "results/tables/phylosignal_results.xlsx",
    script:
        "../scripts/tables/phylosignal_table.R"

rule ripseq_table:
    input:
        "results/ripseq/unr_ripseq_features_attta_sites_in_region.tsv.gz",
        "results/ripseq/unr_ripseq_features_au_content_in_region.tsv.gz",
        "results/ripseq/unr_ripseq.tsv.gz",
    output:
        xlsx= "results/tables/ripseq_table.xlsx",
    script:
        "../scripts/tables/ripseq_table.R"

# rule our_kd_stats:
#     input:
#         gg_pirna_in_kds = rules.plot_pirna_genes_in_our_kd_all.output.rds,
#     output:
#         json = "results/stats/our_kd_stats.json"
#     script:
#         "../scripts/stats/our_kd_stats.R"


# rule table_te_regulator_chip_prox:
#     input:
#         "results/resources/pirna_pathway.tsv",
#         "results/resources/gene_symbol_lookup.tsv.gz",
#         "results/pirna/encode_peaks_dist_to_pirna.gr.rds",
#     output:
#         docx= "results/tables/table_te_regulator_chip_prox.docx",
#     script:
#         "../scripts/tables/known_TE_regulators_in_chips_table.v2.R"

# rule table_te_regulator_in_kds:
#     input:
#         "results/resources/pirna_pathway.tsv",
#         "results/resources/gene_symbol_lookup.tsv.gz",
#         rules.this_study_kd_deseq2.output.grs,
#     output:
#         xlsx= "results/tables/table_known_te_regulators_in_knockdown.xlsx",
#     script:
#         "../scripts/tables/known_TE_regulators_in_knockdowns_table.R"






rule tables:
    input:
        "results/tables/descriptive_lms.xlsx",
        "results/tables/lms_hit_enrichment.xlsx",
        "results/tables/significant_pairs.xlsx",
        "results/tables/metacell_results.xlsx",
        "results/tables/kds_te_enrichment.xlsx",
        "results/tables/motif_results.xlsx",
        "results/tables/ripseq_table.xlsx",
        "results/tables/csem_mosaics_regioner_table.xlsx",
        "results/tables/phylosignal_results.xlsx",