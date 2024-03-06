rule plot_n_models_per_filtering_step:
    input:
        mods = config.get("MERGED_MODELS"),
    output:
        gg = "results/exploratory_and_descriptive/n_models_per_filtering_step.gg.rds",
    script:
        "../scripts/exploratory_and_descriptive/plot_n_models_per_filtering_step.R"

rule coex_metric_distributions:
    input:
        mods = config.get("MERGED_MODELS"),
    output:
        rds = "results/exploratory_and_descriptive/coex_metric_distributions.rds",
    script:
        "../scripts/exploratory_and_descriptive/coex_metric_distributions.R"

# rule TE_expression_distribution:
#     input:
#         mods = config.get("MERGED_MODELS"),
#         se = config.get("THIS_STUDY_DGRP_RNA"),
#     output:
#         rds = "results/exploratory_and_descriptive/TE_expression_distribution.rds",
#     script:
#         "../scripts/exploratory_and_descriptive/TE_expression_distribution.R"

rule plot_n_features_scatter:
    input:
        mods = config.get("MERGED_MODELS"),
    output:
        n_tes_scatter = "results/exploratory_and_descriptive/n_tes_scatter.gg.rds",
        n_genes_scatter = "results/exploratory_and_descriptive/n_genes_scatter.gg.rds",
    script:
        "../scripts/exploratory_and_descriptive/plot_n_features_scatter.R"

rule plot_variance_explained:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds="results/exploratory_and_descriptive/variance_explained.gg.rds",
    script:
        "../scripts/exploratory_and_descriptive/plot_variance_explained.R"

rule plot_mf_shared_coex_hits:
    input:
        mods = config.get("MERGED_MODELS"),
    output:
        gg = "results/exploratory_and_descriptive/mf_shared_coex_hits.gg.rds",
    script:
        "../scripts/exploratory_and_descriptive/plot_mf_shared_coex_hits.R"