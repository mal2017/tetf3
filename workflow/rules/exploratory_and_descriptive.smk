rule n_models_per_filtering_step:
    input:
        mods = config.get("MERGED_MODELS"),
    output:
        rds = "results/exploratory_and_descriptive/n_models_per_filtering_step.rds",
    script:
        "../scripts/exploratory_and_descriptive/n_models_per_filtering_step.R"

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

rule ncoex_scatter_and_hist:
    input:
        mods = config.get("MERGED_MODELS"),
    output:
        rds = "results/exploratory_and_descriptive/ncoex_scatter_and_hist.rds",
    script:
        "../scripts/exploratory_and_descriptive/ncoex_scatter_and_hist.R"

rule pirna_genes_in_lms:
    input:
        mods = config.get("MERGED_MODELS"),
        pirna = rules.make_pirna_gene_list.output.tsv
    output:
        rds = "results/pirna/pirna_genes_in_lms.rds",
    script:
        "../scripts/exploratory_and_descriptive/piRNA_genes_in_lms.R"        


rule pirna_genes_heatmap:
    input:
        mods = config.get("MERGED_MODELS"),
        pirna = rules.make_pirna_gene_list.output.tsv
    output:
        rds = "results/pirna/pirna_genes_heatmap.rds",
    script:
        "../scripts/exploratory_and_descriptive/piRNA_genes_heatmap.R"        

rule pirna_coex_w_te_prop:
    input:
        rds =  rules.pirna_genes_in_lms.output.rds, # just gets 'valid' models and subsets columns
    output:
        rds = "results/pirna/pirna_coex_w_te_prop.rds",
    script:
        "../scripts/exploratory_and_descriptive/piRNA_coex_w_te_prop.R"

rule plot_variance_explained:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds="results/exploratory_and_descriptive/g_variance_explained.rds",
    script:
        "../scripts/exploratory_and_descriptive/plot_variance_explained.R"


rule exploratory_and_descriptive:
    input:
        rules.n_models_per_filtering_step.output,
        rules.coex_metric_distributions.output,
        rules.ncoex_scatter_and_hist.output,
        rules.pirna_genes_in_lms.output,
        rules.pirna_genes_heatmap.output,
        rules.pirna_coex_w_te_prop.output,
        rules.plot_variance_explained.output,