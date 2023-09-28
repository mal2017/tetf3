rule independent_dataset_correlation:
    input:
        mods = config.get("MERGED_MODELS"),
        indep = config.get("INDEPENDENT_DATASET")
    output:
        rds = "results/coexpression_replication/intermediate/replicate_dataset_correlation.rds",
        json = "results/coexpression_replication/intermediate/replicate_dataset_correlation.stats.json"
    script:
        "../scripts/coexpression_replication/replicate_dataset_correlation.R"

rule mf_dataset_correlation:
    input:
        mods = config.get("MERGED_MODELS"),
    output:
        rds = "results/coexpression_replication/intermediate/mf_dataset_correlation.rds",
        json = "results/coexpression_replication/intermediate/mf_dataset_correlation.stats.json"
    script:
        "../scripts/coexpression_replication/mf_dataset_correlation.R"



rule coexpression_replication:
    input:
        rules.independent_dataset_correlation.output,
        rules.mf_dataset_correlation.output,