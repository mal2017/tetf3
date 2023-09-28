rule get_coex_distance:
    input:
        mods = config.get("MERGED_MODELS"),
        zads = rules.get_zad_genes.output.tsv,
        tfs = config.get("TFS"),
    output:
        rds = "results/te_sequence_similarity/coex_dist_df.rds",
        male_dist = "results/te_sequence_similarity/te_male-coex_dist.rds",
        female_dist = "results/te_sequence_similarity/te_female-coex_dist.rds",
    script:
        "../scripts/te_sequence_similarity/get_coex_distance_combined_sexes.R"


rule get_sketch_dist:
    """
    currently dashing v2.1.19
    emits similarity by default - if changing args to a different metric, make sure to change script process_sketch_dist.R accordingly
    """
    input:
        fa =  config.get("TE_FA"),
    output:
        sketch_dist = "results/te_sequence_similarity/te_sketch_dist.txt",
    params:
        args = config.get("DASHING_ARGS"),
        path = config.get("DASHING_PATH"),
    threads:
        10
    shell:
        """
        {params.path} dist {params.args} --parse-by-seq --square {input.fa} > {output.sketch_dist}
        """
        

rule process_sketch_dist: 
    """
    yields a dist object
    """
    input:
        mods = config.get("MERGED_MODELS"),
        txt = rules.get_sketch_dist.output.sketch_dist,
    output:
        dist = "results/te_sequence_similarity/te_sketch_dist.rds",
    params:
        dashing_type = config.get("DASHING_TYPE"),
    script:
        "../scripts/te_sequence_similarity/process_sketch_dist.R"


rule make_te_tree_from_dist:
    input:
        dist = "results/te_sequence_similarity/te_{dist_type}_dist.rds",
        te_classes = config.get("TE_CLASSES"),
    output:
        tidytree = "results/te_sequence_similarity/te_{dist_type}_tidytree.rds"
    script:
        "../scripts/te_sequence_similarity/make_te_tree_from_dist.R"


rule coex_vs_seq_similarity:
    """
    generate a tbl relating coexpression distance to sequence distance
    """
    input:
        coex_dist = rules.get_coex_distance.output.rds,
        seq_dist = rules.process_sketch_dist.output.dist,
        te_classes = config.get("TE_CLASSES"),
    output:
        rds = "results/te_sequence_similarity/coex_vs_seq_similarity.rds"
    script:
        "../scripts/te_sequence_similarity/coex_vs_seq_similarity.R"

rule te_sequence_similarity:
    input:
        rules.process_sketch_dist.output.dist,
        expand("results/te_sequence_similarity/te_{dist_type}_tidytree.rds",dist_type=["male-coex","female-coex","sketch"]),
        rules.coex_vs_seq_similarity.output.rds,

