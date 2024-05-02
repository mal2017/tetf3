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


rule te_sequence_similarity:
    input:
        rules.process_sketch_dist.output.dist,
        expand("results/te_sequence_similarity/te_{dist_type}_tidytree.rds",dist_type=["male-coex","female-coex","sketch"]),

