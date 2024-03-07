rule csem_peaks_regioner:
    """
    currently pulls directly from amarel, but will fix this later
    """
    input:
        ref_ins = config["REF_INS"],
        fixed_ref_ins = rules.putative_fixed_insertions.output.rds,
        coex_json = rules.coexpressed_tes_json.output.json,
    output:
        rds = "results/csem_mosaics/regioner.rds",
        segmentation = "results/csem_mosaics/regioner.segmentation.rds",
    script:
        "../scripts/csem_mosaics/regioner_v2.R"
