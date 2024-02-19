rule calderon22_reanalysis:
    input:
        "upstream/calderon22_supercells.rds",
        "resources/Drosophila_melanogaster_TF.txt",
        "upstream/te_element_lookup.json",
        "upstream/final-models.collected-info.tsv.gz",
    output:
        df = "results/calderon22/calderon22_reanalysis_correlations.rds",
    script:
        "../scripts/scrna/calderon22.v2.R"