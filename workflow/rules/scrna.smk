rule fca_heads_reanalysis:
    input:
        "upstream/fca_supercells.rds",
        "resources/Drosophila_melanogaster_TF.txt",
        "upstream/te_element_lookup.json",
        "upstream/final-models.collected-info.tsv.gz",
    output:
        df = "results/calderon22/fca_reanalysis_correlations.rds",
        supercell = "results/calderon22/fca_reanalysis_supercell.rds",
    script:
        "../scripts/scrna/fca_supercell.R"