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

rule pan_unr_summed_te_counts:
    input:
        sce = config.get('FCA_METACELLS_SCE'),
        json = rules.coexpressed_tes_json.output.json,
        ripseq = rules.unr_ripseq_enrichment.output.tsv,
    output:
        tsv = "results/calderon22/pan_unr_summed_te_counts.rds",
    script:
        "../scripts/scrna/pan_unr_summed_te_counts.R"