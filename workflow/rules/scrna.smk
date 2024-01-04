rule calderon22_reanalysis:
    input:
        "upstream/calderon22_supercells.rds",
        "resources/Drosophila_melanogaster_TF.txt",
        "upstream/te_element_lookup.json",
        "upstream/final-models.collected-info.tsv.gz",
    output:
        df = "results/calderon22/calderon22_reanalysis_correlations.rds",
        gg_spqn = "results/calderon22/gg_spqn.rds",
        g_supercell_size = "results/calderon22/g_supercell_size.rds",
        g_all_cell_corr_pvals = "results/calderon22/g_all_cell_corr_pvals.rds",
        g_all_cell_overlapping_features_coexpressed = "results/calderon22/g_all_cell_overlapping_features_coexpressed.rds",
        g_pan_highly_corr_with_tes = "results/calderon22/g_pan_highly_corr_with_tes.rds",
        g_poscon_and_hits_panel = "results/calderon22/g_poscon_and_hits_panel.rds",
        g_negcon_panel = "results/calderon22/g_negcon_panel.rds",
        g_pan_vs = "results/calderon22/g_pan_vs.rds",
    script:
        "../scripts/scrna/calderon22.R"