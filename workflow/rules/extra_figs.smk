rule figure3_supp_01:
    """
    s2r+
    """"
    input:
        "results/deg/s2rplus.res.tsv.gz",
        rules.pirna_enrichment_in_kd.output.de_pirna_fisher,
        rules.make_gene_symbol_lookup.output,
        rules.s2rplus_coex_te_gsea_by_de.output,
    output:
        pdf="results/figures2/figure3_supp_01.s2rplus.pdf"
    script:
        "../scripts/figures2/s2rplus-all-te-and-our-kd-factor-specific-kd-enrichment.R"


rule figure4_supp_05:
    """
    validate csem peaks
    """
    input:
        "results/motifs/csem_peak_sea.pan.tsv.gz"
    output:
        pdf="results/figures2/figure4_supp_05.known-pan-motifs-in-csem-peaks.pdf",
    script:
        "../scripts/figures2/known-pan-motifs-in-csem-peaks.R"

rule figure4_supp_08:
    input:
        "resources/putatively_bound_insertions.rds",
        "workflow/scripts/utils/plotting.R"
    output:
        pdf="results/figures2/figure4_supp_08.exemplary-bound-locus-1.pdf",
    script:
        "../scripts/figures2/exemplary-bound-locus-1.R"


rule figure3_supp_03:
    input:
        "results/signatures/ourKD_gsea.rds",
        "upstream/final-models.collected-info.tsv.gz"
    output:
        pdf="results/figures2/figure3_supp_03.kd-vs-coex-direction-table.pdf"
    script:
        "../scripts/figures2/kd-vs-coex-direction-table.R"


rule figure5_supp_02:
    input:
        rules.encode_peaks_dist_to_pirna.output.rds,
    output:
        pdf="results/figures2/figure5_supp_02.te-regulators-in-kd-and-encode.pdf"
    script:
        "../scripts/figures2/te-regulators-in-kd-and-encode-peak-prox.R"


rule figure4_supp_02:
    """
    "upstream/csem_mosaics/bigwigs/",
    """
    input:
        repetitiveness = rules.chip_repetitiveness.output.rds,
    output:
        pdf ="results/figures2/figure4_supp_06.h3k9me3-profile.pdf",
    script:
        "../scripts/figures2/h3k9me3-profile-repetitiveness.R"



rule figure4_supp_05:
    """
    pan tracks and qc
    """
    input:
        "upstream/csem_mosaics/bigwigs/",
        rules.plot_quality_by_visual_pericent_inspection_status.output,
        rules.plot_quality_by_visual_pericent_inspection_status.output,
    output:
        pdf="results/figures/panels/figure4_supp_05.csem-tracks-and-qc-pan-profile.pdf",
    params:
        figtitle="Supplement 05 to Figure 4"
    script:
        "../scripts/figures/csem-tracks-and-qc-pan-profile.R"


rule figure5_supp_03:
    input:
        rules.calderon22_reanalysis.output,
    output:
        pdf="results/figures/panels/figure5_supp_03.calderon22.pdf",
    script:
        "../scripts/figures/calderon22.R"


rule s2rplus_coex_te_gsea_by_de:
    input:
        deg = rules.s2rplus_limma.output.tsv,
        coex = config.get("MERGED_MODELS")
    output:
        rds = "results/signatures/s2rplus_te_gsea.rds",
    script:
        "../scripts/signatures/tfrnai_gsea_de.R"


rule plot_motif_and_coex_on_tree:
    input:
        tidytree = "results/te_sequence_similarity/te_sketch_tidytree.rds",
        coex = config.get("MERGED_MODELS"),
        fimo = rules.fimo_denovo_motifs_tes.output.odir, 
    output:
        rds = "results/integrative/motif_and_coex_on_tree.{tf}.plot.rds",
    params:
        tf = "{tf}",
    script:
        "../scripts/integrative/plot_motif_and_coex_on_tree2.R"


rule calderon22_reanalysis:
    input:
        "upstream/calderon22_supercells.rds",
        "resources/Drosophila_melanogaster_TF.txt",
        "upstream/te_element_lookup.json",
        "upstream/final-models.collected-info.tsv.gz",
    output:
        df = "results/calderon22/calderon22_reanalysis_correlations.rds",
        supercell = "results/calderon22/calderon22_reanalysis_supercell.rds",
    script:
        "../scripts/scrna/calderon22.v2.R"

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