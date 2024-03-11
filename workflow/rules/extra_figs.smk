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


