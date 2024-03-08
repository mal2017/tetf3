rule chip_qual_assessment:
    input:
        ss = config.get("ENCODE_SAMPLE_SHEET"),
        fingerprint = config.get("ENCODE_FINGERPRINT_DEEPTOOLS")
    output:
        rds = "results/repetitiveness/chip_qual_assessment.rds"
    script:
        "../scripts/repetitiveness/chip_qual_assessment.R"

rule chip_repetitiveness:
    input:
        qc = rules.chip_qual_assessment.output.rds,
        te_fa = config.get("TE_FA"),
        idxstats = config.get("ENCODE_MASKED_IDXSTATS"),
        ss = config.get("ENCODE_SAMPLE_SHEET"),
    output:
        rds = "results/repetitiveness/chip_repetitiveness.rds"
    script:
        "../scripts/repetitiveness/chip_repetitiveness.R"

rule plot_repetitiveness_by_visual_pericent_inspection_status:
    input:
        rds = rules.chip_repetitiveness.output.rds,
        json = "resources/pericent_enriched_pan_chips_by_inspection.json"
    output:
        gg = "results/repetitiveness/repetitiveness_by_visual_pericent_inspection_status.gg.rds"
    script:
        "../scripts/repetitiveness/plot_repetitiveness_by_visual_pericent_inspection_status.R"


rule plot_quality_by_visual_pericent_inspection_status:
    input:
        json = "resources/pericent_enriched_pan_chips_by_inspection.json",
        rep = rules.chip_repetitiveness.output.rds,
        qc = rules.chip_qual_assessment.output.rds
    output:
        gg = "results/repetitiveness/quality_by_visual_pericent_inspection_status.gg.rds"
    script:
        "../scripts/repetitiveness/plot_quality_by_visual_pericent_inspection_status.R"