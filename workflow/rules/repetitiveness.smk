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

rule repetitiveness:
    input:
        rules.chip_repetitiveness.output.rds,



