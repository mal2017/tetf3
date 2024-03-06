rule s2rplus_limma:
  input:
    se = config.get("S2RPLUS_RNAI_SE"),
    runselect = config.get("S2RPLUS_RNAI_RUNSELECTOR"),
    batch = config.get("S2RPLUS_RNAI_BATCH"),
    mods = config.get("MERGED_MODELS"),
  output:
    tsv = "results/deg/s2rplus.res.tsv.gz",
  script:
    "../scripts/differential_expression/basic_full_limma.R"


rule this_study_kd_deseq2:
  input:
    se = config.get("THIS_STUDY_RNAI"),
  output:
    grs = "results/deg/ourKD.de.grs.rds",
    dds = "results/deg/ourKD.dds.list.rds",
    df = "results/deg/ourKD.de.df.rds",
  script:
    "../scripts/differential_expression/ourKD_deseq2.v3.R"

rule plot_de_volcanos:
  input:
    grs = rules.this_study_kd_deseq2.output.grs,
    pirna = rules.make_pirna_gene_list.output.tsv,
  output:
    gg = "results/deg/de_volcanos.gg.rds",
  script:
    "../scripts/differential_expression/plot_de_volcanos.R"