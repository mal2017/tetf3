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

rule plot_check_kds_by_chip_prox:
  input:
    res="results/deg/ourKD.de.df.rds",
    gr = "results/pirna/encode_peaks_dist_to_pirna.gr.rds",
  output:
    gg="results/deg/check_kds_by_chip_prox.gg.rds",
  script:
    "../scripts/differential_expression/plot_check_kds_by_chip_prox.R"