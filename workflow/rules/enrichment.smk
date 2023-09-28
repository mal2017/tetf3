rule go_gsea:
    input:
        tsv =  "results/rankings/{ranking}.tsv.gz",
    output:
        rds =  "results/enrichment/{ranking}.go_gsea.rds",
    script:
        "../scripts/enrichment/fb_go_gsea.R"

rule gg_gsea:
    input:
        tsv =  "results/rankings/{ranking}.tsv.gz",
        pirna = rules.make_pirna_gene_list.output.tsv,
        zad = rules.get_zad_genes.output.tsv,
        animaltfdb_tfs = config.get("TFS"),
    output:
        rds =  "results/enrichment/{ranking}.gg_gsea.rds",
    script:
        "../scripts/enrichment/fb_gg_gsea.R"


rnks = expand("{t}_estimate_qnorm",t=["mean_abs","max_abs","sum_abs","mean","max","sum","min"]) + expand("{t}_var_exp_x",t=["mean_signed","max_signed","mean","max","sum"])
enrichment_targets = expand("{s}_{d}_{m}_{r}",s=["sig","nofilt"], d = ["main","reps"],m=["male","female"], r=rnks)
#enrichment_targets = ["sig_main_female_max_abs_estimate_qnorm","sig_main_male_max_abs_estimate_qnorm", "nofilt_main_female_mean_abs_estimate_qnorm","nofilt_main_male_mean_abs_estimate_qnorm","sig_main_female_mean_abs_estimate_qnorm","sig_main_male_mean_abs_estimate_qnorm"]
rule enrichment:
    input:
        expand("results/enrichment/{r}.{gt}_gsea.rds",gt=["gg","go"],r=enrichment_targets)
