rule unr_ripseq_enrichment:
    input:
        tsv = config.get("UNR_RIPSEQ_COUNTS"),
        gtf= config.get("COMBINED_GTF"),
        te_json="results/resources/coexpressed_tes.json",
    output:
        dds = "results/ripseq/unr_ripseq.dds.rds",
        tsv = "results/ripseq/unr_ripseq.tsv.gz"
    params:
        alpha = config.get("UNR_RIPSEQ_ALPHA")
    script:
        "../scripts/ripseq/unr_enrichment.R"

rule transcript_au_content:
    input:
        te_fa = config.get("TE_FA"),
        tx_fa = config.get("TX_FA"),
        te_json="results/resources/coexpressed_tes.json",
        rip_fl = "results/ripseq/unr_ripseq.tsv.gz"
    output:
        tsv = "results/ripseq/unr_ripseq_features_au_content.tsv.gz"
    script:
        "../scripts/ripseq/transcript_au_content.R"

rule transcript_au_content_in_region:
    input:
        te_fa = config.get("TE_FA"),
        tx_fa = config.get("TX_FA"),
        te_json="results/resources/coexpressed_tes.json",
        rip_fl = "results/ripseq/unr_ripseq.tsv.gz"
    params:
        relpos = config.get("UNR_RIPSEQ_TX_RELATIVE_POSITION")
    output:
        tsv = "results/ripseq/unr_ripseq_features_au_content_in_region.tsv.gz"
    script:
        "../scripts/ripseq/transcript_au_content_in_region.R"

rule transcript_attta_sites_in_region:
    input:
        te_fa = config.get("TE_FA"),
        tx_fa = config.get("TX_FA"),
        te_json="results/resources/coexpressed_tes.json",
        rip_fl = "results/ripseq/unr_ripseq.tsv.gz"
    params:
        relpos = config.get("UNR_RIPSEQ_TX_RELATIVE_POSITION")
    output:
        tsv = "results/ripseq/unr_ripseq_features_attta_sites_in_region.tsv.gz"
    script:
        "../scripts/ripseq/transcript_attta_sites_in_region.R"

rule unr_bound_tx_in_kd:
    input:
        lkup = "results/resources/gene_symbol_lookup.tsv.gz",
        kd = "results/deg/ourKD.de.df.rds",
        rip = "results/ripseq/unr_ripseq.tsv.gz",
    output:
        rds = "results/ripseq/unr_bound_tx_in_kd.gsea.rds"
    script:
        "../scripts/ripseq/unr_bound_tx_in_kd.R"

rule unr_ripseq_phylosignal:
    input:
        coex = config.get("MERGED_MODELS"),
        dist = "results/te_sequence_similarity/te_sketch_dist.rds",
        rip = "results/ripseq/unr_ripseq.tsv.gz",
        at_content = "results/ripseq/unr_ripseq_features_au_content.tsv.gz",
        attta = "results/ripseq/unr_ripseq_features_attta_sites.tsv.gz",
    params:
        relpos = config.get("UNR_RIPSEQ_TX_RELATIVE_POSITION")    
    output:
        tbl= "results/ripseq/unr_ripseq_phylosignal.tbl.rds",
        p4d = "results/ripseq/unr_ripseq_phylosignal.p4d.rds",
        tree = "results/ripseq/unr_ripseq_phylosignal.tree.rds",
    script:
        "../scripts/phylosignal/unr_ripseq_phylosignal.R"

rule get_ripseq_te_sequences:
    input:
        te_fa = config.get("TE_FA"),
        rip_fl = "results/ripseq/unr_ripseq.tsv.gz"
    output:
        fa = "results/ripseq/unr_ripseq_te_sequences.fa",
        non_bound_fa = "results/ripseq/unr_ripseq_non_bound_te_sequences.fa",
    script:
        "../scripts/ripseq/get_ripseq_te_sequences.R"



rule meme_unr_ripseq:
    """
    - -neg "{input.dir}/Unr/other.fasta" 
    """
    input:
        fa = "results/ripseq/unr_ripseq_te_sequences.fa",
        non_bound_fa = "results/ripseq/unr_ripseq_non_bound_te_sequences.fa",
    output:
        odir = directory("results/ripseq/unr_meme/")
    threads:
        6
    singularity:
        "docker://memesuite/memesuite:5.5.4"
    shell:
        """
        meme {input.fa} \
            -neg {input.non_bound_fa} \
            -oc '{output.odir}' \
            -nmotifs 10 -minw 5 -maxw 8 -dna -mod anr -p {threads} \
            -objfun se
        """



rule unr_ripseq_analysis:
    input:
        dds = "results/ripseq/unr_ripseq.dds.rds",
        tsv = "results/ripseq/unr_ripseq.tsv.gz",
        au = "results/ripseq/unr_ripseq_features_au_content.tsv.gz",
        au_in_region = "results/ripseq/unr_ripseq_features_au_content_in_region.tsv.gz",
        attta_in_region = "results/ripseq/unr_ripseq_features_attta_sites_in_region.tsv.gz",
        gsea = "results/ripseq/unr_bound_tx_in_kd.gsea.rds",
        p4d = "results/ripseq/unr_ripseq_phylosignal.p4d.rds",


