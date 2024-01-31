rule make_gene_symbol_lookup:
    output:
        tsv = "results/resources/gene_symbol_lookup.tsv.gz"
    script:
        "../scripts/make_resources/gene_symbol_lookup.R"

rule make_pirna_gene_list:
    input:
        lkup = rules.make_gene_symbol_lookup.output.tsv,
        handler = "resources/handler2013_supp2.xlsx",
        czech = "resources/czech2013_supp2.xlsx",
    output:
        tsv = "results/resources/pirna_pathway.tsv"
    script:
        "../scripts/make_resources/get_piRNA_genes.R"

rule ref_preproc:
    """
    make a enome fasta with stripped names
    """
    params:
        genome_fa = config.get("GENOME_FA"),
    output:
        genome_fa = "results/resources/genome.fasta",
    script:
        "../scripts/make_resources/ref_preprocessing.R"

rule make_txdb:
    """
    make a reloadable txdb a
    """
    params:
        gtf = config.get("GTF"),
    output:
        txdb = "results/resources/txdb",
    script:
        "../scripts/make_resources/make_txdb.R"

rule get_zad_genes:
    output:
        tsv = "results/resources/zad_genes.tsv"
    script:
        "../scripts/make_resources/get_zad_genes.R"


rule te_sequences_split:
    """
    Split the TE insertions into their own file
    """
    input:
        fa = config.get("TE_FA")
    output:
        odir = directory("results/resources/te_sequences_split"),
    singularity:
        "docker://quay.io/biocontainers/seqkit:2.5.1--h9ee0642_0"
    shell:
        """
        seqkit split -i {input.fa} -O {output.odir}
        """

rule coexpressed_tes_json:
    """
    json with  top level elements all, male, female
    each of which has a separate entry for every gene_id and every gene_symbol
    holding the character vector of coexpressed tes for that gene
    """
    input:
        mods = config.get("MERGED_MODELS"),
    output:
        json = "results/resources/coexpressed_tes.json"
    script:
        "../scripts/make_resources/coexpressed_te_json.R"

rule putative_fixed_insertions:
    """
    find insertions fixed in all assayed DGRP lines and the reference
    """
    input:
        ref_ins = config.get("REF_INS"),
        uniq_tidal = config.get("PENETRANCE"),
    output:
        rds = "results/resources/putative_fixed_insertions.rds",
    script:
        "../scripts/make_resources/putative_fixed_insertions.R"

rule make_resources:
    """
    convenience rule to get all resources in one place
    """
    input:
        rules.make_gene_symbol_lookup.output.tsv,
        rules.make_pirna_gene_list.output.tsv,
        rules.ref_preproc.output.genome_fa,
        rules.make_txdb.output.txdb,
        rules.get_zad_genes.output.tsv,
        rules.coexpressed_tes_json.output.json,