rule get_cg16779_testis_de_tes:
    input:
        fa = "results/motifs/bbmask_consensus_tes/consensus_tes.masked.fasta",
        res = "results/deg/ourKD.de.df.rds",
    output:
        fa = "results/motifs_cg16779/testis_de_tes.fa",
        nfa = "results/motifs_cg16779/testis_nonsig_tes.fa"
    script:
        "../scripts/motifs_cg16779/get_cg16779_de_tes_fa.R"

rule cg16779_testis_de_meme:
    input:
        fa = "results/motifs_cg16779/testis_de_tes.fa",
        nfa = "results/motifs_cg16779/testis_nonsig_tes.fa"
    output:
        odir = directory("results/motifs_cg16779/meme/")
    threads:
        6
    singularity:
        "docker://memesuite/memesuite:5.5.4"
    shell:
        """
        meme {input.fa} \
            -neg {input.nfa} \
            -oc '{output.odir}' \
            -nmotifs 15 -dna -mod anr -minw 6 -maxw 18 -p {threads} \
            -objfun se -revcomp || true
        """
