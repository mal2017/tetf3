rule get_peptides:
    params:
        fasta = "http://ftp.flybase.net/releases/FB2021_04/dmel_r6.41/fasta/dmel-all-translation-r6.41.fasta.gz"
    output:
        "results/peptides/all.fasta.gz"
    shell:
        "wget -O {output} {params.fasta}"

rule longest_peptide:
    input:
        fasta=rules.get_peptides.output
    output:
        fasta = "results/peptides/longest.fasta"
    script:
        "../scripts/zad_phylo/longest_peptide.R"

rule get_hmm:
    params:
        hmm = "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{pf}?annotation=hmm"
    output:
        "results/zad/hmm/{pf}.hmm.gz"
    shell:
        "wget -O {output} {params.hmm}"

rule get_domain_hits:
    """
    http://cryptogenomicon.org/extracting-hmmer-results-to-sequence-files-easel-miniapplications.html
    """
    input:
        fasta=rules.longest_peptide.output,
        hmm=rules.get_hmm.output
    output:
        sto = temp("results/zad/hmmer/domain-hits/{pf}.sto"),
        fasta = "results/zad/hmmer/domain-hits/{pf}.fasta",
    singularity:
        "docker://quay.io/biocontainers/hmmer:3.3.2--h87f3376_2"
    shell:
        """
        hmmsearch -A {output.sto} {input.hmm} {input.fasta} &&
        esl-reformat fasta {output.sto} > {output.fasta}
        """

rule filter_domain_hits:
    """
    for now, just takes the n most N-terminal hits
    """
    input:
        fasta=rules.get_domain_hits.output.fasta
    output:
        fasta = "results/zad/hmmer/domain-hits/{pf}.filtered.fasta"
    params:
        n = lambda wc: config.get("ALIGN_N_HITS").get(wc.pf)
    script:
        "../scripts/zad_phylo/filter_domain_hits.R"


rule align_domain_hits:
    """
    equiv to linsi
    """
    input:
        fasta=rules.filter_domain_hits.output.fasta
    output:
        "results/zad/mafft/domain-hits/{pf}.afa"
    singularity:
        "docker://quay.io/biocontainers/mafft:7.515--hec16e2b_0"
    threads:
        8
    shell:
        "mafft --thread {threads} --maxiterate 1000 --localpair {input.fasta} > {output}"

rule clipkit_domain_hits:
    input:
        afa=rules.align_domain_hits.output
    output:
        afa = "results/zad/mafft/domain-hits/{pf}.clipkit.afa"
    singularity:
        "docker://quay.io/biocontainers/clipkit:1.3.0--pyhdfd78af_0"
    params:
        trimmer = lambda wc: config.get("CLIPKIT_MODEL").get(wc.pf)
    shell:
        "clipkit {input.afa} -l -o {output.afa} -m {params.trimmer}"

rule iqtree_domain_hits:
    input:
        afa=rules.clipkit_domain_hits.output.afa,
    output:
        directory("results/zad/iqtree/domain-hits/{pf}")
    singularity:
        "docker://quay.io/biocontainers/iqtree:2.2.0.3--hb97b32f_1"
    params:
        model = lambda wc: config.get("IQTREE_MODEL").get(wc.pf)
    threads:
        2
    log:
        "results/zad/logs/iqtree_domain_hits/{pf}.log"
    shell:
        """
        mkdir -p {output} &&
        iqtree -T {threads} -B 1000 -s {input.afa} --wbtl \
            -m {params.model} --prefix {output}/{wildcards.pf} 2> {log}
        """

rule zad_phylo:
    input:
        "results/peptides/longest.fasta",
        "results/zad/hmm/PF07776.hmm.gz",
        "results/zad/hmmer/domain-hits/PF07776.fasta",
        "results/zad/mafft/domain-hits/PF07776.afa",
        expand("results/zad/iqtree/domain-hits/{pf}/",pf=["PF07776"])