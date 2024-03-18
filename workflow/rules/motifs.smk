TFSOI=["pan","NfI","CG16779","vvl"]


# -----------------------------------------------------------------------------
# Prep for motif analysis
# -----------------------------------------------------------------------------
rule mask_cons_tes:
    input:
        tes = config.get("TE_FA"),
    output:
        masked = "results/motifs/bbmask_consensus_tes/consensus_tes.masked.fasta",
    params:
        entropy = config.get("BBMASK_ENTROPY"),
        w = config.get("BBMASK_W"),
    singularity:
        "docker://quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0"
    shell:
        """
        bbmask.sh \
            in={input.tes} \
            out={output.masked} \
            w={params.w} \
            entropy={params.entropy}
        """

checkpoint split_cons_tes_per_tf:
    input:
        tfs = config.get("TFS"),
        mods = config.get("MERGED_MODELS"),
        tes = rules.mask_cons_tes.output.masked
    output:
        odir = directory("results/motifs/consensus_tes_per_tf/")
    script:
        "../scripts/motifs/split_consensus_tes_per_tf.R"

checkpoint split_cons_tes_per_tf_unmasked:
    input:
        tfs = config.get("TFS"),
        mods = config.get("MERGED_MODELS"),
        tes = config.get("TE_FA"),
    output:
        odir = directory("results/motifs/consensus_tes_per_tf_unmasked/")
    script:
        "../scripts/motifs/split_consensus_tes_per_tf.R"

rule shuffle_cons:
    input:
        dir = rules.split_cons_tes_per_tf.output.odir,
    output:
        p = "results/motifs/shuffle_cons/{tf}/coex.{shuf_rep}.fasta",
    params:
        fasta = rules.split_cons_tes_per_tf.output.odir + "/{tf}/coex.fasta",
        seed = "{shuf_rep}"
    script:
        "../scripts/motifs/create_shuffled_seqs.R"

rule mask_cons_tes_shuf:
    input:
        coex = rules.shuffle_cons.output.p,
    output:
        coex = "results/motifs/shuffle_cons/{tf}/coex.{shuf_rep}.masked.fasta",
    params:
        entropy = config.get("BBMASK_ENTROPY"),
        w = config.get("BBMASK_W"),
    singularity:
        "docker://quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0"
    shell:
        """
        bbmask.sh \
            in={input.coex} \
            out={output.coex} \
            w={params.w} \
            entropy={params.entropy}
        """

# -----------------------------------------------------------------------------
# MEME
# -----------------------------------------------------------------------------


rule meme_per_tf:
    """
    -
    """
    input:
        dir = rules.split_cons_tes_per_tf.output.odir,
    output:
        odir = directory("results/motifs/meme_per_tf/{tf}/")
    params:
        dir = rules.split_cons_tes_per_tf.output.odir.replace("(","\\(").replace(")","\\)"),
    threads:
        6
    #singularity:
        #"docker://memesuite/memesuite:5.5.4"
    shell:
        """
        meme "{params.dir}/{wildcards.tf}/coex.fasta" \
            -neg "{params.dir}/{wildcards.tf}/other.fasta" \
            -oc '{output.odir}' \
            -nmotifs 15 -dna -mod anr -minw 6 -maxw 18 -p {threads} \
            -objfun se -revcomp || true
        """


# -----------------------------------------------------------------------------
# STREME
# -----------------------------------------------------------------------------

rule streme_per_tf_shuffled:
    """
    """
    input:
        p = rules.mask_cons_tes_shuf.output.coex,
        dir = rules.split_cons_tes_per_tf.output.odir,
    output:
        odir = directory("results/motifs/streme_per_tf_shuffled/{tf}/{shuf_rep}")
    threads:
        1
    singularity:
        "docker://memesuite/memesuite:5.5.4"
    params: 
        args = config.get("STREME_ARGS")
    shell:
        """
        streme --oc '{output.odir}' \
            --p '{input.p}' \
            --n '{input.dir}/{wildcards.tf}/other.fasta' \
            --thresh 1 --nmotif 1 \
            {params.args}
        """

rule streme_per_tf:
    input:
        dir = rules.split_cons_tes_per_tf.output.odir,
    output:
        odir = directory("results/motifs/streme_per_tf/{tf}/")
    params:
        dir = rules.split_cons_tes_per_tf.output.odir.replace("(","\\(").replace(")","\\)"),
        args = config.get("STREME_ARGS")
    threads:
        1
    singularity:
        "docker://memesuite/memesuite:5.5.4"
    shell:
        """
        streme --oc '{output.odir}' \
            --p "{params.dir}/{wildcards.tf}/coex.fasta" \
            --n '{input.dir}/{wildcards.tf}/other.fasta' \
            {params.args}
        """

rule empirical_fdr_for_streme:
    input:
        streme = rules.streme_per_tf.output.odir,
        shuf = expand("results/motifs/streme_per_tf_shuffled/{{tf}}/{shuf_rep}",shuf_rep=list(range(1, 1001))),
    output:
        tsv = "results/motifs/streme_per_tf_empirical_fdr/{tf}_empirical_fdr.tsv"
    script:
        "../scripts/motifs/empirical_fdr.R"
        



# -----------------------------------------------------------------------------
# homer
# -----------------------------------------------------------------------------

rule homer_per_tf:
    """
    theres an environment file, but this still depends on homer being installed in
    the working environment
    """
    input:
        dir = rules.split_cons_tes_per_tf.output.odir,
    output:
        odir = directory("results/motifs/homer_per_tf/{tf}/")
    threads:
        2
    conda:
        "../envs/homer.yaml"
    shell:
        """
        findMotifs.pl "{input.dir}/{wildcards.tf}/coex.fasta" fasta \
            '{output.odir}' -fasta "{input.dir}/{wildcards.tf}/other.fasta" \
            -p {threads}
        """



# -----------------------------------------------------------------------------
# comparison to known motifs
# -----------------------------------------------------------------------------

rule archbold2motifs:
    """
    get representative motifs from archbold 
    (highly similar sequences clustered and made into meme files)
    """
    output:
        meme = "results/motifs/known_motifs/archbold_degenerate.meme",
        sils = "results/motifs/known_motifs/archbold_degenerate_clustering.rds",
        seqs = "results/motifs/known_motifs/archbold_degenerate_seqs.fasta",
    script:
        "../scripts/motifs/archbold_degenerate_2_motifs.R"

rule get_known_motifs:
    """
    motifs id'd by downloading the meme docker image, starting a shell, and running the following shell command:
    `grep -E -e "FBgn0085432|pan" /opt/meme/share/meme-5.5.0/db/motif_databases/*/*.meme`

    'Archbold' refers to archbold 2014 plos genetics

    cisbp motif, flyreg, flyfactor, etc are all the same or very similar to jaspar1.1, so I don't include them here.
    i add the awk command to remove the space before the A in the nucleotide freqs record. this breaks universalmotif::read_meme.
    """
    input:
        arch = rules.archbold2motifs.output.meme,
    output:
        jaspar = "results/motifs/known_motifs/jaspar_pan.meme",
        jaspar2 = "results/motifs/known_motifs/jaspar2_pan.meme",
        archbold_hmg = "results/motifs/known_motifs/archbold_hmg.meme",
        archbold_helper = "results/motifs/known_motifs/archbold_helper.meme",
        combined_pan = "results/motifs/known_motifs/combined_pan.meme",
        all_known = "results/motifs/known_motifs/all_known.meme"
    threads:
        1
    singularity:
        "docker://memesuite/memesuite:5.5.4"
    shell:
        """
        meme-get-motif -id MA0237.1 /opt/meme/share/meme-5.5.4/db/motif_databases/JASPAR/JASPAR2022_CORE_insects_redundant_v2.meme > {output.jaspar} &&
        meme-get-motif -id MA0237.2 /opt/meme/share/meme-5.5.4/db/motif_databases/JASPAR/JASPAR2022_CORE_insects_redundant_v2.meme > {output.jaspar2} &&
        iupac2meme -dna GCCGCCR > {output.archbold_helper} &&
        iupac2meme -dna SCTTTGWSW > {output.archbold_hmg} &&
        meme2meme {output.jaspar} {output.jaspar2} {output.archbold_hmg} {output.archbold_helper} | \
            awk '/^ A 0.25/{{sub(/^ /, "", $0)}}1' > {output.combined_pan}
        meme2meme {input.arch} {output.archbold_hmg} {output.archbold_helper} /opt/meme/share/meme-5.5.4/db/motif_databases/JASPAR/JASPAR2022_CORE_insects_redundant_v2.meme | \
            awk '/^ A 0.25/{{sub(/^ /, "", $0)}}1' > {output.all_known}
        """

rule compare_motifs:
    input:
        denovo = "results/motifs/{motif_program}_per_tf/{tf}/",
        archbold = rules.archbold2motifs.output.meme,
        known_meme= rules.get_known_motifs.output.all_known,
    params:
        motif_program = "{motif_program}",
    output:
        motif_comparison = "results/motifs/comparison/{tf}_denovo_comparison.{motif_program}.rds",
        motif_similarity = "results/motifs/comparison/{tf}_denovo_similarity.{motif_program}.rds",
    script:
        "../scripts/motifs/compare_motifs.R"

rule plot_motif_comparison:
    """
    in practice only used for pan, but adaptable to others
    by adding a wildcard
    """
    input:
        comparison = "results/motifs/comparison/pan_denovo_comparison.{motif_program}.rds",
    output:
        gg = "results/motifs/comparison/pan_denovo_comparison.{motif_program}.gg_df.rds",
    script:
        "../scripts/motifs/plot_motif_comparison.R"

rule compare_denovo_pan_motifs:
    input:
        homer = "results/motifs/homer_per_tf/pan/",
        streme = "results/motifs/streme_per_tf/pan/",
        meme = "results/motifs/meme_per_tf/pan/",
        known_meme = rules.get_known_motifs.output.all_known,
        comparison = "results/motifs/comparison/pan_denovo_comparison.meme.rds",
    output:
        motif_comparison = "results/motifs/comparison/pan_within_denovo_comparison.rds",
        motifs_um = "results/motifs/comparison/pan_within_denovo.universal_motif.rds",
        gg = "results/motifs/comparison/pan_within_denovo.gg.rds",
    script:
        "../scripts/motifs/compare_denovo_pan_motifs.R"


rule sea_known_motifs_on_tes:
    input:
        dir = "results/motifs/consensus_tes_per_tf_unmasked/",
        meme = rules.get_known_motifs.output.combined_pan,
        meme2 = rules.archbold2motifs.output.meme,
    output:
        odir = directory("results/motifs/sea_known_motifs_on_tes/pan")
    singularity:
        "docker://memesuite/memesuite:5.5.4"
    shell:
        """
        sea -p '{input.dir}/pan/coex.fasta' -n '{input.dir}/pan/other.fasta' -m '{input.meme}' -m '{input.meme2}' -oc '{output.odir}'
        """

csem_libraries, = glob_wildcards('upstream/csem_mosaics/sea/sea/pan_{lib}/sea.tsv')

# rule get_csem_pan_peaks:
#     input:
#         bed = "/home/mlawlor/amarel-matt/tetf/subworkflows/tetf_csem_mosaics/results/csem_mosaics/mosaics/pan_{library}/pan_{library}.mosaics.bed",
#         fa = config.get("GENOME_FA")
#     output:
#         fa = "results/motifs/csem_peaks/pan/{library}.fasta",
#         g = temp("results/motifs/csem_peaks/genome.tmp.{library}.fa"),
#     conda:
#         "../envs/bedtools.yaml"
#     shell:
#         """
#         gunzip -c {input.fa} > {output.g} &&
#         bedtools getfasta -fi {output.g} -bed {input.bed} -fo {output.fa}
#         """

# rule denovo_sea_csem_peaks:
#     input:
#         fa = "results/motifs/csem_peaks/pan/{library}.fasta",
#         meme = rules.meme_per_tf.output.odir,
#     output:
#         odir = directory("results/motifs/sea_csem_peaks/pan/{library}")
#     singularity:
#         "docker://memesuite/memesuite:5.5.4"
#     shell:
#         """
#         sea -p {input.fa} -m {input.meme}/meme.txt' --order 0 -oc '{output.odir}'
#         """

rule collect_csem_peak_sea:
    input:
        seas = expand("upstream/csem_mosaics/sea/sea/pan_{library}/", library=csem_libraries)
    output:
        tsv = "results/motifs/csem_peak_sea.known.pan.tsv.gz"
    script:
        "../scripts/motifs/collect_csem_peak_sea.R"

rule plot_upstream_csem_known_pan_sea:
    input:
        tsv = rules.collect_csem_peak_sea.output.tsv
    output:
        gg = "results/motifs/upstream_csem_known_pan_sea.gg.rds"
    script:
        "../scripts/motifs/plot_upstream_csem_known_pan_sea.R"

rule fimo_denovo_motifs_tes:
    input:
        tes = config.get("TE_FA"),
        meme = rules.meme_per_tf.output.odir,
    output:
        odir = directory("results/motifs/fimo_on_tes/denovo/{tf}"),
        tmp = temp("results/motifs/fimo_on_tes/denovo/tmp/all_tes_for_{tf}_fimo.fa")
    singularity:
        "docker://memesuite/memesuite:5.5.4"
    params:
        meme = rules.meme_per_tf.output.odir.replace("(","\(").replace(")","\)"),
        tmp = "results/motifs/fimo_on_tes/denovo/tmp/all_tes_for_{tf}_fimo.fa".replace("(","\(").replace(")","\)"),
    shell:
        """
        gunzip -c {input.tes} > "{params.tmp}" &&
        fimo --bfile --motif-- --oc "{output.odir}" \
            "{params.meme}/meme.txt" \
            "{params.tmp}"
        """

rule n_denovo_vs_sig_coef:
    input:
        coex = config.get("MERGED_MODELS"),
        fimo = rules.fimo_denovo_motifs_tes.output.odir,
    output:
        rds = "results/motifs/n_denovo_vs_sig_coef.{tf}.rds",
    params:
        tf= "{tf}"
    script:
        "../scripts/motifs/n_denovo_vs_sig_coef.R"

def aggregate_fimo_on_tes(wildcards):
    lms_checkpoint_output = checkpoints.split_cons_tes_per_tf.get(**wildcards).output.odir
    #print(checkpoint_output)
    wc_path = os.path.join(lms_checkpoint_output, "{tf}/coex.fasta")
    tfs = glob_wildcards(wc_path).tf

    return expand("results/motifs/fimo_on_tes/{tf}", tf=tfs)


rule motifs:
    input:
        expand("results/motifs/streme_per_tf_empirical_fdr/{tf}_empirical_fdr.tsv",tf=["pan"]),
        expand("results/motifs/meme_per_tf/{tf}/", tf=TFSOI),
        expand("results/motifs/comparison/{tf}_denovo_comparison.{p}.gg_df.rds", tf=["pan"],p=["meme","streme","homer"]),
        "results/motifs/n_denovo_vs_sig_coef.pan.rds",
        "results/motifs/csem_peak_sea.known.pan.tsv.gz",
        "results/motifs/upstream_csem_known_pan_sea.gg.rds",
        "results/motifs/sea_known_motifs_on_tes/pan",
        #"results/motifs/csem_peak_sea.pan.tsv.gz",
        #"results/motifs/sea_denovo_motifs_on_tes/pan/",
        #expand("results/motifs/fimo_on_tes/denovo/{tf}", tf=TFSOI), #aggregate_fimo_on_tes,
        