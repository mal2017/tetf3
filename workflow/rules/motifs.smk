TFSOI=["pan"]


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
        4
    #singularity:
        #"docker://memesuite/memesuite:5.5.4"
    shell:
        """
        meme "{params.dir}/{wildcards.tf}/coex.fasta" \
            -neg "{params.dir}/{wildcards.tf}/other.fasta" \
            -oc '{output.odir}' \
            -nmotifs 15 -dna -mod anr -minw 6 -maxw 18 -p {threads} \
            -objfun se -revcomp
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
        4
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

rule get_known_motifs:
    """
    motifs id'd by downloading the meme docker image, starting a shell, and running the following shell command:
    `grep -E -e "FBgn0085432|pan" /opt/meme/share/meme-5.5.0/db/motif_databases/*/*.meme`

    'Archbold' refers to archbold 2014 plos genetics

    cisbp motif, flyreg, flyfactor, etc are all the same or very similar to jaspar1.1, so I don't include them here.
    i add the awk command to remove the space before the A in the nucleotide freqs record. this breaks universalmotif::read_meme.
    """
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
        meme2meme {output.archbold_hmg} {output.archbold_helper} /opt/meme/share/meme-5.5.4/db/motif_databases/JASPAR/JASPAR2022_CORE_insects_redundant_v2.meme | \
            awk '/^ A 0.25/{{sub(/^ /, "", $0)}}1' > {output.all_known}
        """

rule compare_motifs:
    input:
        meme = rules.streme_per_tf.output.odir,
        known_meme= rules.get_known_motifs.output.all_known,
    output:
        motif_comparison = "results/motifs/comparison/{tf}_denovo_comparison.rds",
        motif_similarity = "results/motifs/comparison/{tf}_denovo_similarity.rds",
    script:
        "../scripts/motifs/compare_motifs.R"

# -----------------------------------------------------------------------------
# Downstreama analysis
# -----------------------------------------------------------------------------



# checkpoint get_remap_peak_seqs:
#     input:
#         bed = rules.annotate_fixed_insertions.output.remap,
#         rpm = config.get("REF_INS"),
#         fa = config.get("GENOME_FA")
#     output:
#         odir = directory("results/motifs/remap_peaks/")
#     script:
#         "../scripts/motifs/get_remap_peak_seqs.R"

# rule sea_remap_peaks:
#     input:
#         dir = rules.get_remap_peak_seqs.output.odir,
#         #xstreme = rules.combine_xstreme_motifs.output.meme # to run on all motifs, comment this out and delete '/combined.meme' from sea command
#         xstreme = rules.xstreme_per_tf.output.odir,
#         #known_motifs = rules.get_pan_motifs.output.combined,
#     output:
#         odir = directory("results/motifs/sea_remap_peaks/{tf}")
#     singularity:
#         "docker://memesuite/memesuite:5.5.3"
#     shell:
#         """
#         sea -p '{input.dir}/{wildcards.tf}.fasta' -m '{input.xstreme}/combined.meme' --order 0 -oc '{output.odir}'
#         """

# def aggregate_sea(wildcards):
#     lms_checkpoint_output = checkpoints.split_cons_tes_per_tf.get(**wildcards).output.odir
#     remap_checkpoint_output = checkpoints.get_remap_peak_seqs.get(**wildcards).output.odir
#     #print(checkpoint_output)
#     wc_path = os.path.join(lms_checkpoint_output, "{tf}/coex.fasta")
#     tfs = glob_wildcards(wc_path).tf
#     filter_wc_path = os.path.join(remap_checkpoint_output, "{tf}.fasta")
#     filters = glob_wildcards(filter_wc_path).tf
#     filters = ["pan"]
#     return expand("results/motifs/sea_remap_peaks/{tf}", tf=[x for x in tfs if x in filters])


# rule collect_remap_peak_sea:
#     input:
#         seas = expand("results/motifs/sea_remap_peaks/{tf}", tf=["pan"])
#     output:
#         tsv = "results/motifs/remap_peak_sea.tsv.gz"
#     script:
#         "../scripts/motifs/collect_remap_peak_sea.R"

rule fimo_denovo_motifs_tes:
    input:
        tes = config.get("TE_FA"),
        streme = rules.streme_per_tf.output.odir,
    output:
        odir = directory("results/motifs/fimo_on_tes/denovo/{tf}"),
        tmp = temp("results/motifs/fimo_on_tes/denovo/tmp/all_tes_for_{tf}_fimo.fa")
    singularity:
        "docker://memesuite/memesuite:5.5.3"
    params:
        streme = rules.streme_per_tf.output.odir.replace("(","\(").replace(")","\)"),
        tmp = "results/motifs/fimo_on_tes/denovo/tmp/all_tes_for_{tf}_fimo.fa".replace("(","\(").replace(")","\)"),
    shell:
        """
        gunzip -c {input.tes} > "{params.tmp}" &&
        fimo --oc "{output.odir}" \
            --thresh 0.1  --qv-thresh \
            "{params.streme}/streme.txt" \
            "{params.tmp}"
        """

def aggregate_fimo_on_tes(wildcards):
    lms_checkpoint_output = checkpoints.split_cons_tes_per_tf.get(**wildcards).output.odir
    #print(checkpoint_output)
    wc_path = os.path.join(lms_checkpoint_output, "{tf}/coex.fasta")
    tfs = glob_wildcards(wc_path).tf

    return expand("results/motifs/fimo_on_tes/{tf}", tf=tfs)

rule motifs:
    input:
        expand("results/motifs/streme_per_tf/{tf}/", tf=TFSOI),
        expand("results/motifs/streme_per_tf_empirical_fdr/{tf}_empirical_fdr.tsv",tf=TFSOI),
        expand("results/motifs/meme_per_tf/{tf}/", tf=TFSOI),
        expand("results/motifs/homer_per_tf/{tf}/", tf=TFSOI),
        #rules.collect_remap_peak_sea.output,
        #expand("results/motifs/fimo_on_tes/denovo/{tf}", tf=TFSOI), #aggregate_fimo_on_tes,
        #expand("results/motifs/fimo_genome_wide/{tf}", tf=TFSOI),
        #expand("results/motifs/comparison/{tf}_denovo_comparison.rds", tf=TFSOI),
        #expand("results/motifs/fimo_on_tes/known/{known_motif_set}", known_motif_set=["combined_pan","all_known"]),
        