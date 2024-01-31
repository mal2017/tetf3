# -----------------------------------------------------------------------------
# MEME
# -----------------------------------------------------------------------------


rule meme_unr:
    """
    - -neg "{input.dir}/Unr/other.fasta" 
    """
    input:
        dir = rules.split_cons_tes_per_tf.output.odir,
    output:
        odir = directory("results/motifs_unr/meme/")
    threads:
        6
    #singularity:
    #    "docker://memesuite/memesuite:5.5.4"
    shell:
        """
        meme "{input.dir}/Unr/coex.fasta" \
            -oc '{output.odir}' \
            -nmotifs 15 -minw 8 -maxw 18 -dna -mod anr -p {threads} \
            -objfun se
        """
