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

# # -----------------------------------------------------------------------------
# # STREME
# # -----------------------------------------------------------------------------

# rule streme_per_tf_shuffled:
#     """
#     """
#     input:
#         p = rules.mask_cons_tes_shuf.output.coex,
#         dir = rules.split_cons_tes_per_tf.output.odir,
#     output:
#         odir = directory("results/motifs/streme_per_tf_shuffled/{tf}/{shuf_rep}")
#     threads:
#         1
#     singularity:
#         "docker://memesuite/memesuite:5.5.4"
#     params: 
#         args = config.get("STREME_ARGS")
#     shell:
#         """
#         streme --oc '{output.odir}' \
#             --p '{input.p}' \
#             --n '{input.dir}/{wildcards.tf}/other.fasta' \
#             --thresh 1 --nmotif 1 \
#             {params.args}
#         """

# rule streme_per_tf:
#     input:
#         dir = rules.split_cons_tes_per_tf.output.odir,
#     output:
#         odir = directory("results/motifs/streme_per_tf/{tf}/")
#     params:
#         dir = rules.split_cons_tes_per_tf.output.odir.replace("(","\\(").replace(")","\\)"),
#         args = config.get("STREME_ARGS")
#     threads:
#         1
#     singularity:
#         "docker://memesuite/memesuite:5.5.4"
#     shell:
#         """
#         streme --oc '{output.odir}' \
#             --p "{params.dir}/{wildcards.tf}/coex.fasta" \
#             --n '{input.dir}/{wildcards.tf}/other.fasta' \
#             {params.args}
#         """

# rule empirical_fdr_for_streme:
#     input:
#         streme = rules.streme_per_tf.output.odir,
#         shuf = expand("results/motifs/streme_per_tf_shuffled/{{tf}}/{shuf_rep}",shuf_rep=list(range(1, 1001))),
#     output:
#         tsv = "results/motifs/streme_per_tf_empirical_fdr/{tf}_empirical_fdr.tsv"
#     script:
#         "../scripts/motifs/empirical_fdr.R"
        



# # -----------------------------------------------------------------------------
# # homer
# # -----------------------------------------------------------------------------

# rule homer_per_tf:
#     """
#     theres an environment file, but this still depends on homer being installed in
#     the working environment
#     """
#     input:
#         dir = rules.split_cons_tes_per_tf.output.odir,
#     output:
#         odir = directory("results/motifs/homer_per_tf/{tf}/")
#     threads:
#         4
#     conda:
#         "../envs/homer.yaml"
#     shell:
#         """
#         findMotifs.pl "{input.dir}/{wildcards.tf}/coex.fasta" fasta \
#             '{output.odir}' -fasta "{input.dir}/{wildcards.tf}/other.fasta" \
#             -p {threads}
#         """



# # -----------------------------------------------------------------------------
# # comparison to known motifs
# # -----------------------------------------------------------------------------

# rule archbold2motifs:
#     output:
#         meme = "results/motifs/known_motifs/archbold_degenerate.meme",
#         sils = "results/motifs/known_motifs/archbold_degenerate_clustering.rds",
#         seqs = "results/motifs/known_motifs/archbold_degenerate_seqs.fasta",
#     script:
#         "../scripts/motifs/archbold_degenerate_2_motifs.R"

# rule get_known_motifs:
#     """
#     motifs id'd by downloading the meme docker image, starting a shell, and running the following shell command:
#     `grep -E -e "FBgn0085432|pan" /opt/meme/share/meme-5.5.0/db/motif_databases/*/*.meme`

#     'Archbold' refers to archbold 2014 plos genetics

#     cisbp motif, flyreg, flyfactor, etc are all the same or very similar to jaspar1.1, so I don't include them here.
#     i add the awk command to remove the space before the A in the nucleotide freqs record. this breaks universalmotif::read_meme.
#     """
#     input:
#         arch = rules.archbold2motifs.output.meme,
#     output:
#         jaspar = "results/motifs/known_motifs/jaspar_pan.meme",
#         jaspar2 = "results/motifs/known_motifs/jaspar2_pan.meme",
#         archbold_hmg = "results/motifs/known_motifs/archbold_hmg.meme",
#         archbold_helper = "results/motifs/known_motifs/archbold_helper.meme",
#         combined_pan = "results/motifs/known_motifs/combined_pan.meme",
#         all_known = "results/motifs/known_motifs/all_known.meme"
#     threads:
#         1
#     singularity:
#         "docker://memesuite/memesuite:5.5.4"
#     shell:
#         """
#         meme-get-motif -id MA0237.1 /opt/meme/share/meme-5.5.4/db/motif_databases/JASPAR/JASPAR2022_CORE_insects_redundant_v2.meme > {output.jaspar} &&
#         meme-get-motif -id MA0237.2 /opt/meme/share/meme-5.5.4/db/motif_databases/JASPAR/JASPAR2022_CORE_insects_redundant_v2.meme > {output.jaspar2} &&
#         iupac2meme -dna GCCGCCR > {output.archbold_helper} &&
#         iupac2meme -dna SCTTTGWSW > {output.archbold_hmg} &&
#         meme2meme {output.jaspar} {output.jaspar2} {output.archbold_hmg} {output.archbold_helper} | \
#             awk '/^ A 0.25/{{sub(/^ /, "", $0)}}1' > {output.combined_pan}
#         meme2meme {input.arch} {output.archbold_hmg} {output.archbold_helper} /opt/meme/share/meme-5.5.4/db/motif_databases/JASPAR/JASPAR2022_CORE_insects_redundant_v2.meme | \
#             awk '/^ A 0.25/{{sub(/^ /, "", $0)}}1' > {output.all_known}
#         """

# rule compare_motifs:
#     input:
#         denovo = "results/motifs/{motif_program}_per_tf/{tf}/",
#         archbold = rules.archbold2motifs.output.meme,
#         known_meme= rules.get_known_motifs.output.all_known,
#     params:
#         motif_program = "{motif_program}",
#     output:
#         motif_comparison = "results/motifs/comparison/{tf}_denovo_comparison.{motif_program}.rds",
#         motif_similarity = "results/motifs/comparison/{tf}_denovo_similarity.{motif_program}.rds",
#     script:
#         "../scripts/motifs/compare_motifs.R"

# rule sea_denovo_motifs_on_tes:
#     input:
#         dir = rules.split_cons_tes_per_tf.output.odir,
#         meme = rules.meme_per_tf.output.odir,
#     output:
#         odir = directory("results/motifs/sea_denovo_motifs_on_tes/{tf}")
#     singularity:
#         "docker://memesuite/memesuite:5.5.3"
#     shell:
#         """
#         sea -p '{input.dir}/{wildcards.tf}/coex.fasta' -n '{input.dir}/{wildcards.tf}/other.fasta' -m '{input.meme}/meme.txt' -oc '{output.odir}'
#         """

# rule sea_known_motifs_on_tes:
#     input:
#         dir = "results/motifs/consensus_tes_per_tf/pan/",
#         meme = rules.get_known_motifs.output.combined_pan,
#         meme2 = rules.archbold2motifs.output.meme,
#     output:
#         odir = directory("results/motifs/sea_known_motifs_on_tes/pan")
#     singularity:
#         "docker://memesuite/memesuite:5.5.3"
#     shell:
#         """
#         sea -p '{input.dir}/coex.fasta' -n '{input.dir}/other.fasta' -m '{input.meme}' -m '{input.meme2}' -oc '{output.odir}'
#         """

# csem_libraries,zz = glob_wildcards('/home/mlawlor/amarel-matt/tetf/subworkflows/tetf_csem_mosaics/results/csem_mosaics/mosaics/pan_{lib}/{lib2}.bed')

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

# rule sea_csem_peaks:
#     input:
#         fa = "results/motifs/csem_peaks/pan/{library}.fasta",
#         meme = rules.get_known_motifs.output.jaspar,
#         meme2 = rules.get_known_motifs.output.jaspar2,
#     output:
#         odir = directory("results/motifs/sea_csem_peaks/pan/{library}")
#     singularity:
#         "docker://memesuite/memesuite:5.5.3"
#     shell:
#         """
#         sea -p {input.fa} -m '{input.meme}' -m '{input.meme2}' --order 0 -oc '{output.odir}'
#         """


# rule collect_csem_peak_sea:
#     input:
#         seas = expand("results/motifs/sea_csem_peaks/pan/{library}/", library=csem_libraries)
#     output:
#         tsv = "results/motifs/csem_peak_sea.pan.tsv.gz"
#     script:
#         "../scripts/motifs/collect_csem_peak_sea.R"


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
#         meme = rules.get_known_motifs.output.jaspar,
#         meme2 = rules.get_known_motifs.output.jaspar2,
#     output:
#         odir = directory("results/motifs/sea_remap_peaks/pan")
#     singularity:
#         "docker://memesuite/memesuite:5.5.3"
#     shell:
#         """
#         sea -p '{input.dir}/pan.fasta' -m '{input.meme}' -m '{input.meme2}' --order 0 -oc '{output.odir}'
#         """


rule motifs_unr:
    input:
        "results/motifs_unr/meme",