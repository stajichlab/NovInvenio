// SUMMARIZE_TBLASTN — build a protein x genome TBLASTN hit summary from per-genome TSVs,
// expanding cluster-representative hits out to every family member. Shared by
// novelty_discovery.nf (vs DISCOVERY_OUT genomes) and novelty_screen.nf (vs BROAD_OUTGROUP genomes) —
// both search the same family representative sequences (families_rep_seq.fasta), just
// against different genome panels, so this is the same summarization either way.
process SUMMARIZE_TBLASTN {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(tblastn_tsvs)
    path(cluster_tsv)
    val(summary_name)

    output:
    path("${summary_name}"), emit: tsv

    script:
    // Guard on tblastn_tsvs actually having entries: zero genomes is always a real
    // config error for VALIDATE's own outgroup_dna (pairwise/mmseqs direction), but a
    // legitimate, documented case for NOVELTY_SCREEN's BROAD_OUTGROUP (a config with no
    // BROAD_OUTGROUP rows degrades gracefully -- see workflows/novelty_screen.nf).
    // bin/summarize_tblastn.py's --hits is required/nargs='+', so an empty list
    // interpolates to "--hits " (no value) and argparse hard-fails before the script's
    // own "no genomes" check ever runs -- touch the (empty) output instead, same
    // "empty means no evidence" pattern ANNOTATE_MATRIX already uses for candidates_fa.
    """
    if [ -n "${tblastn_tsvs}" ]; then
        summarize_tblastn.py \
            --hits ${tblastn_tsvs} \
            --cluster_tsv ${cluster_tsv} \
            --evalue ${params.evalue} \
            --output ${summary_name}
    else
        touch ${summary_name}
    fi
    """
}
