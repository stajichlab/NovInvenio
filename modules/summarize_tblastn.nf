// SUMMARIZE_TBLASTN — build a protein x genome TBLASTN hit summary from per-genome TSVs,
// expanding cluster-representative hits out to every family member. Shared by
// novelty_discovery.nf (vs DISC_OUT genomes) and novelty_screen.nf (vs BROAD_OUT genomes) —
// both search the same family representative sequences (families_rep_seq.fasta), just
// against different genome panels, so this is the same summarization either way.
process SUMMARIZE_TBLASTN {
    label 'low_cpu'
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(tblastn_tsvs)
    path(cluster_tsv)
    val(summary_name)

    output:
    path("${summary_name}"), emit: tsv

    script:
    """
    summarize_tblastn.py \
        --hits ${tblastn_tsvs} \
        --cluster_tsv ${cluster_tsv} \
        --evalue ${params.evalue} \
        --output ${summary_name}
    """
}
