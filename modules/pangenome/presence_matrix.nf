// PRESENCE_MATRIX — build the family x strain presence/absence matrix from
// the tier-1 cluster TSV. Producer of presence_matrix.tsv, consumed by the
// rescue pass, frequency binning, and co-occurrence steps.
process PRESENCE_MATRIX {
    label 'low_cpu'
    tag "presence_matrix"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(cluster_tsv)
    path(samplesheet)

    output:
    path("presence_matrix.tsv"), emit: matrix
    path("presence_matrix.copy_number.tsv"), optional: true, emit: copy_number

    script:
    """
    pangenome_build_presence_matrix.py \
        --cluster_tsv ${cluster_tsv} \
        --config ${samplesheet} \
        --groups '${params.pangenome_ingroup_label},${params.pangenome_outgroup_label}' \
        --id_sep '${params.pangenome_id_sep}' \
        --output presence_matrix.tsv
    """
}
