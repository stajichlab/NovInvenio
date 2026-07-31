// Build a diamond database once per proteome (storeDir-cached) so that the
// pairwise search does not rebuild the same target DB for every query pairing.
process DIAMOND_MAKEDB {
    label 'high_cpu'
    tag "${meta.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    storeDir { "${params.outdir}/${Helpers.projectName(params)}/search_cache" }

    input:
    tuple val(meta), path(proteome_fa)

    output:
    tuple val(meta), path("${meta.id}.dmnd")

    script:
    """
    diamond makedb --in ${proteome_fa} --db ${meta.id} --quiet
    """
}

process DIAMOND_SEARCH {
    label 'high_cpu'
    tag "${meta_q.id}_vs_${meta_t.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    storeDir { "${params.outdir}/${Helpers.projectName(params)}/search_cache" }

    input:
    tuple val(meta_q), path(query_fa), val(meta_t), path(target_db)

    output:
    tuple val(meta_pair), path("${meta_q.id}_vs_${meta_t.id}.diamond.tsv.gz")

    script:
    meta_pair = [
        query_id:     meta_q.id,
        target_id:    meta_t.id,
        query_group:  meta_q.group,
        target_group: meta_t.group,
        tool:         'diamond'
    ]
    def prefix = "${meta_q.id}_vs_${meta_t.id}"
    """
    diamond blastp \
        --query ${query_fa} \
        --db ${target_db.baseName} \
        --outfmt 6 qseqid sseqid evalue bitscore \
        --evalue ${params.parse_evalue} \
        --threads ${task.cpus} \
        --quiet \
        --out ${prefix}.diamond.tsv
    gzip ${prefix}.diamond.tsv
    """
}
