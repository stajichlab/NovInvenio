// Build a blast protein database once per proteome (storeDir-cached) so that the
// pairwise search does not rebuild the same target DB for every query pairing.
process BLAST_MAKEDB {
    label 'high_cpu'
    tag "${meta.id}"

    storeDir { "${params.outdir}/${Helpers.projectName(params)}/search_cache" }

    input:
    tuple val(meta), path(proteome_fa)

    output:
    tuple val(meta), path("${meta.id}.blast_db.*")

    script:
    """
    makeblastdb -in ${proteome_fa} -dbtype prot -out ${meta.id}.blast_db -parse_seqids 2>/dev/null
    """
}

process BLAST_SEARCH {
    label 'high_cpu'
    tag "${meta_q.id}_vs_${meta_t.id}"

    storeDir { "${params.outdir}/${Helpers.projectName(params)}/search_cache" }

    input:
    tuple val(meta_q), path(query_fa), val(meta_t), path(target_db)

    output:
    tuple val(meta_pair), path("${meta_q.id}_vs_${meta_t.id}.blast.tsv.gz")

    script:
    meta_pair = [
        query_id:     meta_q.id,
        target_id:    meta_t.id,
        query_group:  meta_q.group,
        target_group: meta_t.group,
        tool:         'blast'
    ]
    def prefix = "${meta_q.id}_vs_${meta_t.id}"
    """
    blastp \
        -query ${query_fa} \
        -db ${meta_t.id}.blast_db \
        -outfmt "6 qseqid sseqid evalue bitscore" \
        -evalue ${params.parse_evalue} \
        -num_threads ${task.cpus} \
        -out ${prefix}.blast.tsv
    gzip ${prefix}.blast.tsv
    """
}
