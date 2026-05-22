process BLAST_SEARCH {
    label 'high_cpu'
    tag "${meta_q.id}_vs_${meta_t.id}"

    storeDir "${params.outdir}/${params.project}/search_cache"

    input:
    tuple val(meta_q), path(query_fa), val(meta_t), path(target_fa)

    output:
    tuple val(meta_pair), path("${meta_q.id}_vs_${meta_t.id}.blast.tsv")

    script:
    meta_pair = [
        query_id:     meta_q.id,
        target_id:    meta_t.id,
        query_group:  meta_q.group,
        target_group: meta_t.group,
        tool:         'blast'
    ]
    """
    makeblastdb -in ${target_fa} -dbtype prot -out blast_db -parse_seqids 2>/dev/null

    blastp \
        -query ${query_fa} \
        -db blast_db \
        -outfmt "6 qseqid sseqid evalue bitscore" \
        -evalue ${params.evalue} \
        -num_threads ${task.cpus} \
        -out ${meta_q.id}_vs_${meta_t.id}.blast.tsv
    """
}
