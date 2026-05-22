process PHMMER_SEARCH {
    label 'high_cpu'
    tag "${meta_q.id}_vs_${meta_t.id}"

    storeDir "${params.outdir}/${params.project}/search_cache"

    input:
    tuple val(meta_q), path(query_fa), val(meta_t), path(target_fa)

    output:
    tuple val(meta_pair), path("${meta_q.id}_vs_${meta_t.id}.phmmer.tblout")

    script:
    meta_pair = [
        query_id:     meta_q.id,
        target_id:    meta_t.id,
        query_group:  meta_q.group,
        target_group: meta_t.group,
        tool:         'phmmer'
    ]
    """
    phmmer \
        --cpu ${task.cpus} \
        --tblout ${meta_q.id}_vs_${meta_t.id}.phmmer.tblout \
        --noali \
        ${query_fa} ${target_fa} \
        > /dev/null
    """
}
