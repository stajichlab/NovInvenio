process PHMMER_SEARCH {
    label 'med_cpu'
    tag "${meta_q.id}_vs_${meta_t.id}"

    storeDir { "${params.outdir}/${Helpers.projectName(params)}/search_cache" }

    input:
    tuple val(meta_q), path(query_fa), val(meta_t), path(target_fa)

    output:
    tuple val(meta_pair), path("${meta_q.id}_vs_${meta_t.id}.phmmer.tblout.gz")

    script:
    meta_pair = [
        query_id:     meta_q.id,
        target_id:    meta_t.id,
        query_group:  meta_q.group,
        target_group: meta_t.group,
        tool:         'phmmer'
    ]
    def prefix = "${meta_q.id}_vs_${meta_t.id}"
    """
    phmmer \
        --cpu ${task.cpus} \
        --tblout ${prefix}.phmmer.tblout \
        --noali \
        ${query_fa} ${target_fa} \
        > /dev/null
    gzip ${prefix}.phmmer.tblout
    """
}
