process DIAMOND_SEARCH {
    label 'high_cpu'
    tag "${meta_q.id}_vs_${meta_t.id}"

    storeDir "${params.outdir}/${params.project}/search_cache"

    input:
    tuple val(meta_q), path(query_fa), val(meta_t), path(target_fa)

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
    diamond makedb --in ${target_fa} --db target_db --quiet

    diamond blastp \
        --query ${query_fa} \
        --db target_db \
        --outfmt 6 qseqid sseqid evalue bitscore \
        --evalue ${params.evalue} \
        --threads ${task.cpus} \
        --quiet \
        --out ${prefix}.diamond.tsv
    gzip ${prefix}.diamond.tsv
    """
}
