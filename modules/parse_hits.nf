process PARSE_HITS {
    label 'low_cpu'
    tag "${meta_pair.query_id}_vs_${meta_pair.target_id}"

    input:
    tuple val(meta_pair), path(hits_file)

    output:
    tuple val(meta_pair), path("${meta_pair.query_id}_vs_${meta_pair.target_id}.parsed.tsv")

    script:
    def out = "${meta_pair.query_id}_vs_${meta_pair.target_id}.parsed.tsv"
    """
    parse_hits.py \
        --input ${hits_file} \
        --format ${meta_pair.tool} \
        --evalue ${params.evalue} \
        --query-proteome ${meta_pair.query_id} \
        --target-proteome ${meta_pair.target_id} \
        --output ${out}
    """
}
