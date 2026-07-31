nextflow.enable.dsl=2

process PARSE_SELF_HITS {
    label 'low_cpu'
    tag "${meta.id}_self"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/self_hits" }, mode: 'copy'

    input:
    tuple val(meta), path(hits_file)

    output:
    tuple val(meta), path("${meta.id}.paralog_cutoffs.tsv"), emit: tsv

    script:
    """
    parse_self_hits.py \
        --input ${hits_file} \
        --format ${params.run_tool} \
        --proteome ${meta.id} \
        --output ${meta.id}.paralog_cutoffs.tsv
    """
}
