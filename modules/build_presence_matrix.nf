process BUILD_PRESENCE_MATRIX {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(parsed_tsvs)         // collected list of all parsed hit TSVs
    path(paralog_cutoffs)     // collected list of per-species paralog_cutoffs.tsv files
    path(config_csv)
    val(query_group)          // 'IN' (default, novelty direction) or 'OUT' (loss direction)
    val(min_frac)             // presence threshold within query_group
    val(other_max_frac)       // max fraction of the other group a candidate may still be present in
    val(matrix_name)          // output matrix filename
    val(candidates_name)      // output candidates filename
    val(evalues_name)         // output e-value sidecar filename (report evidence only)
    val(targets_name)         // output target-protein-ID sidecar filename (report evidence only)

    output:
    path("${matrix_name}"),     emit: matrix
    path("${candidates_name}"), emit: candidates
    path("${evalues_name}"),    emit: evalues
    path("${targets_name}"),    emit: targets
    // Only written when params.other_coverage_floor_qcov is set (issue #158).
    path("${matrix_name.replaceAll(/\.tsv$/, '')}.coverage_floor_rejections.tsv"), emit: floor_rejections, optional: true

    script:
    def floor = params.other_coverage_floor_qcov
    def floor_args = floor ? "--other-coverage-floor-qcov ${floor} --output-coverage-floor-rejections ${matrix_name.replaceAll(/\.tsv$/, '')}.coverage_floor_rejections.tsv" : ''
    """
    build_presence_matrix.py \
        --hits ${parsed_tsvs} \
        --paralog-cutoffs ${paralog_cutoffs} \
        --config ${config_csv} \
        --ingroup-min-frac ${min_frac} \
        --query-group ${query_group} \
        --other-max-frac ${other_max_frac} \
        --paralog-competition-scope ${params.paralog_competition_scope} \
        --paralog-rescue-evalue ${params.paralog_rescue_evalue ?: 0} \
        ${params.paralog_rescue_delta != null ? "--paralog-rescue-delta ${params.paralog_rescue_delta}" : ''} \
        ${floor_args} \
        --output-matrix ${matrix_name} \
        --output-candidates ${candidates_name} \
        --output-evalues ${evalues_name} \
        --output-targets ${targets_name}
    """
}
