process BUILD_PRESENCE_MATRIX {
    label 'low_cpu'
    publishDir "${params.outdir}/${params.project}", mode: 'copy'

    input:
    path(parsed_tsvs)   // collected list of all parsed hit TSVs
    path(config_csv)

    output:
    path("presence_matrix.tsv"), emit: matrix
    path("candidates.txt"),      emit: candidates

    script:
    """
    build_presence_matrix.py \
        --hits ${parsed_tsvs} \
        --config ${config_csv} \
        --ingroup-min-frac ${params.ingroup_min_frac} \
        --output-matrix presence_matrix.tsv \
        --output-candidates candidates.txt
    """
}
