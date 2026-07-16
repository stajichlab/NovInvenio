nextflow.enable.dsl=2

// Render the interactive HTML report (report.html) from the annotated presence
// matrix, the TBLASTN summary and the per-species novelties tables.
// The page is self-contained — no network access is needed to open it.

workflow REPORT {
    take:
    annotated_matrix   // path: presence_matrix.function.tsv
    tblastn_summary    // path: tblastn_summary.tsv
    novelties          // path(s): novelties.<SHORT>.tsv
    candidates_fa      // path: candidates.fa (protein sequences)
    config_csv         // path: analysis CSV

    main:
    MAKE_REPORT(annotated_matrix, tblastn_summary, novelties, candidates_fa, config_csv)

    emit:
    report = MAKE_REPORT.out.report
}

process MAKE_REPORT {
    label 'low_cpu'
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(annotated_matrix)
    path(tblastn_summary)
    path(novelties)
    path(candidates_fa)
    path(config_csv)

    output:
    path("report.html"), emit: report

    script:
    """
    make_report.py \
        --matrix ${annotated_matrix} \
        --config ${config_csv} \
        --tblastn_summary ${tblastn_summary} \
        --novelties ${novelties} \
        --candidates_fa ${candidates_fa} \
        --project ${Helpers.projectName(params)} \
        --ingroup_min_frac ${params.ingroup_min_frac} \
        --sequences ${params.report_sequences} \
        --output report.html
    """
}
