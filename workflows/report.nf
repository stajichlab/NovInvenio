nextflow.enable.dsl=2

// Render the interactive HTML reports:
//   - novelties.html — ingroup-specific candidates (annotated matrix + TBLASTN
//     summary + per-species novelties tables)
//   - core.html — near-universally conserved genes (annotated matrix +
//     cluster_tsv only, no new search/annotation step)
//   - losses.html — candidate lineage-specific gene losses (loss-search
//     annotated matrix + loss TBLASTN summary + loss cluster_tsv, all produced
//     by the LOSS_SEARCH/LOSS_CLUSTER/LOSS_VALIDATE/LOSS_ANNOTATE mirror
//     pipeline in main.nf)
// All three pages are self-contained — no network access is needed to open them.

workflow REPORT {
    take:
    annotated_matrix       // path: presence_matrix.function.tsv
    tblastn_summary        // path: tblastn_summary.tsv
    novelties               // path(s): novelties.<SHORT>.tsv
    candidates_fa           // path: candidates.fa (protein sequences)
    cluster_tsv              // path: mmseqs *_cluster.tsv (rep -> member) for gene-family grouping
    loss_annotated_matrix    // path: loss_presence_matrix.function.tsv
    loss_tblastn_summary     // path: loss_tblastn_summary.tsv
    loss_cluster_tsv         // path: loss mmseqs *_cluster.tsv
    config_csv               // path: analysis CSV

    main:
    MAKE_REPORT(annotated_matrix, tblastn_summary, novelties, candidates_fa, cluster_tsv, config_csv)
    MAKE_CORE_REPORT(annotated_matrix, cluster_tsv, config_csv)
    MAKE_LOSSES_REPORT(loss_annotated_matrix, loss_tblastn_summary, loss_cluster_tsv, config_csv)

    // Final step: gather the three reports under view/<project>/ with a report.html
    // landing page describing the run (ingroup/outgroup, tool, thresholds).
    COLLATE_REPORTS(
        MAKE_REPORT.out.report,
        MAKE_CORE_REPORT.out.report,
        MAKE_LOSSES_REPORT.out.report,
        config_csv,
    )

    emit:
    report        = MAKE_REPORT.out.report
    core_report   = MAKE_CORE_REPORT.out.report
    losses_report = MAKE_LOSSES_REPORT.out.report
    index_report  = COLLATE_REPORTS.out.index
}

process MAKE_REPORT {
    label 'low_cpu'
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(annotated_matrix)
    path(tblastn_summary)
    path(novelties)
    path(candidates_fa)
    path(cluster_tsv)
    path(config_csv)

    output:
    path("novelties.html"), emit: report

    script:
    """
    make_report.py \
        --matrix ${annotated_matrix} \
        --config ${config_csv} \
        --tblastn_summary ${tblastn_summary} \
        --novelties ${novelties} \
        --candidates_fa ${candidates_fa} \
        --cluster_tsv ${cluster_tsv} \
        --project ${Helpers.projectName(params)} \
        --ingroup_min_frac ${params.ingroup_min_frac} \
        --sequences ${params.report_sequences} \
        --output novelties.html
    """
}

process MAKE_CORE_REPORT {
    label 'low_cpu'
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(annotated_matrix)
    path(cluster_tsv)
    path(config_csv)

    output:
    path("core.html"), emit: report

    script:
    """
    make_core_report.py \
        --matrix ${annotated_matrix} \
        --config ${config_csv} \
        --cluster_tsv ${cluster_tsv} \
        --project ${Helpers.projectName(params)} \
        --core_min_frac ${params.core_min_frac} \
        --output core.html
    """
}

process MAKE_LOSSES_REPORT {
    label 'low_cpu'
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(loss_annotated_matrix)
    path(loss_tblastn_summary)
    path(loss_cluster_tsv)
    path(config_csv)

    output:
    path("losses.html"), emit: report

    script:
    """
    make_losses_report.py \
        --matrix ${loss_annotated_matrix} \
        --config ${config_csv} \
        --tblastn_summary ${loss_tblastn_summary} \
        --cluster_tsv ${loss_cluster_tsv} \
        --project ${Helpers.projectName(params)} \
        --outgroup_min_frac ${params.outgroup_min_frac} \
        --loss_ingroup_max_frac ${params.loss_ingroup_max_frac} \
        --output losses.html
    """
}

// Collate the three reports into view/<project>/ and generate a report.html
// landing page. The three HTML files are re-published here (copied) so the whole
// result set lives in one shareable folder alongside its index.
process COLLATE_REPORTS {
    label 'low_cpu'
    publishDir { "view/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(novelties_html)
    path(core_html)
    path(losses_html)
    path(config_csv)

    output:
    path("report.html"),     emit: index
    path(novelties_html),    emit: novelties
    path(core_html),         emit: core
    path(losses_html),       emit: losses

    script:
    """
    make_index_report.py \
        --config ${config_csv} \
        --project ${Helpers.projectName(params)} \
        --reports_dir . \
        --run_tool ${params.run_tool} \
        --ingroup_min_frac ${params.ingroup_min_frac} \
        --outgroup_min_frac ${params.outgroup_min_frac} \
        --loss_ingroup_max_frac ${params.loss_ingroup_max_frac} \
        --core_min_frac ${params.core_min_frac} \
        --output report.html
    """
}
