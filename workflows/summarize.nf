nextflow.enable.dsl=2

// Produce per-species novelty candidate files (novelties.<SHORT>.tsv).
// Combines the annotated presence matrix with the TBLASTN summary: present in the ingroup
// and absent from all outgroup proteomes. TBLASTN hits against outgroup genomes are
// reported in the tblastn_outgroup_hits column but are NOT used as a filter (MAKE_NOVELTIES
// passes --skip_tblastn_filter), so a candidate can have outgroup genomic hits.

workflow SUMMARIZE {
    take:
    annotated_matrix   // path: presence_matrix.function.tsv
    tblastn_summary    // path: tblastn_summary.tsv (protein × genome hit matrix)
    cluster_tsv        // path: mmseqs *_cluster.tsv (rep -> member) for gene-family grouping
    config_csv         // path: analysis CSV
    other_evidence     // path: *.other_evidence.tsv.gz, or a 0-byte stub = not measured (issue #208)
    tblastn_coverage   // path: tblastn_summary.coverage.tsv.gz, or a 0-byte stub (issue #208)

    main:
    MAKE_NOVELTIES(annotated_matrix, tblastn_summary, cluster_tsv, config_csv,
                   other_evidence, tblastn_coverage)

    emit:
    novelties = MAKE_NOVELTIES.out.novelties
}

process MAKE_NOVELTIES {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(annotated_matrix)
    path(tblastn_summary)
    path(cluster_tsv)
    path(config_csv)
    path(other_evidence, stageAs: 'other_evidence.in.tsv.gz')
    path(tblastn_coverage, stageAs: 'tblastn_coverage.in.tsv.gz')

    output:
    path("novelties.*.tsv"), emit: novelties

    script:
    // 0-byte inputs are stubs ("not measured"): omit the flag so the columns stay absent.
    def ev_arg  = other_evidence.size() > 0 ? "--other_evidence ${other_evidence}" : ''
    def cov_arg = tblastn_coverage.size() > 0 ? "--tblastn_coverage ${tblastn_coverage}" : ''
    """
    make_novelties.py \
        --matrix ${annotated_matrix} \
        --tblastn_summary ${tblastn_summary} \
        --cluster_tsv ${cluster_tsv} \
        --config ${config_csv} \
        --ingroup_min ${params.ingroup_min_frac} \
        --other_signal_qcov ${params.other_signal_qcov} \
        ${ev_arg} ${cov_arg} \
        --output_dir . \
        --skip_tblastn_filter
    """
}
