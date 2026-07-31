// EMPTY_EVALUES_STUB — placeholder e-value sidecar (issue #44) for presence-matrix
// producers that don't (yet) emit hit e-value evidence: the mmseqs/PROFILE_SEARCH family-
// profile pathway. lib/report_data.py's read_evalues() treats a missing/empty file as "no
// e-value evidence available" and the report simply omits that field — this keeps
// MAKE_REPORT's input contract uniform (always a real file) across every --cluster_tool
// without forcing every pathway to implement e-value tracking in this first pass.
process EMPTY_EVALUES_STUB {
    label 'low_cpu'
    tag "empty_evalues_stub"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    output:
    path("empty_evalues.tsv"), emit: evalues

    shell:
    '''
    : > empty_evalues.tsv
    '''
}
