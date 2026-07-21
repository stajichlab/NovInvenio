// PROFILE_PRESENCE_MATRIX — turn per-proteome family hmmsearch results + cluster
// membership into the standard presence_matrix.tsv + candidates.txt (ADR-0002).
// Emits the SAME contract as BUILD_PRESENCE_MATRIX so the downstream CLUSTER / VALIDATE
// / ANNOTATE / REPORT chain is reused unchanged.
process PROFILE_PRESENCE_MATRIX {
    label 'low_cpu'
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(domtblouts)
    path(cluster_tsv)
    path(families_tsv)
    path(protein_map)
    path(config_csv)

    output:
    path("presence_matrix.tsv"), emit: matrix
    path("candidates.txt"),      emit: candidates

    script:
    """
    profile_to_matrix.py \
        --domtblout ${domtblouts} \
        --cluster-tsv ${cluster_tsv} \
        --families ${families_tsv} \
        --protein-map ${protein_map} \
        --config ${config_csv} \
        --evalue ${params.hmm_presence_evalue} \
        --min-coverage ${params.hmm_presence_cov} \
        --ingroup-min-frac ${params.ingroup_min_frac} \
        --query-group IN \
        --other-max-frac 0.0 \
        --output-matrix presence_matrix.tsv \
        --output-candidates candidates.txt
    """
}
