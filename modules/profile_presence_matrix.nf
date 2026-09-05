// PROFILE_PRESENCE_MATRIX — turn per-proteome family hmmsearch results + cluster
// membership into the standard presence_matrix.tsv + candidates.txt (ADR-0002).
// Emits the SAME contract as BUILD_PRESENCE_MATRIX so the downstream CLUSTER / VALIDATE
// / ANNOTATE / REPORT chain is reused unchanged.
process PROFILE_PRESENCE_MATRIX {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(domtblouts)
    path(cluster_tsv)
    path(families_tsv)
    path(protein_map)
    path(config_csv)
    val(query_group)      // 'IN' (novelty) | 'OUT' (loss)
    val(query_min_frac)   // presence threshold within the query group
    val(other_max_frac)   // max fraction of the other group a candidate may survive in
    val(out_prefix)       // '' (novelty) | 'loss_' — distinct output filenames per direction

    output:
    path("${out_prefix}presence_matrix.tsv"), emit: matrix
    path("${out_prefix}candidates.txt"),      emit: candidates

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
        --min-covered-residues ${params.hmm_presence_min_residues} \
        --ingroup-min-frac ${query_min_frac} \
        --query-group ${query_group} \
        --other-max-frac ${other_max_frac} \
        --output-matrix ${out_prefix}presence_matrix.tsv \
        --output-candidates ${out_prefix}candidates.txt
    """
}
