// PROFILE_CANDIDATE_CLUSTERS — family-as-cluster (ADR-0002 Q7). Reuse the profile
// pathway's existing gene families (restricted to candidate-containing families) as the
// VALIDATE/ANNOTATE/REPORT cluster unit, instead of re-running mmseqs on the candidates
// (the generic CLUSTER workflow). Emits the SAME three outputs CLUSTER does, so main.nf
// can wire either interchangeably.
process PROFILE_CANDIDATE_CLUSTERS {
    label 'low_cpu'
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(candidates_txt)
    path(family_cluster_tsv)   // families_cluster.tsv (rep<TAB>member)
    path(family_reps_fa)       // families_rep_seq.fasta
    path(ingroup_fa)           // concatenated ingroup proteomes
    val(candidates_fa_name)    // 'candidates.fa' (matches the pairwise pathway)

    output:
    path("${candidates_fa_name}"),          emit: candidates_fa
    path("candidate_family_reps.fasta"),    emit: representatives
    path("candidate_families_cluster.tsv"), emit: cluster_tsv

    script:
    """
    candidate_families.py \
        --candidates ${candidates_txt} \
        --family-cluster-tsv ${family_cluster_tsv} \
        --family-reps ${family_reps_fa} \
        --ingroup-fasta ${ingroup_fa} \
        --out-cluster-tsv candidate_families_cluster.tsv \
        --out-representatives candidate_family_reps.fasta \
        --out-candidates-fa ${candidates_fa_name}
    """
}
