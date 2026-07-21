// MMSEQS_FAMILY_CLUSTER — cluster the concatenated ingroup proteome into gene
// families (ADR-0002 Q2/Q3). Distinct from MMSEQS_CLUSTER (ADR-0001), which clusters
// the already-filtered *candidate* set; this clusters the whole ingroup up front to
// define the families that seed the profile search. Sensitive cascaded clustering
// (`-s`), coverage-guarded to avoid promiscuous-domain over-merging.
process MMSEQS_FAMILY_CLUSTER {
    label 'high_cpu'
    tag "family_cluster"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/families" }, mode: 'copy'

    input:
    path(ingroup_fa)

    output:
    path("families_cluster.tsv"),    emit: cluster_tsv
    path("families_rep_seq.fasta"),  emit: representatives

    shell:
    '''
    if [ ! -s !{ingroup_fa} ]; then
        # Empty ingroup FASTA — emit empty outputs rather than let mmseqs fail.
        touch families_cluster.tsv families_rep_seq.fasta
        exit 0
    fi

    tmpdir=${SCRATCH:-${TMPDIR:-/tmp}}
    mmseqs easy-cluster \
        !{ingroup_fa} \
        families \
        ${tmpdir}/mmseqs_fam_$$ \
        --threads !{task.cpus} \
        -s !{params.family_sensitivity} \
        --min-seq-id !{params.family_min_seq_id} \
        -c !{params.family_cov} \
        --cov-mode !{params.family_cov_mode} \
        --cluster-mode !{params.family_cluster_mode}
    '''
}
