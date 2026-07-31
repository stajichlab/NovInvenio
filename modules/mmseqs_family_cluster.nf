// MMSEQS_FAMILY_CLUSTER — cluster the concatenated seed-group proteome into gene
// families (ADR-0002 Q2/Q3). The seed group is the ingroup for the novelty direction and
// the outgroup for the loss mirror (out_prefix distinguishes their outputs). Distinct from
// MMSEQS_CLUSTER (ADR-0001), which clusters the already-filtered *candidate* set; this
// clusters the whole seed group up front to define the families that seed the profile
// search. Sensitive cascaded clustering (`-s`), coverage-guarded against over-merging.
process MMSEQS_FAMILY_CLUSTER {
    label 'high_cpu'
    tag "family_cluster"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/${out_prefix}families" }, mode: 'copy'

    input:
    path(seed_fa)      // concatenated seed-group proteomes (ingroup for novelty, outgroup for loss)
    val(out_prefix)    // '' (novelty → families/) | 'loss_' (loss → loss_families/)

    output:
    path("families_cluster.tsv"),    emit: cluster_tsv
    path("families_rep_seq.fasta"),  emit: representatives

    shell:
    '''
    if [ ! -s !{seed_fa} ]; then
        # Empty seed FASTA — emit empty outputs rather than let mmseqs fail.
        touch families_cluster.tsv families_rep_seq.fasta
        exit 0
    fi

    tmpdir=${SCRATCH:-${TMPDIR:-/tmp}}
    mmseqs easy-cluster \
        !{seed_fa} \
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
