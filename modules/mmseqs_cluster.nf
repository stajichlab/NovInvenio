process MMSEQS_CLUSTER {
    label 'high_cpu'
    tag "mmseqs_cluster"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/clusters" }, mode: 'copy'

    input:
    path(candidates_fa)   // FASTA of candidate proteins
    val(out_prefix)       // 'clusters' (novelty direction) or 'loss_clusters' (loss
                           // direction) — both publish into the same clusters/ subdir
                           // under distinct filenames, so a second invocation for the
                           // loss-search direction does not overwrite the first.

    output:
    path("${out_prefix}_rep_seq.fasta"),  emit: representatives
    path("${out_prefix}_all_seqs.fasta"), emit: all_seqs
    path("${out_prefix}_cluster.tsv"),    emit: tsv

    shell:
    '''
    if [ ! -s !{candidates_fa} ]; then
        # No candidates (e.g. a single-ingroup-species run where nothing
        # cleared the "absent from all outgroups" filter) — mmseqs
        # createdb fails ungracefully on an empty FASTA, so skip clustering
        # and emit empty outputs instead.
        touch !{out_prefix}_rep_seq.fasta !{out_prefix}_all_seqs.fasta !{out_prefix}_cluster.tsv
        exit 0
    fi

    tmpdir=${SCRATCH:-${TMPDIR:-/tmp}}
    mmseqs easy-cluster \
        !{candidates_fa} \
        !{out_prefix} \
        ${tmpdir}/mmseqs_$$ \
        --threads !{task.cpus} \
        --min-seq-id 0.3 \
        -c 0.8 \
        --cov-mode 0
    '''
}
