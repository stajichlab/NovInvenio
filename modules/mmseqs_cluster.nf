process MMSEQS_CLUSTER {
    label 'high_cpu'
    tag "mmseqs_cluster"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/clusters" }, mode: 'copy'

    input:
    path(candidates_fa)   // FASTA of candidate proteins

    output:
    path("clusters_rep_seq.fasta"),  emit: representatives
    path("clusters_all_seqs.fasta"), emit: all_seqs
    path("clusters_cluster.tsv"),    emit: tsv

    shell:
    '''
    if [ ! -s !{candidates_fa} ]; then
        # No candidates (e.g. a single-ingroup-species run where nothing
        # cleared the "absent from all outgroups" filter) — mmseqs
        # createdb fails ungracefully on an empty FASTA, so skip clustering
        # and emit empty outputs instead.
        touch clusters_rep_seq.fasta clusters_all_seqs.fasta clusters_cluster.tsv
        exit 0
    fi

    tmpdir=${SCRATCH:-${TMPDIR:-/tmp}}
    mmseqs easy-cluster \
        !{candidates_fa} \
        clusters \
        ${tmpdir}/mmseqs_$$ \
        --threads !{task.cpus} \
        --min-seq-id 0.3 \
        -c 0.8 \
        --cov-mode 0
    '''
}
