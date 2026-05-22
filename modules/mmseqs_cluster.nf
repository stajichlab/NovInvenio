process MMSEQS_CLUSTER {
    label 'high_cpu'
    tag "mmseqs_cluster"
    publishDir "${params.outdir}/${params.project}/clusters", mode: 'copy'

    input:
    path(candidates_fa)   // FASTA of candidate proteins

    output:
    path("clusters_rep_seq.fasta"),  emit: representatives
    path("clusters_all_seqs.fasta"), emit: all_seqs
    path("clusters_cluster.tsv"),    emit: tsv

    script:
    """
    mmseqs easy-cluster \
        ${candidates_fa} \
        clusters \
        tmp \
        --threads ${task.cpus} \
        --min-seq-id 0.3 \
        -c 0.8 \
        --cov-mode 0
    """
}
