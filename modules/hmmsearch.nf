process HMMSEARCH {
    label 'high_cpu'
    tag "${cluster_id}"
    publishDir "${params.outdir}/${params.project}/hmm_results", mode: 'copy'

    input:
    tuple val(cluster_id), path(hmm_file)
    path(target_db)   // SwissProt / UniProt FASTA

    output:
    tuple val(cluster_id), path("${cluster_id}.hmmsearch.tblout")

    script:
    """
    hmmsearch \
        --cpu ${task.cpus} \
        --tblout ${cluster_id}.hmmsearch.tblout \
        --noali \
        -E ${params.evalue} \
        ${hmm_file} ${target_db} \
        > /dev/null
    """
}
