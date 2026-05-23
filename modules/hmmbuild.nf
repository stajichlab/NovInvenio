process HMMBUILD {
    label 'low_cpu'
    tag "${cluster_id}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/hmms" }, mode: 'copy'

    input:
    tuple val(cluster_id), path(cluster_msa)   // MSA of a cluster (e.g. from muscle/mafft)

    output:
    tuple val(cluster_id), path("${cluster_id}.hmm")

    script:
    """
    hmmbuild --cpu ${task.cpus} ${cluster_id}.hmm ${cluster_msa}
    """
}
