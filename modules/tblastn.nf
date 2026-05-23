process TBLASTN {
    label 'high_cpu'
    tag "${meta_genome.id}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/tblastn" }, mode: 'copy'

    input:
    tuple val(meta_genome), path(genome_fa)   // outgroup genome (DNA)
    path(reps_fa)                              // all cluster representative proteins

    output:
    tuple val(meta_genome), path("${meta_genome.id}.tblastn.tsv"), emit: tsv

    script:
    """
    makeblastdb -in ${genome_fa} -dbtype nucl -out genome_db 2>/dev/null

    tblastn \
        -query ${reps_fa} \
        -db genome_db \
        -outfmt "6 qseqid sseqid evalue bitscore pident length qstart qend sstart send" \
        -evalue ${params.evalue} \
        -num_threads ${task.cpus} \
        -out ${meta_genome.id}.tblastn.tsv
    """
}
