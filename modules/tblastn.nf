process TBLASTN {
    label 'high_cpu'
    tag "${meta_rep.id}_vs_${meta_genome.id}"
    publishDir "${params.outdir}/${params.project}/tblastn", mode: 'copy'

    input:
    tuple val(meta_rep), path(rep_fa)         // cluster representative protein FASTA
    tuple val(meta_genome), path(genome_fa)   // outgroup genome (DNA)

    output:
    tuple val(meta_pair), path("${meta_rep.id}_vs_${meta_genome.id}.tblastn.tsv")

    script:
    meta_pair = [ query_id: meta_rep.id, target_id: meta_genome.id ]
    """
    makeblastdb -in ${genome_fa} -dbtype nucl -out genome_db 2>/dev/null

    tblastn \
        -query ${rep_fa} \
        -db genome_db \
        -outfmt "6 qseqid sseqid evalue bitscore pident length qstart qend sstart send" \
        -evalue ${params.evalue} \
        -num_threads ${task.cpus} \
        -out ${meta_rep.id}_vs_${meta_genome.id}.tblastn.tsv
    """
}
