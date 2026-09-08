// Build a nucleotide blast database once per outgroup genome (storeDir-cached)
// so the genome DB is not rebuilt on every pipeline run.
process TBLASTN_MAKEDB {
    label 'low_cpu'
    tag "${meta_genome.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    storeDir { "${params.outdir}/${Helpers.projectName(params)}/tblastn/db" }

    input:
    tuple val(meta_genome), path(genome_fa)

    output:
    tuple val(meta_genome), path("${meta_genome.id}.genome_db.*")

    script:
    """
    makeblastdb -in ${genome_fa} -dbtype nucl -out ${meta_genome.id}.genome_db 2>/dev/null
    """
}

process TBLASTN {
    label 'high_cpu'
    tag "${meta_genome.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/tblastn" }, mode: 'copy'

    input:
    tuple val(meta_genome), path(genome_db)    // pre-built outgroup genome DB
    path(reps_fa)                              // all cluster representative proteins

    output:
    tuple val(meta_genome), path("${meta_genome.id}.tblastn.tsv"), emit: tsv

    script:
    // Columns 0-9 (qseqid..send) are the original contract that
    // bin/summarize_tblastn.py's presence/absence matrix relies on -- new
    // columns are appended, never inserted, so any positional reader of the
    // first 10 fields is unaffected. sframe/qseq/sseq feed
    // bin/build_alignment_shards.py's pairwise-alignment archive: tblastn is
    // gapped by default, so qseq/sseq are always equal length (padded with
    // '-') and sframe disambiguates the minus-strand case where sstart > send.
    """
    tblastn \
        -query ${reps_fa} \
        -db ${meta_genome.id}.genome_db \
        -outfmt "6 qseqid sseqid evalue bitscore pident length qstart qend sstart send sframe qseq sseq" \
        -evalue ${params.evalue} \
        -num_threads ${task.cpus} \
        -out ${meta_genome.id}.tblastn.tsv
    """
}
