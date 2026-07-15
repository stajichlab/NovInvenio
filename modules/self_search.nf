nextflow.enable.dsl=2

// Each process searches a proteome against itself requesting at most 2 hits
// with a relaxed e-value, so that the best non-self (rank-2) hit is captured
// even when it would normally fail the significance threshold.
// Results are stored in search_cache alongside regular pairwise hits.

process PHMMER_SELF {
    label 'med_cpu'
    tag "${meta.id}_self"
    storeDir { "${params.outdir}/${Helpers.projectName(params)}/search_cache" }

    input:
    tuple val(meta), path(proteome_fa)

    output:
    tuple val(meta), path("${meta.id}_vs_${meta.id}.phmmer.tblout.gz"), emit: hits

    script:
    def prefix = "${meta.id}_vs_${meta.id}"
    """
    phmmer \
        --cpu ${task.cpus} \
        --tblout ${prefix}.phmmer.tblout \
        --noali \
        -E 100 \
        ${proteome_fa} ${proteome_fa} \
        > /dev/null
    gzip ${prefix}.phmmer.tblout
    """
}

process DIAMOND_SELF {
    label 'high_cpu'
    tag "${meta.id}_self"
    storeDir { "${params.outdir}/${Helpers.projectName(params)}/search_cache" }

    input:
    tuple val(meta), path(proteome_fa)

    output:
    tuple val(meta), path("${meta.id}_vs_${meta.id}.diamond.tsv.gz"), emit: hits

    script:
    def prefix = "${meta.id}_vs_${meta.id}"
    """
    diamond makedb --in ${proteome_fa} --db self_db --quiet

    diamond blastp \
        --query ${proteome_fa} \
        --db self_db \
        --outfmt 6 qseqid sseqid evalue bitscore \
        --max-target-seqs 2 \
        --evalue 100 \
        --threads ${task.cpus} \
        --quiet \
        --out ${prefix}.diamond.tsv
    gzip ${prefix}.diamond.tsv
    """
}

process BLAST_SELF {
    label 'high_cpu'
    tag "${meta.id}_self"
    storeDir { "${params.outdir}/${Helpers.projectName(params)}/search_cache" }

    input:
    tuple val(meta), path(proteome_fa)

    output:
    tuple val(meta), path("${meta.id}_vs_${meta.id}.blast.tsv.gz"), emit: hits

    script:
    def prefix = "${meta.id}_vs_${meta.id}"
    """
    makeblastdb -in ${proteome_fa} -dbtype prot -out self_db -parse_seqids 2>/dev/null

    blastp \
        -query ${proteome_fa} \
        -db self_db \
        -outfmt "6 qseqid sseqid evalue bitscore" \
        -max_target_seqs 2 \
        -evalue 100 \
        -num_threads ${task.cpus} \
        -out ${prefix}.blast.tsv
    gzip ${prefix}.blast.tsv
    """
}
