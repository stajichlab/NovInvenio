nextflow.enable.dsl=2

// Each process searches a proteome against itself requesting at most 2 hits
// with a relaxed e-value, so that the best non-self (rank-2) hit is captured
// even when it would normally fail the significance threshold.
// Results are stored in search_cache alongside regular pairwise hits.
//
// DIAMOND_SELF runs at --very-sensitive (not diamond's default fast mode used
// by the pairwise DIAMOND_SEARCH/DIAMOND_SELF-equivalent cross-species search):
// confirmed on real data (pezizo_set1, N. crassa) that default/--sensitive/
// --more-sensitive modes miss the true HEX-1/eIF-5A within-genome paralog pair
// entirely (zero hits even at -e 100), so the paralog-competition filter had no
// paralog to compare against and a spurious HEX-1 "presence" call (from its
// distant cross-hit to an outgroup's eIF-5A) went uncaught. --very-sensitive
// finds the real pair (E=2.8e-06). Only 11 self-searches per study (one per
// proteome, not O(species^2)), so the extra sensitivity cost is negligible.

process PHMMER_SELF {
    label 'med_cpu'
    tag "${meta.id}_self"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
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
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
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
        --very-sensitive \
        --threads ${task.cpus} \
        --quiet \
        --out ${prefix}.diamond.tsv
    gzip ${prefix}.diamond.tsv
    """
}

process BLAST_SELF {
    label 'high_cpu'
    tag "${meta.id}_self"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
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
