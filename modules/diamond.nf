// Build a diamond database once per proteome (storeDir-cached) so that the
// pairwise search does not rebuild the same target DB for every query pairing.
process DIAMOND_MAKEDB {
    label 'high_cpu'
    tag "${meta.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    storeDir { "${params.outdir}/${Helpers.projectName(params)}/search_cache" }

    input:
    tuple val(meta), path(proteome_fa)

    output:
    tuple val(meta), path("${meta.id}.dmnd")

    script:
    """
    diamond makedb --in ${proteome_fa} --db ${meta.id} --quiet
    """
}

// One job per QUERY genome, searching every target in a loop inside the task,
// instead of one job per (query, target) PAIR -- for an N-ingroup x M-other-genome
// study this cuts job count from N*M to N. Real per-pair diamond compute is tiny
// (~1-3s even for this study's largest proteomes, measured directly) so the job
// COUNT was the real cost driver, not diamond's own runtime -- the same "trade job
// count for job size" pattern nextflow.config's hmm_search_chunk_size comment
// already documents for FAMILY_HMMSEARCH, just not previously applied here.
//
// Deliberately NOT storeDir (unlike DIAMOND_MAKEDB above): storeDir's skip-if-
// exists check is defined against a single named output, not a variable-length
// glob -- with a batched job's output being "however many targets this query had,"
// Nextflow can't reliably verify a partial cache hit before running. Using normal
// resume-based task caching + publishDir instead: still skips unchanged work on a
// `-resume` within the same run, and as a bonus, this is no longer invisible to
// nextflow's own trace files the way storeDir was (see the cluster-vs-pairwise-
// sensitivity investigation's Analysis 3, which found Tier P's true search cost
// unmeasurable for exactly this reason).
process DIAMOND_SEARCH {
    label 'high_cpu'
    tag "${meta_q.id}_vs_all"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    publishDir { "${params.outdir}/${Helpers.projectName(params)}/search_cache" }, mode: 'copy'

    input:
    tuple val(meta_q), path(query_fa), val(meta_t_list), path(target_db_list)

    output:
    path("*.diamond.tsv.gz")

    script:
    def searches = [meta_t_list, target_db_list].transpose().collect { meta_t, target_db ->
        def prefix = "${meta_q.id}_vs_${meta_t.id}"
        """
        diamond blastp \\
            --query ${query_fa} \\
            --db ${target_db.baseName} \\
            --outfmt 6 qseqid sseqid evalue bitscore length pident qcovhsp scovhsp qlen slen \\
            --evalue ${params.parse_evalue} \\
            --threads ${task.cpus} \\
            ${params.diamond_sensitivity ? "--${params.diamond_sensitivity}" : ''} \\
            --quiet \\
            --out ${prefix}.diamond.tsv
        gzip ${prefix}.diamond.tsv
        """
    }.join('\n')
    """
    ${searches}
    """
}
