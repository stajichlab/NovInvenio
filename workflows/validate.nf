nextflow.enable.dsl=2

include { TBLASTN_MAKEDB } from '../modules/tblastn'
include { TBLASTN        } from '../modules/tblastn'

workflow VALIDATE {
    take:
    representatives_fa   // path: all cluster representative proteins (single FASTA)
    genome_dna_ch        // [meta, genome_fa] — genome sequences to search against
                          //   (outgroup genomes for the novelty direction, ingroup
                          //   genomes for the loss direction)
    cluster_tsv           // path: mmseqs *_cluster.tsv (rep → member mapping)
    summary_name          // val: output filename, e.g. 'tblastn_summary.tsv' or
                           //   'loss_tblastn_summary.tsv'
    candidates_file        // path: candidates.txt or loss_candidates.txt — the actual
                           //   novelty/loss candidate list, used to filter which TBLASTN
                           //   hits get archived as alignment shards (issue #72)
    alignments_dir_name   // val: output directory name, e.g. 'alignments' or
                           //   'loss_alignments' (kept separate so the two directions'
                           //   docs/-published shards never collide)

    main:
    // Build each target genome DB once (storeDir-cached, keyed by meta.id so the
    // two directions never collide), then run one TBLASTN job per genome with
    // all reps searched together.
    genome_db_ch = TBLASTN_MAKEDB(genome_dna_ch)
    TBLASTN(genome_db_ch, representatives_fa)
    tblastn_tsv_ch = TBLASTN.out.tsv.map { meta, tsv -> tsv }.collect()

    // Summarise TBLASTN hits — expand rep-level hits to cluster members
    SUMMARIZE_TBLASTN(tblastn_tsv_ch, cluster_tsv, summary_name)

    // Per-genome alignment shards for the report's TBLASTN alignment popup
    // (docs/-only feature — see CLAUDE.md's report constraints and issue #72).
    BUILD_ALIGNMENT_SHARDS(tblastn_tsv_ch, candidates_file, cluster_tsv, alignments_dir_name)

    emit:
    tblastn_hits    = TBLASTN.out.tsv
    summary         = SUMMARIZE_TBLASTN.out.tsv
    alignments      = BUILD_ALIGNMENT_SHARDS.out.dir
}

process SUMMARIZE_TBLASTN {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(tblastn_tsvs)
    path(cluster_tsv)
    val(output_name)

    output:
    path("${output_name}"), emit: tsv

    script:
    """
    summarize_tblastn.py \
        --hits ${tblastn_tsvs} \
        --cluster_tsv ${cluster_tsv} \
        --evalue ${params.evalue} \
        --output ${output_name}
    """
}

// Per-genome gzip JSON alignment shards for the report's TBLASTN alignment
// popup (see bin/build_alignment_shards.py's module docstring). Published to
// docs/<project>/ (the Pages-served copy) rather than results/<project>/ --
// this data is a build artifact meant for GitHub-Release-asset publishing in
// the analysis repo, not a results/ deliverable and not meant to be
// git-committed (see NovInvenio's issue #73 discussion).
process BUILD_ALIGNMENT_SHARDS {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${Helpers.docsDir(params, workflow.launchDir)}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(tblastn_tsvs)
    path(candidates_file)
    path(cluster_tsv)
    val(outdir_name)

    output:
    path("${outdir_name}"), emit: dir

    script:
    """
    build_alignment_shards.py \
        --hits ${tblastn_tsvs} \
        --candidates ${candidates_file} \
        --cluster_tsv ${cluster_tsv} \
        --evalue ${params.evalue} \
        --project ${Helpers.projectName(params)} \
        --outdir ${outdir_name}
    """
}
