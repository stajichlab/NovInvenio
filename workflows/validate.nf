nextflow.enable.dsl=2

include { TBLASTN_MAKEDB } from '../modules/tblastn'
include { TBLASTN        } from '../modules/tblastn'

workflow VALIDATE {
    take:
    representatives_fa   // path: all cluster representative proteins (single FASTA)
    outgroup_dna_ch      // [meta, genome_fa] — outgroup genome sequences
    cluster_tsv          // path: mmseqs *_cluster.tsv (rep → member mapping)

    main:
    // Build each outgroup genome DB once (storeDir-cached), then run one
    // TBLASTN job per genome with all reps searched together.
    genome_db_ch = TBLASTN_MAKEDB(outgroup_dna_ch)
    TBLASTN(genome_db_ch, representatives_fa)

    // Summarise TBLASTN hits — expand rep-level hits to cluster members
    SUMMARIZE_TBLASTN(
        TBLASTN.out.tsv.map { meta, tsv -> tsv }.collect(),
        cluster_tsv
    )

    emit:
    tblastn_hits    = TBLASTN.out.tsv
    summary         = SUMMARIZE_TBLASTN.out.tsv
}

process SUMMARIZE_TBLASTN {
    label 'low_cpu'
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(tblastn_tsvs)
    path(cluster_tsv)

    output:
    path("tblastn_summary.tsv"), emit: tsv

    script:
    """
    summarize_tblastn.py \
        --hits ${tblastn_tsvs} \
        --cluster_tsv ${cluster_tsv} \
        --evalue ${params.evalue} \
        --output tblastn_summary.tsv
    """
}
