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

    main:
    // Build each target genome DB once (storeDir-cached, keyed by meta.id so the
    // two directions never collide), then run one TBLASTN job per genome with
    // all reps searched together.
    genome_db_ch = TBLASTN_MAKEDB(genome_dna_ch)
    TBLASTN(genome_db_ch, representatives_fa)

    // Summarise TBLASTN hits — expand rep-level hits to cluster members
    SUMMARIZE_TBLASTN(
        TBLASTN.out.tsv.map { meta, tsv -> tsv }.collect(),
        cluster_tsv,
        summary_name
    )

    emit:
    tblastn_hits    = TBLASTN.out.tsv
    summary         = SUMMARIZE_TBLASTN.out.tsv
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
