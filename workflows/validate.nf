nextflow.enable.dsl=2

include { TBLASTN } from '../modules/tblastn'

workflow VALIDATE {
    take:
    representatives_fa   // path: cluster representative proteins
    outgroup_dna_ch      // [meta, genome_fa] — outgroup genome sequences

    main:
    // Split representative FASTA into one entry per representative for parallel TBLASTN
    rep_ch = SPLIT_REPS(representatives_fa)
        .flatten()
        .map { fa -> [ [id: fa.baseName], fa ] }

    // Cross every representative against every outgroup genome
    pairs_ch = rep_ch.combine(outgroup_dna_ch)
        .map { meta_rep, rep_fa, meta_g, genome_fa -> [ meta_rep, rep_fa, meta_g, genome_fa ] }

    TBLASTN(
        pairs_ch.map { mr, rf, mg, gf -> [ mr, rf ] },
        pairs_ch.map { mr, rf, mg, gf -> [ mg, gf ] }
    )

    // Summarise TBLASTN hits — flag candidates with hits in outgroup genomes
    SUMMARIZE_TBLASTN(TBLASTN.out.collect { meta_pair, tsv -> tsv })

    emit:
    tblastn_hits    = TBLASTN.out
    summary         = SUMMARIZE_TBLASTN.out
}

process SPLIT_REPS {
    label 'low_cpu'

    input:
    path(reps_fa)

    output:
    path("rep_*.fa")

    script:
    """
    split_fasta.py --input ${reps_fa} --prefix rep_ --suffix .fa
    """
}

process SUMMARIZE_TBLASTN {
    label 'low_cpu'
    publishDir "${params.outdir}/${params.project}", mode: 'copy'

    input:
    path(tblastn_tsvs)

    output:
    path("tblastn_summary.tsv")

    script:
    """
    summarize_tblastn.py \
        --hits ${tblastn_tsvs} \
        --evalue ${params.evalue} \
        --output tblastn_summary.tsv
    """
}
