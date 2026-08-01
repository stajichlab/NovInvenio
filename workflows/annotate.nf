nextflow.enable.dsl=2

include { ANNOTATE_PFAM; EMPTY_PFAM_TBLOUT } from '../modules/annotate_pfam'

// Functional annotation of candidate proteins:
//   1. ANNOTATE_PFAM: hmmscan vs Pfam-A, chunked scatter-gather (only when params.pfam_hmm
//      is set) — see modules/annotate_pfam.nf for why this is chunked rather than one
//      single-task hmmscan.
//   2. diamond blastp vs SwissProt  (only when params.swissprot_dmnd is set)
//   3. annotate_presence_matrix.py merges all hits into the presence matrix,
//      using params.modelorgs_config YAML for gene name lookups

workflow ANNOTATE {
    take:
    candidates_fa    // path: candidates.fa (all candidate proteins)
    matrix           // path: presence_matrix.tsv
    pfam_hmm         // val: absolute path to Pfam-A.hmm, or ''
    swissprot_dmnd   // val: absolute path to SwissProt .dmnd, or ''
    morgs_config     // val: absolute path to modelorgs YAML, or ''
    output_prefix    // val: '' (novelty direction) or 'loss_' (loss direction) —
                      //   prefixes all three output filenames so a second ANNOTATE
                      //   call does not overwrite the first

    main:
    if (pfam_hmm) {
        ANNOTATE_PFAM(candidates_fa, pfam_hmm, output_prefix)
        pfam_tblout_ch = ANNOTATE_PFAM.out.tblout
    } else {
        EMPTY_PFAM_TBLOUT(output_prefix)
        pfam_tblout_ch = EMPTY_PFAM_TBLOUT.out.tblout
    }

    ANNOTATE_MATRIX(candidates_fa, matrix, pfam_tblout_ch, pfam_hmm, swissprot_dmnd,
                     morgs_config, output_prefix)

    emit:
    annotated_matrix = ANNOTATE_MATRIX.out.matrix
    pfam_tblout      = pfam_tblout_ch
    swissprot_tsv    = ANNOTATE_MATRIX.out.swissprot_tsv
}

// SwissProt search + final matrix annotation. The Pfam scan happens upstream in
// ANNOTATE_PFAM (chunked scatter-gather); this process just consumes its merged tblout.
// DB files and the modelorgs YAML are referenced by absolute (or CWD-relative)
// path from params; they are not staged into the work directory.
process ANNOTATE_MATRIX {
    label 'med_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(candidates_fa)
    path(matrix)
    path(pfam_tblout)
    val(pfam_hmm)        // absolute path to Pfam-A.hmm, or '' — only used for the --pfam_hits flag
    val(swissprot_dmnd)  // absolute path to SwissProt .dmnd, or ''
    val(morgs_config)    // absolute path to modelorgs YAML, or ''
    val(output_prefix)   // '' or 'loss_'

    output:
    path("${output_prefix}presence_matrix.function.tsv"), emit: matrix
    path("${output_prefix}candidates.swissprot.tsv"),     emit: swissprot_tsv

    script:
    def sprot_out = "${output_prefix}candidates.swissprot.tsv"

    def sprot_search = swissprot_dmnd ? """\
        diamond blastp \\
            --query ${candidates_fa} \\
            --db ${swissprot_dmnd} \\
            --outfmt 6 qseqid sseqid stitle evalue bitscore \\
            --evalue 1e-5 \\
            --max-target-seqs 1 \\
            --threads ${task.cpus} \\
            --out ${sprot_out}
        """ : "touch ${sprot_out}"

    def pfam_flag  = pfam_hmm      ? "--pfam_hits ${pfam_tblout}"   : ''
    def sprot_flag = swissprot_dmnd ? "--swissprot_hits ${sprot_out}" : ''
    def morgs_flag = morgs_config  ? "--modelorgs_config ${morgs_config} --launch_dir ${projectDir}" : ''
    """
    ${sprot_search}
    annotate_presence_matrix.py \\
        --matrix ${matrix} \\
        --candidates_fa ${candidates_fa} \\
        ${morgs_flag} \\
        ${pfam_flag} \\
        ${sprot_flag} \\
        --output ${output_prefix}presence_matrix.function.tsv
    """
}
