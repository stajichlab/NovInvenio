nextflow.enable.dsl=2

// Functional annotation of candidate proteins:
//   1. hmmsearch vs Pfam-A  (only when params.pfam_hmm is set)
//   2. diamond blastp vs SwissProt  (only when params.swissprot_dmnd is set)
//   3. annotate_presence_matrix.py merges all hits into the presence matrix,
//      using params.modelorgs_config YAML for gene name lookups

workflow ANNOTATE {
    take:
    candidates_fa    // path: candidates.fa (all ingroup candidate proteins)
    matrix           // path: presence_matrix.tsv
    pfam_hmm         // val: absolute path to Pfam-A.hmm, or ''
    swissprot_dmnd   // val: absolute path to SwissProt .dmnd, or ''
    morgs_config     // val: absolute path to modelorgs YAML, or ''

    main:
    ANNOTATE_MATRIX(candidates_fa, matrix, pfam_hmm, swissprot_dmnd, morgs_config)

    emit:
    annotated_matrix = ANNOTATE_MATRIX.out.matrix
    pfam_tblout      = ANNOTATE_MATRIX.out.pfam_tblout
    swissprot_tsv    = ANNOTATE_MATRIX.out.swissprot_tsv
}

// Single process handling Pfam search, SwissProt search, and matrix annotation.
// DB files and the modelorgs YAML are referenced by absolute (or CWD-relative)
// path from params; they are not staged into the work directory.
process ANNOTATE_MATRIX {
    label 'high_cpu'
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(candidates_fa)
    path(matrix)
    val(pfam_hmm)        // absolute path to Pfam-A.hmm, or ''
    val(swissprot_dmnd)  // absolute path to SwissProt .dmnd, or ''
    val(morgs_config)    // absolute path to modelorgs YAML, or ''

    output:
    path("presence_matrix.function.tsv"), emit: matrix
    path("candidates.pfam.tblout"),       emit: pfam_tblout
    path("candidates.swissprot.tsv"),     emit: swissprot_tsv

    script:
    def n_mpi           = params.hmm_mpi_tasks ?: params.max_cpus
    def mpi_cmd         = params.hmm_mpi ? "mpirun -np ${n_mpi}" : ""
    def pfam_cpu        = params.hmm_mpi ? "--mpi" : "--cpu ${task.cpus}"
    def diamond_threads = params.hmm_mpi ? n_mpi : task.cpus

    def pfam_search = pfam_hmm ? """\
        ${mpi_cmd} hmmsearch ${pfam_cpu} \\
            --tblout candidates.pfam.tblout \\
            --noali \\
            ${pfam_hmm} \\
            ${candidates_fa} > /dev/null
        """ : "touch candidates.pfam.tblout"

    def sprot_search = swissprot_dmnd ? """\
        diamond blastp \\
            --query ${candidates_fa} \\
            --db ${swissprot_dmnd} \\
            --outfmt 6 qseqid sseqid stitle evalue bitscore \\
            --evalue 1e-5 \\
            --max-target-seqs 1 \\
            --threads ${diamond_threads} \\
            --out candidates.swissprot.tsv
        """ : "touch candidates.swissprot.tsv"

    def pfam_flag  = pfam_hmm      ? '--pfam_hits candidates.pfam.tblout'           : ''
    def sprot_flag = swissprot_dmnd ? '--swissprot_hits candidates.swissprot.tsv'    : ''
    def morgs_flag = morgs_config  ? "--modelorgs_config ${morgs_config} --launch_dir ${projectDir}" : ''
    """
    ${pfam_search}
    ${sprot_search}
    annotate_presence_matrix.py \\
        --matrix ${matrix} \\
        --candidates_fa ${candidates_fa} \\
        ${morgs_flag} \\
        ${pfam_flag} \\
        ${sprot_flag} \\
        --output presence_matrix.function.tsv
    """
}
