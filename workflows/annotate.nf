nextflow.enable.dsl=2

// Functional annotation of candidate proteins:
//   1. hmmscan vs Pfam-A  (only when params.pfam_hmm is set)
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
    ANNOTATE_MATRIX(candidates_fa, matrix, pfam_hmm, swissprot_dmnd, morgs_config, output_prefix)

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
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(candidates_fa)
    path(matrix)
    val(pfam_hmm)        // absolute path to Pfam-A.hmm, or ''
    val(swissprot_dmnd)  // absolute path to SwissProt .dmnd, or ''
    val(morgs_config)    // absolute path to modelorgs YAML, or ''
    val(output_prefix)   // '' or 'loss_'

    output:
    path("${output_prefix}presence_matrix.function.tsv"), emit: matrix
    path("${output_prefix}candidates.pfam.tblout"),       emit: pfam_tblout
    path("${output_prefix}candidates.swissprot.tsv"),     emit: swissprot_tsv

    script:
    def n_mpi           = params.hmm_mpi_tasks ?: params.max_cpus
    def mpi_cmd         = params.hmm_mpi ? "mpirun -np ${n_mpi}" : ""
    def pfam_cpu        = params.hmm_mpi ? "--mpi" : "--cpu ${task.cpus}"
    def diamond_threads = params.hmm_mpi ? n_mpi : task.cpus
    def pfam_out        = "${output_prefix}candidates.pfam.tblout"
    def sprot_out        = "${output_prefix}candidates.swissprot.tsv"

    // hmmscan, not hmmsearch: candidates.fa (thousands of query proteins) is scanned
    // against the pressed Pfam-A DB (target). hmmsearch would instead treat each of
    // Pfam's ~20,800 profiles as a query against a small candidate set — the wrong
    // direction for multithreading (which parallelizes over the target), which is why
    // this step was observed at ~2.8 real cores no matter how many were requested
    // (issue: ANNOTATE_MATRIX runtime profiling). Pfam-A.hmm must be hmmpress'd (it is,
    // see db/pfam/).
    def pfam_search = pfam_hmm ? """\
        ${mpi_cmd} hmmscan ${pfam_cpu} \\
            --tblout ${pfam_out} \\
            --noali \\
            ${pfam_hmm} \\
            ${candidates_fa} > /dev/null
        """ : "touch ${pfam_out}"

    def sprot_search = swissprot_dmnd ? """\
        diamond blastp \\
            --query ${candidates_fa} \\
            --db ${swissprot_dmnd} \\
            --outfmt 6 qseqid sseqid stitle evalue bitscore \\
            --evalue 1e-5 \\
            --max-target-seqs 1 \\
            --threads ${diamond_threads} \\
            --out ${sprot_out}
        """ : "touch ${sprot_out}"

    def pfam_flag  = pfam_hmm      ? "--pfam_hits ${pfam_out}"   : ''
    def sprot_flag = swissprot_dmnd ? "--swissprot_hits ${sprot_out}" : ''
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
        --output ${output_prefix}presence_matrix.function.tsv
    """
}
