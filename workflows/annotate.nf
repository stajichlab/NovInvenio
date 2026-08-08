nextflow.enable.dsl=2

// Functional annotation of candidate proteins:
//   1. hmmsearch vs Pfam-A  (only when params.pfam_hmm is set)
//   2. diamond blastp vs SwissProt  (only when params.swissprot_dmnd is set)
//   3. annotate_presence_matrix.py merges all hits into the presence matrix,
//      using params.modelorgs_config YAML for gene name lookups
//
// hmmsearch, not hmmscan: benchmarked both directions on a real ~500-protein subset vs the
// full Pfam-A DB. hmmscan (candidates as query, Pfam pressed DB as target) badly under-used
// requested cores here (94% CPU / 14m55s reading Pfam off GPFS; 165% / 6m46s even after
// copying the DB to node-local scratch) -- apparently GPFS-latency-bound re-scanning the
// full pressed DB per query, and a chunked hmmscan scatter-gather was built and tested to
// work around that. Once SLURM_CPU_BIND=none was identified as the actual fix for MPI/multi-
// task srun launches on this cluster (see nextflow-hpcc skill -- --cpu-bind is an srun-only
// flag, not accepted by sbatch, so it has to be exported as an env var, not passed via
// clusterOptions), `mpirun -np N hmmsearch --mpi` (Pfam profiles dispatched across N-1
// workers, candidates as the small fixed target set) hit 860% CPU / 1m54s on the same
// subset -- ~4-8x faster than either hmmscan variant, and structurally simpler (one process,
// no chunking/staging machinery).

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

    // SLURM_CPU_BIND=none: MPI (-n > 1) launches fail outright on this cluster without it
    // ("CPU binding outside of job step allocation") -- --cpu-bind is an srun-only flag
    // (sbatch rejects it), so it's exported here rather than passed via clusterOptions.
    // No-op in thread mode. See nextflow-hpcc skill / modules/hmmsearch.nf.
    def bind_env = params.hmm_mpi ? "export SLURM_CPU_BIND=none" : ""

    // Guard on candidates_fa actually having content: an empty candidate set is a real,
    // legitimate outcome (e.g. a direction/lineage with zero clade-specific losses), not
    // an error -- but hmmsearch/hmmscan/diamond all hard-fail on an empty target/query
    // file rather than reporting zero hits, so touch the (empty) output instead of
    // invoking the search tool at all.
    def pfam_search = pfam_hmm ? """\
        if [ -s ${candidates_fa} ]; then
            ${bind_env}
            ${mpi_cmd} hmmsearch ${pfam_cpu} \\
                --tblout ${pfam_out} \\
                --noali \\
                ${pfam_hmm} \\
                ${candidates_fa} > /dev/null
        else
            touch ${pfam_out}
        fi
        """ : "touch ${pfam_out}"

    def sprot_search = swissprot_dmnd ? """\
        if [ -s ${candidates_fa} ]; then
            diamond blastp \\
                --query ${candidates_fa} \\
                --db ${swissprot_dmnd} \\
                --outfmt 6 qseqid sseqid stitle evalue bitscore \\
                --evalue 1e-5 \\
                --max-target-seqs 1 \\
                --threads ${diamond_threads} \\
                --out ${sprot_out}
        else
            touch ${sprot_out}
        fi
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
