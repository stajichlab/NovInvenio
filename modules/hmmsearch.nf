nextflow.enable.dsl=2

process HMMSEARCH {
    label 'hmmsearch'
    tag "${cluster_id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/hmm_results" }, mode: 'copy'

    input:
    tuple val(cluster_id), path(hmm_file)
    path(target_db)   // SwissProt / UniProt FASTA

    output:
    tuple val(cluster_id), path("${cluster_id}.hmmsearch.tblout"), emit: tblout

    script:
    // MPI mode: mpirun -np N hmmsearch --mpi (N single-threaded MPI workers + 1 master)
    // Thread mode: hmmsearch --cpu N  (single process, N threads)
    def n_mpi    = params.hmm_mpi_tasks ?: params.max_cpus
    def mpi_cmd  = params.hmm_mpi ? "mpirun -np ${n_mpi}" : ""
    def cpu_flag = params.hmm_mpi ? "--mpi" : "--cpu ${task.cpus}"
    // SLURM_CPU_BIND=none: a multi-task MPI launch fails outright on this cluster without
    // it ("CPU binding outside of job step allocation") -- --cpu-bind itself is an srun-only
    // flag (sbatch rejects it), so it can't go in clusterOptions; exporting the env var
    // makes any srun mpirun invokes internally (PMI/PMIx under a SLURM allocation) pick it
    // up as its default. No-op in thread mode. See nextflow-hpcc skill.
    def bind_env = params.hmm_mpi ? "export SLURM_CPU_BIND=none" : ""
    """
    ${bind_env}
    ${mpi_cmd} hmmsearch ${cpu_flag} \
        --tblout ${cluster_id}.hmmsearch.tblout \
        --noali \
        -E ${params.evalue} \
        ${hmm_file} ${target_db} \
        > /dev/null
    """
}
