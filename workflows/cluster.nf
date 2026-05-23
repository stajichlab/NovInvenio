nextflow.enable.dsl=2

include { MMSEQS_CLUSTER } from '../modules/mmseqs_cluster'
include { HMMBUILD       } from '../modules/hmmbuild'
include { HMMSEARCH      } from '../modules/hmmsearch'

workflow CLUSTER {
    take:
    candidates_file   // path: candidates.txt (one protein ID per line)
    ingroup_prot_ch   // [meta, protein_fa] — to extract candidate sequences from
    config_csv        // path: analysis CSV (for short → FASTA filename mapping)

    main:
    // Extract candidate sequences from all ingroup proteomes
    EXTRACT_CANDIDATES(candidates_file, ingroup_prot_ch.map { m, fa -> fa }.collect(), config_csv)

    // Cluster with mmseqs2 (OrthoFinder support is a future addition)
    MMSEQS_CLUSTER(EXTRACT_CANDIDATES.out.candidates_fa)

    // TODO: per-cluster MSA → hmmbuild → hmmsearch against annotation db
    // Requires splitting the cluster output by cluster representative,
    // aligning each cluster with mafft or muscle, then HMMBUILD + HMMSEARCH.
    // Implement once mafft/muscle is added to pixi.toml.

    emit:
    candidates_fa   = EXTRACT_CANDIDATES.out.candidates_fa
    representatives = MMSEQS_CLUSTER.out.representatives
    cluster_tsv     = MMSEQS_CLUSTER.out.tsv
}

process EXTRACT_CANDIDATES {
    label 'low_cpu'
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(candidates_txt)
    path(proteome_fastas)   // all ingroup proteome FASTAs
    path(config_csv)

    output:
    path("candidates.fa"), emit: candidates_fa

    script:
    """
    extract_candidates.py \
        --candidates ${candidates_txt} \
        --fastas ${proteome_fastas} \
        --config ${config_csv} \
        --output candidates.fa
    """
}
