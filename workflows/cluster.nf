nextflow.enable.dsl=2

include { MMSEQS_CLUSTER } from '../modules/mmseqs_cluster'
include { HMMBUILD       } from '../modules/hmmbuild'
include { HMMSEARCH      } from '../modules/hmmsearch'

workflow CLUSTER {
    take:
    candidates_file   // path: candidates.txt (one protein ID per line)
    ingroup_prot_ch   // [meta, protein_fa] — to extract candidate sequences from

    main:
    // Extract candidate sequences from all ingroup proteomes
    EXTRACT_CANDIDATES(candidates_file, ingroup_prot_ch.map { m, fa -> fa }.collect())

    // Cluster with mmseqs2 (OrthoFinder support is a future addition)
    MMSEQS_CLUSTER(EXTRACT_CANDIDATES.out.candidates_fa)

    // TODO: per-cluster MSA → hmmbuild → hmmsearch against annotation db
    // Requires splitting the cluster output by cluster representative,
    // aligning each cluster with mafft or muscle, then HMMBUILD + HMMSEARCH.
    // Implement once mafft/muscle is added to pixi.toml.

    emit:
    representatives = MMSEQS_CLUSTER.out.representatives
    cluster_tsv     = MMSEQS_CLUSTER.out.tsv
}

process EXTRACT_CANDIDATES {
    label 'low_cpu'
    publishDir "${params.outdir}/${params.project}", mode: 'copy'

    input:
    path(candidates_txt)
    path(proteome_fastas)   // all ingroup proteome FASTAs

    output:
    path("candidates.fa"), emit: candidates_fa

    script:
    """
    extract_candidates.py \
        --candidates ${candidates_txt} \
        --fastas ${proteome_fastas} \
        --output candidates.fa
    """
}
