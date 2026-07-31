nextflow.enable.dsl=2

include { MMSEQS_CLUSTER } from '../modules/mmseqs_cluster'
include { HMMBUILD       } from '../modules/hmmbuild'
include { HMMSEARCH      } from '../modules/hmmsearch'

workflow CLUSTER {
    take:
    candidates_file    // path: candidates.txt (one protein ID per line)
    prot_ch            // [meta, protein_fa] — to extract candidate sequences from
                        //   (ingroup for the novelty direction, outgroup for the loss
                        //   direction — whichever group candidates_file's source
                        //   proteomes were drawn from)
    config_csv         // path: analysis CSV (for short → FASTA filename mapping)
    candidates_fa_name // val: output FASTA filename, e.g. 'candidates.fa' or 'loss_candidates.fa'
    cluster_prefix      // val: mmseqs output prefix, e.g. 'clusters' or 'loss_clusters'

    main:
    EXTRACT_CANDIDATES(candidates_file, prot_ch.map { m, fa -> fa }.collect(), config_csv, candidates_fa_name)

    // Cluster with mmseqs2 (OrthoFinder support is a future addition)
    MMSEQS_CLUSTER(EXTRACT_CANDIDATES.out.candidates_fa, cluster_prefix)

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
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(candidates_txt)
    path(proteome_fastas)   // proteome FASTAs to search (see CLUSTER's prot_ch note)
    path(config_csv)
    val(output_name)        // 'candidates.fa' (novelty direction) or 'loss_candidates.fa'

    output:
    path("${output_name}"), emit: candidates_fa

    script:
    """
    extract_candidates.py \
        --candidates ${candidates_txt} \
        --fastas ${proteome_fastas} \
        --config ${config_csv} \
        --output ${output_name}
    """
}
