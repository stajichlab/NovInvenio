// PAIR_CLASSIFICATION — classify each FDR-significant co-occurring family
// pair as physically linked (captain-gene-explained or not) vs. trans
// (candidate non-physical interaction/co-evolution). See
// bin/pangenome_pair_classification.py's module docstring for the full
// label taxonomy and documented partial-implementation caveats.
process PAIR_CLASSIFICATION {
    label 'med_cpu'
    tag "pair_classification"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(cooccurring_pairs)
    path(family_positions)
    path(cluster_tsv)
    path(captain_tblout)   // may be an empty stub file when no captain-gene marker is used

    output:
    path("pair_classification.tsv"), emit: classification

    script:
    """
    pangenome_pair_classification.py \
        --cooccurring_pairs ${cooccurring_pairs} \
        --family_positions ${family_positions} \
        --cluster_tsv ${cluster_tsv} \
        --captain_tblout ${captain_tblout} \
        --k ${params.pangenome_pair_class_k} \
        --physical_threshold ${params.pangenome_pair_class_physical_threshold} \
        --trans_threshold ${params.pangenome_pair_class_trans_threshold} \
        --min_co_carrying ${params.pangenome_pair_class_min_co_carrying} \
        --perm_alpha ${params.pangenome_pair_class_perm_alpha} \
        --min_clades ${params.pangenome_pair_class_min_clades} \
        --output pair_classification.tsv
    """
}
