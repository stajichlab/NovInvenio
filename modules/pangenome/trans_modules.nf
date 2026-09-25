// LEIDEN_MODULES -- collapse the `trans` (physically unlinked, statistically
// significant) co-occurrence pairs from PAIR_CLASSIFICATION into interpretable
// gene-family modules via Leiden community detection. Runs unconditionally
// right after PAIR_CLASSIFICATION (cheap; degrades gracefully to empty output
// when a study has zero `trans` pairs -- see bin/pangenome_detect_trans_modules.py's
// module docstring for why that's a real, expected outcome and not an error).
process LEIDEN_MODULES {
    label 'low_cpu'
    tag "leiden_modules"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(pair_classification)

    output:
    path("family_modules.tsv"), emit: family_modules
    path("module_summary.tsv"), emit: module_summary

    script:
    """
    pangenome_detect_trans_modules.py \
        --pair_classification ${pair_classification} \
        --resolution ${params.pangenome_leiden_resolution} \
        --seed ${params.pangenome_leiden_seed} \
        --output_families family_modules.tsv \
        --output_modules module_summary.tsv
    """
}

// MODULE_DOMAINS -- per-module Pfam domain summary, cross-referencing
// LEIDEN_MODULES' family/module assignments against the same Pfam scan
// FAMILY_PFAM_SCAN already produced for island enrichment. Only runs inside
// the --pangenome_island_pfam_hmm block (needs that domtblout); a study with
// zero trans edges still produces a valid, empty summary (see
// bin/pangenome_module_domains.py).
process MODULE_DOMAINS {
    label 'low_cpu'
    tag "module_domains"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(family_modules)
    path(domtblout)

    output:
    path("module_domains.tsv"), emit: table

    script:
    """
    pangenome_module_domains.py \
        --family_modules ${family_modules} \
        --domtblout ${domtblout} \
        --min_module_size ${params.pangenome_module_min_size} \
        --output module_domains.tsv
    """
}

// MODULE_NEIGHBORHOOD -- View B1 (issue #182): for each Leiden module, how close
// its member genes sit inside each strain's own assembly, against a permutation
// null drawn from the same strain. Pairs on different contigs are excluded as
// uninformative, and the excluded fraction is reported. Positions are never
// compared across strains. Scale, null pool, n_perm and seed are written into
// every output row. See lib/pangenome_neighborhood.py.
process MODULE_NEIGHBORHOOD {
    label 'low_cpu'
    tag "module_neighborhood"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(family_modules)
    path(gene_positions)
    path(cluster_tsv)
    path(frequency_table)

    output:
    path("module_neighborhood.tsv"), emit: table

    script:
    """
    pangenome_neighborhood.py \
        --family_modules ${family_modules} \
        --gene_positions ${gene_positions} \
        --cluster_tsv ${cluster_tsv} \
        --frequency_table ${frequency_table} \
        --max_kb ${params.pangenome_neighborhood_max_kb} \
        --min_gene_gap ${params.pangenome_neighborhood_min_gene_gap} \
        --n_perm ${params.pangenome_neighborhood_n_perm} \
        --seed ${params.pangenome_neighborhood_seed} \
        --null_pool ${params.pangenome_neighborhood_null_pool} \
        --min_module_size ${params.pangenome_module_min_size} \
        --output module_neighborhood.tsv
    """
}
