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
