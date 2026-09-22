// REPORT_TABLES / REPORT_RENDER -- tidy aggregation tables, then
// figures+Markdown, for the pangenome island+Pfam enrichment step. Split
// into two processes (not one) so table aggregation stays testable without
// a matplotlib dependency and independently reusable (e.g. by this repo's
// docs/ publishing pipeline).
process REPORT_TABLES {
    label 'low_cpu'
    tag "report_tables"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome/report_tables" }, mode: 'copy'

    input:
    path(significant_islands)
    path(island_pfam_enrichment)
    path(pair_classification)
    path(presence_matrix)
    path(frequency_table)
    path(domtblout)
    path(cluster_tsv)
    path(gene_positions)

    output:
    path("islands_with_domains.tsv"), emit: islands_with_domains
    path("island_size_distribution.tsv"), emit: size_distribution
    path("classification_counts.tsv"), emit: classification_counts
    path("marker_summary.tsv"), emit: marker_summary
    path("per_strain_summary.tsv"), emit: per_strain_summary

    script:
    """
    pangenome_report_tables.py \
        --significant_islands ${significant_islands} \
        --island_pfam_enrichment ${island_pfam_enrichment} \
        --pair_classification ${pair_classification} \
        --presence_matrix ${presence_matrix} \
        --frequency_table ${frequency_table} \
        --domtblout ${domtblout} \
        --domain_evalue ${params.pangenome_pfam_domain_evalue} \
        --cluster_tsv ${cluster_tsv} \
        --gene_positions ${gene_positions} \
        --id_sep '${params.pangenome_id_sep.replace("'", "'\\''")}' \
        --out_dir .
    """
}

process REPORT_RENDER {
    label 'low_cpu'
    tag "report_render"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(frequency_table)
    path(presence_matrix)
    path(islands_with_domains)
    path(island_size_distribution)
    path(classification_counts)
    path(island_pfam_enrichment)
    path(marker_summary)
    path(per_strain_summary)
    path(diagnostics_banner_md)

    output:
    path("report/report.md"), emit: report
    path("report/pangenome_openness.tsv"), emit: openness
    path("report/figures/*"), emit: figures_png
    path("report/figures_pdf/*"), emit: figures_pdf

    script:
    """
    pangenome_report_render.py \
        --frequency_table ${frequency_table} \
        --presence_matrix ${presence_matrix} \
        --islands_with_domains ${islands_with_domains} \
        --island_size_distribution ${island_size_distribution} \
        --classification_counts ${classification_counts} \
        --island_pfam_enrichment ${island_pfam_enrichment} \
        --marker_summary ${marker_summary} \
        --per_strain_summary ${per_strain_summary} \
        --top_islands_min_strains ${params.pangenome_top_islands_min_strains} \
        --n_permutations ${params.pangenome_accumulation_permutations} \
        --seed ${params.pangenome_accumulation_seed} \
        --diagnostics_banner ${diagnostics_banner_md} \
        --out_dir report
    """
}

// DIAGNOSTICS -- surfaces pipeline diagnostics (issue #134) into
// diagnostics.tsv (machine-readable) and Markdown/HTML banners consumed by
// REPORT_RENDER and ISLAND_SYNTENY, instead of leaving them only in
// RESCUE_PASS's stderr funnel counts. Advisory by default;
// params.pangenome_strict promotes any triggered diagnostic to a run
// failure (see bin/pangenome_diagnostics.py's docstring for the one
// exception: the zero-hit rescue guard, issue #126, which already always
// fails inside RESCUE_PASS itself and is not re-implemented here).
process DIAGNOSTICS {
    label 'low_cpu'
    tag "diagnostics"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(rescue_funnel)

    output:
    path("diagnostics/diagnostics.tsv"), emit: tsv
    path("diagnostics/diagnostics_banner.md"), emit: banner_md
    path("diagnostics/diagnostics_banner.html"), emit: banner_html

    script:
    // rescue_funnel is an empty placeholder file (EMPTY_EVALUES_STUB, same
    // convention as this pipeline's other optional-input stubs) when
    // params.pangenome_rescue_enable is false -- rescue_redundancy is then
    // reported not_computed rather than erroring.
    def strict_arg = params.pangenome_strict ? '--pangenome_strict' : ''
    def funnel_arg = (rescue_funnel.size() > 0) ? "--rescue_funnel ${rescue_funnel}" : ''
    """
    pangenome_diagnostics.py \
        ${funnel_arg} \
        --rescue_redundancy_threshold ${params.pangenome_rescue_redundancy_threshold} \
        ${strict_arg} \
        --out_dir diagnostics
    """
}
