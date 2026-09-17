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
        --id_sep '${params.pangenome_id_sep}' \
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
        --n_permutations ${params.pangenome_accumulation_permutations} \
        --seed ${params.pangenome_accumulation_seed} \
        --out_dir report
    """
}
