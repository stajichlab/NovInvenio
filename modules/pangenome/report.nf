// REPORT_TABLES / REPORT_RENDER -- tidy aggregation tables, then
// figures+Markdown, for the pangenome report. Split into two processes (not
// one) so table aggregation stays testable without a matplotlib dependency
// and independently reusable (e.g. by this repo's docs/ publishing
// pipeline).
//
// Issue #135: both processes now run UNCONDITIONALLY, once per pipeline run
// (not just when `--pangenome_island_pfam_hmm` is set) -- classification
// counts and per-strain summaries only ever needed pair_classification/
// presence_matrix/frequency_table, all computed regardless of the
// islands+Pfam branch, so gating the whole report behind that branch meant
// most runs never got a report at all (see workflows/pangenome_profile.nf's
// header comment and issue #135). `significant_islands`/
// `island_pfam_enrichment`/`domtblout` are still genuinely islands+Pfam-only
// -- on a run where that branch didn't execute, the caller feeds these an
// empty (0-byte) stub file (workflows/pangenome_profile.nf's
// EMPTY_SIGNIFICANT_ISLANDS_STUB/EMPTY_ISLAND_ENRICHMENT_STUB/
// EMPTY_DOMTBLOUT_STUB, same convention as EMPTY_RESCUE_TSV_STUB etc.), and
// the `*_arg` conditionals below (`.size() > 0`, same pattern as
// DIAGNOSTICS's `funnel_arg`) omit the corresponding CLI flag entirely so
// pangenome_report_tables.py/pangenome_report_render.py treat it as "not
// computed" rather than trying to parse an empty file as real data.
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
    def islands_arg     = (significant_islands.size() > 0)     ? "--significant_islands ${significant_islands}"         : ''
    def enrichment_arg  = (island_pfam_enrichment.size() > 0)   ? "--island_pfam_enrichment ${island_pfam_enrichment}"   : ''
    def domtblout_arg   = (domtblout.size() > 0)                ? "--domtblout ${domtblout}"                            : ''
    """
    pangenome_report_tables.py \
        ${islands_arg} \
        ${enrichment_arg} \
        --pair_classification ${pair_classification} \
        --presence_matrix ${presence_matrix} \
        --frequency_table ${frequency_table} \
        ${domtblout_arg} \
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
    path(module_neighborhood)   // View B1 (issue #182)

    output:
    path("report/report.md"), emit: report
    path("report/pangenome_openness.tsv"), emit: openness
    path("report/figures/*"), emit: figures_png
    path("report/figures_pdf/*"), emit: figures_pdf

    script:
    def islands_arg    = (islands_with_domains.size() > 0)     ? "--islands_with_domains ${islands_with_domains}"           : ''
    def size_dist_arg  = (island_size_distribution.size() > 0) ? "--island_size_distribution ${island_size_distribution}"   : ''
    def enrichment_arg = (island_pfam_enrichment.size() > 0)   ? "--island_pfam_enrichment ${island_pfam_enrichment}"       : ''
    def marker_arg     = (marker_summary.size() > 0)           ? "--marker_summary ${marker_summary}"                       : ''
    def diagnostics_arg = (diagnostics_banner_md.size() > 0)   ? "--diagnostics_banner ${diagnostics_banner_md}"            : ''
    """
    pangenome_report_render.py \
        --frequency_table ${frequency_table} \
        --presence_matrix ${presence_matrix} \
        ${islands_arg} \
        ${size_dist_arg} \
        --classification_counts ${classification_counts} \
        ${enrichment_arg} \
        ${marker_arg} \
        --per_strain_summary ${per_strain_summary} \
        --top_islands_min_strains ${params.pangenome_top_islands_min_strains} \
        --n_permutations ${params.pangenome_accumulation_permutations} \
        --seed ${params.pangenome_accumulation_seed} \
        ${diagnostics_arg} \
        --module_neighborhood ${module_neighborhood} \
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
    path(assembly_correlations)

    output:
    path("diagnostics/diagnostics.tsv"), emit: tsv
    path("diagnostics/diagnostics_banner.md"), emit: banner_md
    path("diagnostics/diagnostics_banner.html"), emit: banner_html

    script:
    // rescue_funnel is an empty placeholder file (EMPTY_EVALUES_STUB, same
    // convention as this pipeline's other optional-input stubs) when
    // params.pangenome_rescue_enable is false -- rescue_redundancy is then
    // reported not_computed rather than erroring.
    def strict_arg = Helpers.asBool(params.pangenome_strict) ? '--pangenome_strict' : ''
    def funnel_arg = (rescue_funnel.size() > 0) ? "--rescue_funnel ${rescue_funnel}" : ''
    // assembly_correlations is ASSEMBLY_QUALITY_QC's correlations table
    // (issue #130). Its trip threshold reuses pangenome_qc_rho_warn_threshold
    // so this diagnostic and assembly_quality_report.md's WARNING agree.
    """
    pangenome_diagnostics.py \
        ${funnel_arg} \
        --assembly_correlations ${assembly_correlations} \
        --assembly_rho_threshold ${params.pangenome_qc_rho_warn_threshold} \
        --rescue_redundancy_threshold ${params.pangenome_rescue_redundancy_threshold} \
        ${strict_arg} \
        --out_dir diagnostics
    """
}
