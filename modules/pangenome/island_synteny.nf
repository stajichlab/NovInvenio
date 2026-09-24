// ISLAND_SYNTENY -- per-island presence/absence grid (issue #116), rows =
// strains collapsed into distinct haplotypes, columns = member families in
// locus order, so a deletion breakpoint shows up as a vertical edge.
//
// Runs inside the --pangenome_island_pfam_hmm block because the islands
// table it consumes (REPORT_TABLES' islands_with_domains.tsv) only exists
// there. Degrades to a valid "no islands" page rather than failing when a
// study has none.
//
// --domtblout (Controller Ruling R8) is the same Pfam domtblout
// FAMILY_PFAM_SCAN already produces for DOMAIN_ENRICHMENT/REPORT_TABLES --
// this process just reuses it to colour the per-family glyph strip, no new
// scan.
//
// --config (issue #119) is the analysis samplesheet -- the same CSV
// PRESENCE_MATRIX already receives (modules/pangenome/presence_matrix.nf) --
// threaded through so pangenome_island_synteny.py can build a {Short:
// Species} map for the page's "by species" row sort. It carries a
// Species column keyed by Short; lib/config_parser.py::parse_config()
// already parses it, so this is a channel-wiring change, not a new
// derivation.
process ISLAND_SYNTENY {
    label 'low_cpu'
    tag "island_synteny"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(islands_with_domains)
    path(presence_matrix)
    path(family_positions)
    path(domtblout)
    path(samplesheet)
    path(diagnostics_banner_html)
    path(gene_positions)
    path(cluster_tsv)
    path(rescue_positions)

    output:
    path("island_synteny.html"), emit: page

    script:
    """
    pangenome_island_synteny.py \
        --islands_with_domains ${islands_with_domains} \
        --presence_matrix ${presence_matrix} \
        --family_positions ${family_positions} \
        --domtblout ${domtblout} \
        --domain_evalue ${params.pangenome_pfam_domain_evalue} \
        --project '${Helpers.projectName(params)}' \
        --min_strains ${params.pangenome_top_islands_min_strains} \
        --top_islands ${params.pangenome_viz_top_islands} \
        --config ${samplesheet} \
        --diagnostics_banner ${diagnostics_banner_html} \
        --gene_positions ${gene_positions} \
        --cluster_tsv ${cluster_tsv} \
        --rescue_positions ${rescue_positions} \
        --id_sep '${params.pangenome_id_sep}' \
        --output island_synteny.html
    """
}
