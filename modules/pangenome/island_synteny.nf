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
        --output island_synteny.html
    """
}
