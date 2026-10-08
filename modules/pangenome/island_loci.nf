// ISLAND_LOCI -- island locus view (spec
// docs/superpowers/specs/2026-09-24-island-locus-view-design.md): groups
// islands into loci, picks each locus's exemplar, computes every strain's
// cell states, row class and the breakpoint track, and writes
// island_loci.json for ISLAND_SYNTENY. With --pangenome_clinker it also
// writes island_regions.tsv (clinker strains and rank ranges, spec section 8).
//
// Streams family_positions three times, filtered each time (see
// bin/pangenome_island_loci.py). Measured on the 529-strain Coccidioides run
// (Task 13, 2026-09-26, 200 candidate loci, 50 drawn): 38 s wall, 1.36 GB
// peak RSS, 1.85 MB island_loci.json -- inside low_cpu's 4 GB.
//
// With --pangenome_locus_dna_check (default true) this is pass 2 of the
// DNA presence check (spec section 4b, modules/pangenome/island_dna_check.nf):
// dna_calls are ISLAND_DNA_CHECK's dna_calls_*.tsv plus an empty stub file,
// and dna_check is 'true'. With false, dna_calls is the stub alone.
//
// --assembly_quality is ASSEMBLY_QUALITY_QC's table; it covers the ingroup
// only, so outgroup strains fall back to locus-contig gene count for ties
// (lib/island_locus.py::quality_key, plan Ruling R4).
process ISLAND_LOCI {
    label 'low_cpu'
    tag "island_loci"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(islands_with_domains)
    path(presence_matrix)
    path(family_positions)
    path(frequency_table)
    path(assembly_quality)
    path(samplesheet)
    path(domtblout)
    path(gene_positions)
    path(cluster_tsv)
    path(rescue_positions, stageAs: 'rescue_positions.tsv')   // 0 bytes when rescue is off
    path(dna_calls)
    val(dna_check)

    output:
    path("island_loci.json"), emit: loci
    path("island_regions.tsv"), emit: regions

    script:
    def clinker_n = Helpers.asBool(params.pangenome_clinker) ? params.pangenome_clinker_max_strains : 0
    """
    pangenome_island_loci.py \
        --islands_with_domains ${islands_with_domains} \
        --presence_matrix ${presence_matrix} \
        --family_positions ${family_positions} \
        --frequency_table ${frequency_table} \
        --assembly_quality ${assembly_quality} \
        --config ${samplesheet} \
        --domtblout ${domtblout} \
        --domain_evalue ${params.pangenome_pfam_domain_evalue} \
        --gene_positions ${gene_positions} \
        --cluster_tsv ${cluster_tsv} \
        --id_sep '${params.pangenome_id_sep}' \
        --flank ${params.pangenome_locus_flank} \
        --flank_min ${params.pangenome_locus_flank_min} \
        --k ${params.pangenome_locus_k} \
        --empty_frac ${params.pangenome_locus_empty_frac} \
        --min_column_strains ${params.pangenome_locus_min_column_strains} \
        --island_reference_strains '${params.pangenome_island_reference_strains}' \
        --containment ${params.pangenome_locus_containment} \
        --rank_by ${params.pangenome_locus_rank} \
        --top_loci ${params.pangenome_top_loci} \
        --per_rank ${params.pangenome_locus_per_rank} \
        --poly_min_strains ${params.pangenome_locus_poly_min_strains} \
        --poly_min_frac ${params.pangenome_locus_poly_min_frac} \
        --poly_max_frac ${params.pangenome_locus_poly_max_frac} \
        --fixed_diff ${params.pangenome_locus_fixed_diff} \
        --candidates ${params.pangenome_locus_candidates} \
        --rescue_positions ${rescue_positions} \
        --dna_check ${dna_check} \
        --dna_calls ${dna_calls} \
        --dna_min_id ${params.pangenome_locus_dna_min_id} \
        --dna_min_cov ${params.pangenome_locus_dna_min_cov} \
        --min_strains ${params.pangenome_top_islands_min_strains} \
        --clinker_max_strains ${clinker_n} \
        --regions_out island_regions.tsv \
        --project '${Helpers.projectName(params)}' \
        --output island_loci.json
    """
}
