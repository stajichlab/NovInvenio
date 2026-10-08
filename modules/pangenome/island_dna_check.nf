// DNA presence check of the island locus view (spec
// docs/superpowers/specs/2026-09-24-island-locus-view-design.md, section 4b).
//
// ISLAND_DNA_TARGETS is pass 1 of bin/pangenome_island_loci.py: the same
// inputs and parameters as ISLAND_LOCI, plus --dna_targets_dir. It covers
// every candidate locus (--pangenome_candidates, default 200), not just the
// drawn --pangenome_top_loci: the drawn set is chosen AFTER the DNA check,
// from the DNA-informed ranking, so pass 1 and pass 2 must see identical
// candidate sets (changed 2026-09-27; previously only the pre-check
// --top_loci winners were checked, so a locus could never be promoted into
// view by a confirmed deletion). It writes one work list per
// --pangenome_locus_dna_batch loci: the exemplar's locus DNA and, per
// checked strain (flank intact, a locus column not in place), its DNA from
// the innermost left-flank gene's start to the innermost right-flank gene's
// end (flank genes included).
// An exemplar locus column that is a TBLASTN rescue hit takes its span from
// rescue_positions.tsv (start) and the exemplar's own per-strain tblastn
// output (end of the HSP at that start; plan Ruling R26). The tblastn files
// are staged under rescue_tblastn/ because, without the rescue pass, both
// rescue inputs are the same empty stub file name.
//
// ISLAND_DNA_CHECK runs blastn -task megablast (subject mode) per (locus,
// strain) on one work list. It reads each strain's genome FASTA from the
// study data_dir (passed as `val`, the GENE_POSITIONS convention -- needs
// `--bind /bigdata` under singularity) once per batch. Measured 2026-09-27
// (Task 19, SLURM 29118233, node r41, short partition) on the 529-strain
// Coccidioides run, one batch = all 50 drawn loci, 20159 (locus, strain)
// checks: 10 min 09 s wall at 16 cpus, 338 MB peak RSS -- well under the
// 1-1.5 h HPCC job-sizing target, so the default batch stays 50 (one task)
// and no explicit `time` is set. ISLAND_LOCI then applies the calls (pass 2).
// These numbers predate the 2026-09-27 fix above and cover only 50 loci;
// a real run now covers all --candidates (~4x more), so batch sizing is to
// be re-measured, not guessed.
process ISLAND_DNA_TARGETS {
    label 'low_cpu'
    tag "island_dna_targets"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

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
    path(rescue_positions)
    path(rescue_tblastn, stageAs: 'rescue_tblastn/*')

    output:
    path("dna_targets/batch_*.tsv"), optional: true, emit: batches

    script:
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
        --min_strains ${params.pangenome_top_islands_min_strains} \
        --dna_targets_dir dna_targets \
        --dna_batch ${params.pangenome_locus_dna_batch} \
        --rescue_positions ${rescue_positions} \
        --rescue_tblastn ${rescue_tblastn} \
        --project '${Helpers.projectName(params)}' \
        --output island_loci.pre_dna.json
    """
}

process ISLAND_DNA_CHECK {
    label 'med_cpu'
    tag "${targets.baseName}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome/island_dna_calls" }, mode: 'copy'

    input:
    path(targets)
    path(samplesheet)
    val(data_dir_abs)

    output:
    path("dna_calls_*.tsv"), emit: calls

    script:
    """
    pangenome_island_dna_check.py \
        --targets ${targets} \
        --config ${samplesheet} \
        --data_dir ${data_dir_abs} \
        --min_id ${params.pangenome_locus_dna_min_id} \
        --min_cov ${params.pangenome_locus_dna_min_cov} \
        --cpus ${task.cpus} \
        --output dna_calls_${targets.baseName}.tsv
    """
}
