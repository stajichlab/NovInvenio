// ISLAND_CLINKER -- clinker figure per locus (spec section 8), batched
// --pangenome_clinker_batch loci per task (plan Ruling R10: one locus took
// 56 s wall on 4 cores for 12 strains x ~20 genes, 2026-09-26, so 50 loci
// fit one task of under an hour instead of 50 two-minute SLURM jobs).
//
// clinker is gamcil/clinker 0.0.32 from PyPI (pixi.toml
// [pypi-dependencies]). Bioconda's "clinker" 1.33 is an unrelated RNA-seq
// tool -- do not add it. -gf groups genes by tier-1 family, so clinker's
// groups are the grid's column IDs. bin/pangenome_island_clinker.py drops the
// embedded sequences (5.35 MB -> 1.01 MB on a real locus) unless
// --pangenome_clinker_slim false, and skips a locus whose clinker run fails
// (plan Ruling R15); the page then says "No synteny figure for this locus".
process ISLAND_CLINKER {
    label 'med_cpu'
    tag "island_clinker"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome/clinker" }, mode: 'copy'

    input:
    path(locus_dirs)

    output:
    path("L*.html"), optional: true, emit: html

    script:
    def keep = Helpers.asBool(params.pangenome_clinker_slim) ? '' : '--keep_sequences'
    """
    pangenome_island_clinker.py --locus_dirs ${locus_dirs} --cpus ${task.cpus} --out_dir . ${keep}
    """
}
