// ISLAND_GBK_SLICE -- one GenBank file per (locus, strain) region in
// island_regions.tsv, for ISLAND_CLINKER (spec section 8). One task for all
// drawn loci. Reads each chosen strain's genome FASTA, GFF3 and protein
// FASTA from the study data_dir (passed as `val`, the GENE_POSITIONS
// convention -- needs `--bind /bigdata` under singularity).
// Measured 2026-09-26 on 3 loci x 12 strains of the 529-strain run: 35 s,
// 60 MB peak RSS.
process ISLAND_GBK_SLICE {
    label 'low_cpu'
    tag "island_gbk_slice"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy',
        pattern: 'gbk/island_slices.tsv'

    input:
    path(regions)
    path(samplesheet)
    val(data_dir_abs)
    val(gff3_dir_abs)
    path(gene_positions)
    path(rescue_positions)
    path(cluster_tsv)

    output:
    path("gbk/L*"), optional: true, emit: locus_dirs
    path("gbk/island_slices.tsv"), emit: slices

    script:
    """
    pangenome_island_gbk_slice.py \
        --regions ${regions} \
        --config ${samplesheet} \
        --data_dir ${data_dir_abs} \
        --gff3_dir ${gff3_dir_abs} \
        --gene_positions ${gene_positions} \
        --rescue_positions ${rescue_positions} \
        --cluster_tsv ${cluster_tsv} \
        --id_sep '${params.pangenome_id_sep}' \
        --out_dir gbk
    """
}
