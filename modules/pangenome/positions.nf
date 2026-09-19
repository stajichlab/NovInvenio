// GENE_POSITIONS — per-strain protein_id -> genomic position, parsed from
// each strain's own GFF3 CDS records. `gff3_dir_abs` is passed as a `val`
// (absolute path string), the same convention this repo already uses for
// pfam_hmm/swissprot_dmnd/data_dir (see main.nf/workflows/annotate.nf) --
// under -profile singularity this requires `--bind /bigdata` (already set in
// conf/ucr_hpcc_slurm.config) since Nextflow's autoMounts can't detect a
// path buried inside a `val` string.
process GENE_POSITIONS {
    label 'low_cpu'
    tag "gene_positions"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(samplesheet)
    val(gff3_dir_abs)
    val(protein_dir_abs)

    output:
    path("gene_positions.tsv.zst"), emit: positions

    script:
    """
    pangenome_build_gene_positions.py \
        --config ${samplesheet} --gff3_dir ${gff3_dir_abs} \
        --protein_dir ${protein_dir_abs} \
        --groups '${params.pangenome_ingroup_label},${params.pangenome_outgroup_label}' \
        --output gene_positions.tsv.zst
    """
}

// EXTRACT_RESCUE_POSITIONS — genomic positions for GENOME_ONLY (rescue-pass)
// presence calls, which have no GFF3-annotated protein_id to resolve via
// GENE_POSITIONS alone. Re-parses the same per-strain tblastn outputs
// RESCUE_PASS already consumed.
process EXTRACT_RESCUE_POSITIONS {
    label 'med_cpu'
    tag "extract_rescue_positions"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(rescued_matrix)
    path(tblastn_tsvs)

    output:
    path("rescue_positions.tsv"), emit: positions

    script:
    def tblastn_args = tblastn_tsvs.collect { "--tblastn_tsv ${it}" }.join(' ')
    """
    pangenome_extract_rescue_positions.py \
        --matrix ${rescued_matrix} \
        ${tblastn_args} \
        --min_pident ${params.pangenome_rescue_min_pident} \
        --min_qcov ${params.pangenome_rescue_min_qcov} \
        --processes ${task.cpus} \
        --output rescue_positions.tsv
    """
}

// FAMILY_POSITIONS — join gene_positions.tsv (+ optional rescue_positions.tsv)
// with the tier-1 cluster TSV into per-strain, per-family gene-order RANK
// positions -- the input pair-classification's physical-linkage test consumes.
process FAMILY_POSITIONS {
    label 'low_cpu'
    tag "family_positions"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(gene_positions)
    path(cluster_tsv)
    path(rescue_positions)   // may be a stub/empty file when rescue is disabled

    output:
    path("family_positions.tsv.zst"), emit: positions

    script:
    def rescue_arg = params.pangenome_rescue_enable ? "--rescue_positions ${rescue_positions}" : ''
    """
    pangenome_build_family_positions.py \
        --gene_positions ${gene_positions} --cluster_tsv ${cluster_tsv} \
        ${rescue_arg} \
        --id_sep '${params.pangenome_id_sep}' \
        --output family_positions.tsv.zst
    """
}
