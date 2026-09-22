// ASSEMBLY_QUALITY_QC — assembly-quality vs pangenome-content diagnostic
// (issue #130). A real 529-strain Coccidioides pangenome study found
// fragmented assemblies spuriously inflate accessory-family counts (a gene
// broken across a contig boundary yields two partial protein models that
// fail clustering and become fake strain-private families):
// rho(accessory ~ N50) = -0.53, rho(accessory ~ n_contigs) = +0.53. This
// process makes that confound visible on EVERY pangenome run, unconditionally
// (no gating param), rather than requiring a one-off manual analysis to
// rediscover it each time. Diagnostic only -- it never excludes a strain or
// alters a presence call (see bin/pangenome_assembly_quality_qc.py's module
// docstring for the full mechanism and the exact metrics reported).
//
// Computes its own assembly stats (n_contigs/N50/total_length) directly from
// each ingroup strain's DNA FASTA rather than depending on
// strain_inventory.tsv, so this runs the same whether or not
// `--pangenome_dereplicate` is enabled.
process ASSEMBLY_QUALITY_QC {
    label 'low_cpu'
    tag "assembly_quality_qc"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(samplesheet)
    val(data_dir_abs)
    path(rescued_matrix)
    path(frequency_table)
    path(gene_positions)
    path(cluster_tsv)

    output:
    path("assembly_quality_vs_content.tsv"), emit: table
    path("assembly_quality_correlations.tsv"), emit: correlations
    path("assembly_quality_report.md"), emit: report

    script:
    """
    pangenome_assembly_quality_qc.py \
        --config ${samplesheet} --data_dir ${data_dir_abs} \
        --matrix ${rescued_matrix} --frequency_table ${frequency_table} \
        --gene_positions ${gene_positions} --cluster_tsv ${cluster_tsv} \
        --ingroup_label '${params.pangenome_ingroup_label}' \
        --id_sep '${params.pangenome_id_sep}' \
        --terminus_window_bp ${params.pangenome_qc_terminus_window_bp} \
        --rho_warn_threshold ${params.pangenome_qc_rho_warn_threshold} \
        --out_dir .
    """
}
