// FREQUENCY_BINS — core/soft-core/shell/cloud/singleton frequency binning
// over the ingroup (optionally dereplicated) strain set.
process FREQUENCY_BINS {
    label 'low_cpu'
    tag "frequency_bins"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(matrix)
    path(samplesheet)
    path(inventory)   // may be a stub/empty file when dereplication is disabled

    output:
    path("frequency_table.tsv"), emit: table

    script:
    def inventory_arg = Helpers.asBool(params.pangenome_dereplicate) ? "--inventory ${inventory}" : ''
    """
    pangenome_frequency_bins.py \
        --matrix ${matrix} --config ${samplesheet} \
        --ingroup_label '${params.pangenome_ingroup_label}' \
        ${inventory_arg} \
        --core_cutoff ${params.pangenome_core_cutoff} \
        --softcore_cutoff ${params.pangenome_softcore_cutoff} \
        --shell_cutoff ${params.pangenome_shell_cutoff} \
        --output frequency_table.tsv
    """
}

// COOCCURRENCE — Fisher/BH-FDR co-occurrence screen between shell+cloud
// families, with outgroup-based gain/loss polarization and an exact
// within-clade-stratified test on FDR survivors. Memory/time here should
// scale with (eligible-family-count)^2 (the candidate-pair count driver);
// see conf/ucr_hpcc_slurm.config's withName override for this process.
process COOCCURRENCE {
    label 'high_cpu'
    tag "cooccurrence"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(matrix)
    path(frequency_table)
    path(samplesheet)
    path(inventory)

    output:
    path("cooccurring_pairs.tsv.zst"), emit: pairs

    script:
    def inventory_arg = Helpers.asBool(params.pangenome_dereplicate) ? "--inventory ${inventory}" : ''
    """
    pangenome_cooccurrence.py \
        --matrix ${matrix} --frequency_table ${frequency_table} --config ${samplesheet} \
        --ingroup_label '${params.pangenome_ingroup_label}' \
        --outgroup_label '${params.pangenome_outgroup_label}' \
        ${inventory_arg} \
        --min_strain_count ${params.pangenome_min_strain_count} \
        --fdr_alpha ${params.pangenome_fdr_alpha} \
        --screen_alpha ${params.pangenome_screen_alpha} \
        --polarity_loss_min_frac ${params.pangenome_polarity_loss_min_frac} \
        --polarity_gain_max_frac ${params.pangenome_polarity_gain_max_frac} \
        --output cooccurring_pairs.tsv.zst
    """
}
