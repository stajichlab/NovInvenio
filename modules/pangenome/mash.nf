// MASH_SKETCH — whole-genome Mash sketch + all-vs-all distance matrix for a
// given strain set. A single reusable process (invoked twice by
// workflows/pangenome_profile.nf, aliased MASH_SKETCH_ALL / MASH_SKETCH_INGROUP)
// instead of each downstream script re-sketching independently in Python, per
// NEXTFLOW_MIGRATION_NOTES.md section D's "real duplicated work to fix" note.
//
// NOT run once for a single combined set covering both downstream consumers:
// DEREPLICATE needs the FULL strain set (dedup must catch a duplicate isolate
// regardless of ingroup/outgroup), but ASSIGN_CLADES specifically needs
// INGROUP-ONLY distances -- the originating study found that including even
// one outgroup species dominates the Mash distance structure (a between-
// species distance dwarfs any within-species distance) and silently collapses
// the within-species clade signal to a degenerate "outgroup vs everyone else"
// split. Sketching the full set once and handing that same matrix to both
// would reintroduce exactly that artifact, so this process is called twice
// (over two different strain-set channels) rather than truly once -- still
// removes the real duplication (two independent Python subprocess
// sketch/dist calls become one shared process implementation), just not a
// literal single invocation. See the module docstrings of
// bin/pangenome_dereplicate_strains.py / bin/pangenome_assign_clades.py.
process MASH_SKETCH {
    label 'med_cpu'
    tag "mash_sketch_${out_prefix}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(genome_fastas)
    val(out_prefix)

    output:
    path("${out_prefix}.mash_dist.tsv"), emit: dist_tsv

    script:
    """
    mash sketch -p ${task.cpus} -o ${out_prefix}_sketches ${genome_fastas}
    mash dist -p ${task.cpus} -t ${out_prefix}_sketches.msh ${out_prefix}_sketches.msh \
        > ${out_prefix}.mash_dist.tsv
    """
}

// Strain-inventory dedup (assembly-quality proxy + mash-threshold dedup
// groups) -- consumed by frequency binning / co-occurrence's optional
// dereplicated denominator.
process DEREPLICATE {
    label 'low_cpu'
    tag "dereplicate_strains"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(samplesheet)
    val(data_dir_abs)
    path(mash_dist_tsv)

    output:
    path("strain_inventory.tsv"), emit: inventory

    script:
    """
    pangenome_dereplicate_strains.py \
        --config ${samplesheet} --data_dir ${data_dir_abs} \
        --mash_threshold ${params.pangenome_mash_threshold} \
        --mash_dist_tsv ${mash_dist_tsv} \
        --output strain_inventory.tsv
    """
}

// Mash+PCoA+k-means clade assignment (ingroup-only) -- a working
// population-stratification label (samplesheet TaxonGroup) for
// co-occurrence's clade-permutation null when a study has no better
// DAPC/SNP-based clade scheme.
process ASSIGN_CLADES {
    label 'low_cpu'
    tag "assign_clades"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(samplesheet)
    val(data_dir_abs)
    path(mash_dist_tsv)

    output:
    path("clade_assignments.tsv"), emit: clades

    script:
    def k_fixed_arg = params.pangenome_clades_k_fixed ? "--k_fixed ${params.pangenome_clades_k_fixed}" : ''
    """
    pangenome_assign_clades.py \
        --config ${samplesheet} --data_dir ${data_dir_abs} \
        --groups '${params.pangenome_ingroup_label}' \
        --k_range '${params.pangenome_clades_k_range}' ${k_fixed_arg} \
        --n_pcoa_components ${params.pangenome_clades_n_pcoa} \
        --mash_dist_tsv ${mash_dist_tsv} \
        --output clade_assignments.tsv
    """
}

// Fill the samplesheet's empty TaxonGroup cells from ASSIGN_CLADES' Mash
// clade labels, WITHOUT overwriting any value a study already curated
// (e.g. a hand-entered DAPC clade). Closes a real gap the originating study
// left as a manual step (NEXTFLOW_MIGRATION_NOTES.md section B): without
// this, pangenome_cooccurrence.py's clade-permutation null would silently
// run unstratified for any fresh study (it reads TaxonGroup off the
// samplesheet directly, never off clade_assignments.tsv on its own). See
// bin/pangenome_fill_taxon_group.py's module docstring for the generic
// "fill only what's blank" policy this subworkflow ships with -- a study
// wanting a richer priority scheme (e.g. a published-table cross-reference)
// should fill TaxonGroup itself before this subworkflow runs.
process FILL_TAXON_GROUP {
    label 'low_cpu'
    tag "fill_taxon_group"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(samplesheet)
    path(clade_assignments)

    output:
    path("samplesheet.with_clades.csv"), emit: samplesheet

    script:
    """
    pangenome_fill_taxon_group.py \
        --config ${samplesheet} --clade_assignments ${clade_assignments} \
        --output samplesheet.with_clades.csv
    """
}
