// BUILD_ISLANDS -- accessory-island construction + statistical significance
// gate (pangenome_build_islands.py). Named marker hit tables (0+, e.g.
// captain/sm_backbone) are passed as two PARALLEL lists -- marker_names
// (val) and marker_tblout_files (path, so Nextflow actually stages them
// into this task's work directory) -- rather than one combined list of
// tuples, to sidestep Nextflow's default flattening behavior on
// channel.collect() over tuples (a real bug an earlier draft of this task
// had: [[n1,p1],[n2,p2]].collect() flattens to [n1,p1,n2,p2] unless you
// explicitly ask for a nested list). See
// notes/superpowers/specs/2026-09-16-pangenome-island-pfam-enrichment-design.md.
process BUILD_ISLANDS {
    label 'low_cpu'
    tag "build_islands"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(family_positions)
    path(frequency_table)
    path(pair_classification)
    path(cluster_tsv)
    val(marker_names)             // parallel list of marker names, may be empty
    path(marker_tblout_files)     // parallel list of tblout files, same length, may be empty

    output:
    path("significant_islands.tsv"), emit: islands

    script:
    def marker_pairs = [marker_names, marker_tblout_files].transpose()
    def marker_args = marker_pairs.collect { name, path -> "--marker_tblout ${name}=${path}" }.join(' ')
    """
    pangenome_build_islands.py \
        --family_positions ${family_positions} \
        --frequency_table ${frequency_table} \
        --pair_classification ${pair_classification} \
        --cluster_tsv ${cluster_tsv} \
        --min_island_size ${params.pangenome_island_min_size} \
        ${marker_args} \
        --output significant_islands.tsv
    """
}

// MARKER_HMMSEARCH -- new, small, dedicated process for named marker
// searches (e.g. captain/DUF3435, sm_backbone/PKS-NRPS, or any future
// named marker per notes/superpowers/specs/
// 2026-09-16-pangenome-island-pfam-enrichment-design.md). Deliberately
// NOT an alias of modules/pangenome/captain.nf's CAPTAIN_HMMSEARCH (that
// process's output filename is hardcoded and it can't be invoked more
// than once per workflow) -- this is genuinely new, tiny code: same
// hmmsearch command, `tuple` input/output so Nextflow's channel-based
// multiplicity runs it once per marker via ONE process invocation (Task 8
// wires this with Channel.fromList, not a loop), and a name-templated
// output filename.
process MARKER_HMMSEARCH {
    label 'med_cpu'
    tag { marker_name }
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome/markers" }, mode: 'copy'

    input:
    tuple val(marker_name), path(marker_hmm)
    path(all_strains_fa)

    output:
    tuple val(marker_name), path("${marker_name}_vs_study.tblout"), emit: result

    script:
    """
    hmmsearch --tblout ${marker_name}_vs_study.tblout \
        -E ${params.pangenome_marker_evalue} \
        --cpu ${task.cpus} ${marker_hmm} ${all_strains_fa} > /dev/null
    """
}
