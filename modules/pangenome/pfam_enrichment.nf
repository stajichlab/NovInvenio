// SELECT_BACKGROUND_REPS / FAMILY_PFAM_SCAN / DOMAIN_ENRICHMENT -- Pfam
// functional-enrichment testing for accessory-island member families.
// ONE hmmscan over the full shell+cloud-eligible background (not a
// separate island-vs-background scan pair) -- island/background
// partitioning happens inside DOMAIN_ENRICHMENT itself, matching
// summarize_island_functions.py's actual logic (background must be a
// superset of island members for the Fisher test to be valid).
process SELECT_BACKGROUND_REPS {
    label 'low_cpu'
    tag "select_background_reps"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(rep_fasta)
    path(frequency_table)

    output:
    path("background_reps.fa"), emit: fasta

    script:
    """
    pangenome_select_background_reps.py \
        --rep_fasta ${rep_fasta} \
        --frequency_table ${frequency_table} \
        --output background_reps.fa
    """
}

process FAMILY_PFAM_SCAN {
    label 'med_cpu'
    tag "family_pfam_scan"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(background_reps_fasta)
    path(pfam_hmm)

    output:
    path("pfam.domtblout"), emit: domtblout

    script:
    """
    hmmscan --domtblout pfam.domtblout \
        -E ${params.pangenome_pfam_domain_evalue} \
        --cpu ${task.cpus} \
        ${pfam_hmm} ${background_reps_fasta} > /dev/null
    """
}

process DOMAIN_ENRICHMENT {
    label 'low_cpu'
    tag "domain_enrichment"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(significant_islands)
    path(domtblout)
    path(frequency_table)

    output:
    path("island_pfam_enrichment.tsv"), emit: enrichment

    script:
    """
    pangenome_domain_enrichment.py \
        --significant_islands ${significant_islands} \
        --domtblout ${domtblout} \
        --frequency_table ${frequency_table} \
        --output island_pfam_enrichment.tsv
    """
}
