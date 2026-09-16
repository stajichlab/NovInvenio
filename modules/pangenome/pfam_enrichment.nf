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
    // hmmscan requires the target HMM database to be hmmpress-indexed
    // (.h3f/.h3i/.h3m/.h3p). Nextflow's `path(pfam_hmm)` input only stages
    // the single named file into the task work dir, not any sibling index
    // files that may exist alongside the source path -- so a real,
    // already-pressed Pfam-A.hmm on disk still lands here without its
    // indices. Press it in the task's own work dir on demand rather than
    // requiring the caller to pass the index files through as extra
    // channel inputs (harder to wire, and not every HMM db a study might
    // point this at is guaranteed pre-pressed).
    """
    if [ ! -e ${pfam_hmm}.h3f ]; then
        hmmpress ${pfam_hmm}
    fi
    hmmscan --domtblout pfam.domtblout \
        --domE ${params.pangenome_pfam_domain_evalue} \
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
        --domain_evalue ${params.pangenome_pfam_domain_evalue} \
        --output island_pfam_enrichment.tsv
    """
}
