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

// HMMPRESS_PFAM -- press the HMM database ONCE, up front, and pass the
// index set into every scan task.
//
// hmmscan requires the target HMM database to be hmmpress-indexed
// (.h3f/.h3i/.h3m/.h3p). Nextflow's `path()` input stages only the single
// named file into a task work dir, not any sibling index files that may
// exist alongside the source path -- so a real, already-pressed Pfam-A.hmm
// on disk still lands in a task without its indices. This used to be
// handled by pressing inside FAMILY_PFAM_SCAN itself, which was fine while
// that was one job; now that the scan is scattered across N chunk tasks
// (issue #112), pressing in-task would repeat ~30k profiles N times.
process HMMPRESS_PFAM {
    label 'low_cpu'
    tag "hmmpress_pfam"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(pfam_hmm)

    output:
    path("pressed/*"), emit: db

    script:
    // Copy rather than symlink: hmmpress writes its indices next to the
    // .hmm, and the staged input is a symlink into the caller's read-only
    // source directory.
    """
    mkdir -p pressed
    cp -L ${pfam_hmm} pressed/${pfam_hmm}
    hmmpress pressed/${pfam_hmm}
    """
}

// FAMILY_PFAM_SCAN -- one hmmscan per query chunk (issue #112). The caller
// scatters background_reps.fa with splitFasta(by: pangenome_pfam_chunk_size),
// so Nextflow submits one independent SLURM job per chunk; MERGE_PFAM_DOMTBLOUT
// reassembles the single pfam.domtblout every consumer still expects.
//
// Why chunk: as one monolithic job this step measured 1h50m-3h50m on the
// 529-strain genus_vs_ureesii study (22,061 query sequences vs 30,134 Pfam
// profiles), and twice FAILED at 1h59m against a 2-hour wall-clock cap --
// losing the whole scan each time. Chunked, a failure costs one chunk.
process FAMILY_PFAM_SCAN {
    label 'med_cpu'
    tag { "family_pfam_scan:${chunk_fasta.baseName}" }
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(chunk_fasta)
    path(pfam_db)

    output:
    path("${chunk_fasta.baseName}.domtblout"), emit: domtblout

    script:
    // The pressed set is staged as several files; the .hmm itself is the
    // one without an .h3* suffix.
    """
    HMM=\$(ls | grep -v '\\.h3[fimp]\$' | grep '\\.hmm\$' | head -1)
    hmmscan --domtblout ${chunk_fasta.baseName}.domtblout \
        --domE ${params.pangenome_pfam_domain_evalue} \
        --cpu ${task.cpus} \
        "\$HMM" ${chunk_fasta} > /dev/null
    """
}

// MERGE_PFAM_DOMTBLOUT -- reassemble the chunked scans into the single
// pfam.domtblout that DOMAIN_ENRICHMENT / MODULE_DOMAINS / REPORT_TABLES
// already consume. A plain `cat` is correct: hmmscan's domtblout carries
// its header and footer as `#`-prefixed lines, and lib-side
// parse_domtblout() skips every `#` line, so interleaved per-chunk
// headers are ignored rather than mis-parsed.
process MERGE_PFAM_DOMTBLOUT {
    label 'low_cpu'
    tag "merge_pfam_domtblout"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(chunk_domtblouts)

    output:
    path("pfam.domtblout"), emit: domtblout

    script:
    """
    cat ${chunk_domtblouts} > pfam.domtblout
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
