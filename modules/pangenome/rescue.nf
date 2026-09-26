// Genome-level tblastn rescue pass: catches coverage-failed protein-model
// absences (fragmented/split gene models, draft-assembly annotation gaps)
// that protein-level clustering alone misses.
//
// Per-strain scatter (one BLAST database per strain, restricted to that
// strain's own currently-ABSENT family queries) -- NOT a shared combined
// genome database. The originating study's first version used one combined
// database across all strains with a shared `-max_target_seqs`, which was
// found to silently truncate hits ACROSS strains per query (89.4% of hit
// queries landed at exactly the cap in a real 295-strain run) -- a real
// correctness bug, not a style concern. The per-strain design removes the
// cross-strain competition for hits entirely: each database is tiny (one
// strain's genome), so a generous `--rescue_max_target_seqs` is cheap and
// only needs to be large enough to catch real intra-genome duplication.
// See NEXTFLOW_MIGRATION_NOTES.md section C for the full history.

// Extract, once, the per-strain FASTA of family representatives currently
// ABSENT in that strain -- the per-strain rescue query set.
process EXTRACT_ABSENT_QUERIES {
    label 'low_cpu'
    tag "extract_absent_queries"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(matrix)
    path(rep_fasta)

    output:
    path("per_strain_queries/*.absent.fa"), optional: true, emit: query_fastas
    path("per_strain_queries/manifest.tsv"), emit: manifest

    script:
    """
    pangenome_extract_absent_family_queries.py \
        --matrix ${matrix} --rep_fasta ${rep_fasta} \
        --out_dir per_strain_queries
    """
}

// Build one small nucleotide BLAST database per strain from its own
// Short-prefixed genome FASTA (see modules/pangenome/prefix_and_cluster.nf's
// PREFIX_GENOME) -- storeDir-cached like the novelty/loss pathway's
// TBLASTN_MAKEDB (modules/tblastn.nf), keyed by meta.id.
process MAKE_STRAIN_GENOME_DB {
    label 'low_cpu'
    tag "${meta.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    storeDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome/rescue/db" }

    input:
    tuple val(meta), path(genome_fa)

    output:
    tuple val(meta), path("${meta.id}.genome_db.*")

    script:
    """
    makeblastdb -in ${genome_fa} -dbtype nucl -out ${meta.id}.genome_db 2>/dev/null
    """
}

// tblastn: this strain's absent-family queries vs. this strain's own genome
// database only. Output is zstd-compressed as produced (never touching disk
// uncompressed), matching this repo's large-text-output compression default.
process TBLASTN_PER_STRAIN {
    label 'med_cpu'
    tag "${meta.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome/rescue/per_strain_chunks" }, mode: 'copy'

    input:
    tuple val(meta), path(genome_db), path(query_fa)

    output:
    tuple val(meta), path("${meta.id}.tblastn.tsv.zst"), emit: tsv

    script:
    """
    tblastn \
        -query ${query_fa} \
        -db ${meta.id}.genome_db \
        -outfmt "6 std qcovs" \
        -evalue ${params.pangenome_rescue_evalue} \
        -max_target_seqs ${params.pangenome_rescue_max_target_seqs} \
        -num_threads ${task.cpus} \
        | zstd -T${task.cpus} -o ${meta.id}.tblastn.tsv.zst
    """
}

// Fold every strain's tblastn hits back into the presence matrix, upgrading
// ABSENT calls to GENOME_ONLY wherever a qualifying hit exists.
//
// gene_positions/cluster_tsv/rep_fasta wire in the issue #133 structural
// rescue filter: reject a hit whose span overlaps a predicted gene already
// assigned to a DIFFERENT family in that strain (77.6% of all rescuable
// cells, exhaustively measured -- the same locus double-counted, not a real
// annotation dropout), whose query family rep is implausibly short, or that
// falls in a repeat-hotspot window. `gene_positions` therefore now runs
// BEFORE this process in workflows/pangenome_profile.nf (it always could --
// it only depends on the samplesheet/GFF3s, never on clustering or the
// rescue pass itself).
process RESCUE_PASS {
    label 'low_cpu'
    tag "rescue_pass"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(matrix)
    path(tblastn_tsvs)
    path(gene_positions)
    path(cluster_tsv)
    path(rep_fasta)

    output:
    path("presence_matrix.rescued.tsv"), emit: matrix
    path("presence_matrix.rescued.copy_number.tsv"), optional: true, emit: copy_number
    path("rescue_funnel.tsv"), emit: funnel

    script:
    def tblastn_args = tblastn_tsvs.collect { "--tblastn_tsv ${it}" }.join(' ')
    def structural_args = Helpers.asBool(params.pangenome_rescue_structural_filter)
        ? "--gene_positions ${gene_positions} --cluster_tsv ${cluster_tsv} --rep_fasta ${rep_fasta} " +
          "--rescue_min_rep_length ${params.pangenome_rescue_min_rep_length} " +
          "--rescue_hotspot_window ${params.pangenome_rescue_hotspot_window} " +
          "--rescue_hotspot_min_families ${params.pangenome_rescue_hotspot_min_families}"
        : ''
    """
    pangenome_rescue_pass.py \
        --matrix ${matrix} \
        ${tblastn_args} \
        --min_pident ${params.pangenome_rescue_min_pident} \
        --min_qcov ${params.pangenome_rescue_min_qcov} \
        --id_sep '${params.pangenome_id_sep}' \
        ${structural_args} \
        --funnel_tsv rescue_funnel.tsv \
        --output presence_matrix.rescued.tsv
    """
}
