// EXTRACT_PROTEIN_DESCRIPTIONS — one-time pass over every proteome's peptide FASTA
// headers (UniProt-style: ">id description OS=... GN=... ...") producing a shared
// protein_id -> gene_name/description lookup TSV. Feeds two downstream consumers
// that otherwise only ever see a bare protein ID:
//   - BUILD_ALIGNMENT_SHARDS (workflows/validate.nf): embeds the *query* protein's
//     own name/description into its TBLASTN alignment shard entry.
//   - REPORT (workflows/report.nf): resolves a pairwise search hit's *target*
//     protein ID (bin/build_presence_matrix.py's --output-targets) into a display
//     name next to the existing "Hit e-values" evidence.
// Both are report-only, best-effort additions -- a protein ID with no matching
// FASTA header (or a non-UniProt-style header) just falls back to showing the
// bare ID, never an error.
process EXTRACT_PROTEIN_DESCRIPTIONS {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(pep_fastas)

    output:
    path('descriptions.tsv'), emit: tsv

    script:
    """
    extract_protein_descriptions.py \
        --pep ${pep_fastas} \
        --output descriptions.tsv
    """
}
