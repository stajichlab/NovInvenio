nextflow.enable.dsl=2

// Match one proteome to UniProt records via the library index (bin/uniprot_link.py):
// UniProt accession, then RefSeq ID, then exact sequence (own species first, then any).
// Replaces UNIPROT_XREF (issue #92). Runs for every config proteome when --uniprot_index
// is set; see docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md.
process UNIPROT_LINK {
    label 'low_cpu'
    tag "${meta.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    tuple val(meta), path(protein_fa)
    val(index_dir)

    output:
    path("${meta.id}.uniprot_link.tsv"), emit: tsv

    script:
    def taxid    = meta.taxid ? "--taxid ${meta.taxid}" : ''
    def restrict = meta.uniprot_restrict ? "--restrict-proteome ${meta.uniprot_restrict}" : ''
    """
    uniprot_link.py --index ${index_dir} --protein-fasta ${protein_fa} \\
        --short ${meta.id} --species "${meta.species}" ${taxid} ${restrict} \\
        --output ${meta.id}.uniprot_link.tsv
    """
}
