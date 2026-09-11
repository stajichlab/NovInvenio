nextflow.enable.dsl=2

// Cross-walks an NCBI RefSeq-sourced proteome to its UniProt record via the DR RefSeq
// line (bin/build_uniprot_refseq_xref.py), so ANNOTATE can add uniprot_xrefs (and the
// other uniprot_*-prefixed columns) even when the pipeline's own protein_id is an NCBI
// accession, not a UniProt one. Only invoked for species with a configured
// UniProtDatGz column (main.nf's uniprot_xref_ch) -- a species with no UniProt
// proteome at all is simply absent from that channel, not an error.
//
// No storeDir: unlike the pairwise/self search modules this caches, a single .dat.gz
// parse is a cheap, fast, purely local operation (seconds, not the tool-runtime this
// repo's storeDir convention exists to avoid re-paying) -- ordinary -resume caching on
// the task's own inputs is enough.
process UNIPROT_XREF {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    tuple val(meta), path(protein_fa), path(dat_gz)

    output:
    path("${meta.id}.uniprot_xref.tsv"), emit: tsv

    script:
    """
    build_uniprot_refseq_xref.py \\
        --dat-gz ${dat_gz} \\
        --protein-fasta ${protein_fa} \\
        --short ${meta.id} \\
        --output ${meta.id}.uniprot_xref.tsv
    """
}
