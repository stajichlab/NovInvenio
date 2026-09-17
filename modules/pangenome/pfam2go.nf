nextflow.enable.dsl=2

// PFAM2GO -- optional GO-term annotation of island_pfam_enrichment.tsv via a
// standard pfam2go mapping file. Only invoked when params.pangenome_pfam2go
// is set (an operator-supplied local file -- see pangenome.nf's help text for
// where to obtain it; this pipeline does not auto-fetch it).
process PFAM2GO {
    label 'low_cpu'
    tag "pfam2go"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(island_pfam_enrichment)
    path(pfam2go)

    output:
    path("island_pfam_enrichment.go.tsv"), emit: annotated

    script:
    """
    pangenome_pfam2go.py \
        --island_pfam_enrichment ${island_pfam_enrichment} \
        --pfam2go ${pfam2go} \
        --output island_pfam_enrichment.go.tsv
    """
}
