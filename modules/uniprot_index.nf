nextflow.enable.dsl=2

// UniProt library index build (docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md).
// Three processes in one file: they are only ever used together (workflows/uniprot_index.nf)
// and share one output contract, the index directory at params.uniprot_index.

process UNIPROT_PLAN_CHUNKS {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    val(library)

    output:
    path("chunks.tsv"), emit: chunks

    script:
    """
    uniprot_plan_chunks.py --library ${library} --library-csv ${params.uniprot_library_csv} \\
        --chunk-gb ${params.uniprot_index_chunk_gb} --output chunks.tsv
    """
}

process UNIPROT_PARSE_CHUNK {
    label 'low_cpu'
    tag "${chunk_id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    tuple val(chunk_id), path(chunks_tsv), val(library)

    output:
    path("records/*.records.tsv.zst"), emit: records
    path("stats/*.json"),              emit: stats

    script:
    """
    uniprot_parse_dat.py --library ${library} --chunks ${chunks_tsv} \\
        --chunk-id ${chunk_id} --outdir .
    """
}

process UNIPROT_SEQ_INDEX {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    storeDir { params.uniprot_index }

    input:
    path(records, stageAs: 'in_records/*')
    path(stats,   stageAs: 'in_stats/*')
    val(library)

    output:
    path("manifest.json"),    emit: manifest
    path("seq_index.sqlite"), emit: sqlite
    path("records"),          emit: records

    script:
    """
    mkdir -p records
    cp -L in_records/* records/
    uniprot_build_index.py --records-dir records --stats-dir in_stats \\
        --library ${library} --library-csv ${params.uniprot_library_csv} --output-dir .
    """
}
