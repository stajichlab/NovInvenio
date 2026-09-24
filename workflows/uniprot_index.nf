nextflow.enable.dsl=2

include { UNIPROT_PLAN_CHUNKS; UNIPROT_PARSE_CHUNK; UNIPROT_SEQ_INDEX } from '../modules/uniprot_index'

// One-time build of a UniProt library index (run via: nextflow run main.nf --build_uniprot_index; see main.nf).
workflow UNIPROT_INDEX_BUILD {
    take:
    library   // val: absolute path to the UniProt library directory

    main:
    UNIPROT_PLAN_CHUNKS(library)
    chunk_ids = UNIPROT_PLAN_CHUNKS.out.chunks
        .splitCsv(header: true, sep: '\t')
        .map { it.chunk_id }
        .unique()
    UNIPROT_PARSE_CHUNK(chunk_ids.combine(UNIPROT_PLAN_CHUNKS.out.chunks).combine(library))
    UNIPROT_SEQ_INDEX(UNIPROT_PARSE_CHUNK.out.records.flatten().collect(),
                      UNIPROT_PARSE_CHUNK.out.stats.flatten().collect(),
                      library)

    emit:
    manifest = UNIPROT_SEQ_INDEX.out.manifest
}
