// FAMILY_HMMSEARCH — scan the family HMM database against every proteome (ADR-0002 Q1/Q5),
// as a scatter-gather so the search parallelises across the cluster instead of running one
// unchunked job per proteome (issue #22 — the sibling fix to #14's BUILD_FAMILY_PROFILES
// chunking). At order scale (e.g. Chaetothyriales: 209 proteomes x 45k-62k families) a single
// proteome-vs-whole-db job can take hours, exceeding the `short` queue; splitting the family
// HMM db into chunks turns that into many short (proteome x chunk) jobs instead.
//
//   SPLIT_HMM_DB    (once)  → family_profiles.hmm split into N-HMM chunks
//        │ combine every proteome with every chunk
//   HMMSEARCH_CHUNK (many)  → hmmsearch one chunk against one proteome → partial domtblout
//        │ group by proteome, concatenate
//   MERGE_DOMTBLOUT (once per proteome) → <id>.family.domtblout (same as the old single-job output)
//
// `-Z` (params.hmm_Z) fixes the assumed TARGET (proteome) size for E-value comparability —
// it is per-target, not per-query-db, so splitting the HMM query db across chunks does not
// change its semantics. Concatenating a proteome's per-chunk domtblouts before parsing is
// equivalent to one unsplit run: each family (HMM) appears in exactly one chunk, so
// profile_to_matrix.py's per-(query,target) aggregation sees no duplicate or missing rows.

// Split the family HMM db into fixed-size chunks (by HMM record, i.e. `//`-terminated blocks).
process SPLIT_HMM_DB {
    label 'low_cpu'
    tag "split_hmm_db"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(profiles_hmm)
    val(chunk_size)

    output:
    path("hmm_chunk_*.hmm"), emit: chunks

    shell:
    '''
    if [ ! -s !{profiles_hmm} ]; then
        # Empty family HMM db (e.g. no profiled families) — one empty placeholder chunk so
        # every proteome still gets exactly one (empty) domtblout, matching the unchunked
        # behaviour for this edge case.
        touch hmm_chunk_000000.hmm
        exit 0
    fi
    awk -v chunk_size="!{chunk_size}" '
        BEGIN { idx = 0; n = 0; outfile = sprintf("hmm_chunk_%06d.hmm", idx) }
        {
            print > outfile
            if ($0 == "//") {
                n++
                if (n % chunk_size == 0) {
                    close(outfile)
                    idx++
                    outfile = sprintf("hmm_chunk_%06d.hmm", idx)
                }
            }
        }
        END { close(outfile) }
    ' !{profiles_hmm}
    '''
}

// hmmsearch one HMM-db chunk against one proteome. No storeDir: the output depends on the
// whole family HMM db (all seed-group families + params), which the filename cannot encode.
process HMMSEARCH_CHUNK {
    label 'hmmsearch'
    tag { "${meta.id}_${hmm_chunk.baseName}" }
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    tuple val(meta), path(proteome_fa), path(hmm_chunk)
    val(out_prefix)    // '' (novelty) | 'loss_' — only used to keep tags distinct in logs

    output:
    tuple val(meta), path("${meta.id}.${hmm_chunk.baseName}.domtblout"), emit: domtblout

    script:
    def zflag = params.hmm_Z ? "-Z ${params.hmm_Z}" : ""
    """
    if [ -s ${hmm_chunk} ]; then
        hmmsearch --cpu ${task.cpus} \
            --domtblout ${meta.id}.${hmm_chunk.baseName}.domtblout \
            --noali \
            -E ${params.hmm_report_evalue} ${zflag} \
            ${hmm_chunk} ${proteome_fa} > /dev/null
    else
        touch ${meta.id}.${hmm_chunk.baseName}.domtblout
    fi
    """
}

// Concatenate one proteome's per-chunk domtblouts into the final per-proteome file (same
// name/shape the old single-job FAMILY_HMMSEARCH produced).
process MERGE_DOMTBLOUT {
    label 'low_cpu'
    tag { meta.id }
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/${out_prefix}family_hmmsearch" }, mode: 'copy'

    input:
    tuple val(meta), path(chunk_domtblouts)
    val(out_prefix)    // '' (novelty) | 'loss_' — family_hmmsearch/ vs loss_family_hmmsearch/

    output:
    tuple val(meta), path("${meta.id}.family.domtblout"), emit: domtblout

    script:
    """
    cat ${chunk_domtblouts} > ${meta.id}.family.domtblout
    """
}

// Scatter-gather wrapper — same take/emit contract the single-process version had, so
// workflows/profile_search.nf's FAMILY_HMMSEARCH(...).out.domtblout callers are unchanged.
workflow FAMILY_HMMSEARCH {
    take:
    proteomes_ch    // [meta, protein_fa] — every proteome (both groups)
    profiles_hmm    // path: family_profiles.hmm (value channel)
    out_prefix      // '' (novelty) | 'loss_'

    main:
    def search_chunk_size = params.hmm_search_chunk_size as int
    SPLIT_HMM_DB(profiles_hmm, search_chunk_size)

    // Cross every proteome with every HMM chunk: [meta, proteome_fa, hmm_chunk].
    combo = proteomes_ch.combine(SPLIT_HMM_DB.out.chunks.flatten())

    HMMSEARCH_CHUNK(combo, out_prefix)

    // Group each proteome's chunk results back together and concatenate.
    grouped = HMMSEARCH_CHUNK.out.domtblout.groupTuple()
    MERGE_DOMTBLOUT(grouped, out_prefix)

    emit:
    domtblout = MERGE_DOMTBLOUT.out.domtblout
}
