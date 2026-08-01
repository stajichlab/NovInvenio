// ANNOTATE_PFAM — hmmscan candidates.fa against the pressed Pfam-A DB, as a scatter-gather
// so it parallelises across the cluster instead of running one long single-task hmmscan.
// hmmscan's own --cpu threading was benchmarked against a real ~5.7k-protein candidate set
// and found to badly under-utilise cores here (~1.4 of 8 requested, GPFS-latency-bound
// re-scanning the full pressed Pfam DB per query) -- splitting the query FASTA into chunks
// and running each as its own low-cpu SLURM task sidesteps that bottleneck by getting
// parallelism from many concurrent tasks instead of from hmmscan's internal threading.
// Per-chunk cpus=2: a follow-up local-scratch benchmark (500-protein subset, DB copied to
// node-local disk first) still only reached ~1.65 real cores even with 8 requested, so more
// than 2 buys nothing for a chunk this size.
//
// scratch=true (see conf/ucr_hpcc_slurm.config) copies each task's inputs to node-local
// disk -- benchmarked at 2.2x wall-clock speedup for the hmmscan step itself (94% CPU/14m55s
// reading Pfam straight off GPFS vs 165% CPU/6m46s from node-local scratch, same 500-protein
// subset, same node). The DB copy itself is cheap in isolation (~3s for the ~2.5 GB of
// pressed index files this actually needs -- .h3f/.h3i/.h3m/.h3p; hmmscan does NOT read the
// 2.1 GB flat Pfam-A.hmm text file or the .dat file at all, confirmed by running hmmscan
// against a directory containing only the four pressed files), but that was measured in
// isolation. If N chunks all launch around the same time, N concurrent ~2.5 GB reads is a
// real burst against shared GPFS bandwidth (affects other cluster users too, not just this
// job) -- maxForks below bounds how many chunks copy the DB concurrently.
//
//   SPLIT_CANDIDATES (once)  → candidates.fa split into N-sequence chunks
//        │ combine each chunk with the Pfam-A pressed index files (.h3f/.h3i/.h3m/.h3p only)
//   PFAM_HMMSCAN_CHUNK (many)  → hmmscan one chunk vs Pfam-A → partial tblout
//        │ collect
//   MERGE_PFAM_TBLOUT (once)  → candidates.pfam.tblout (same shape the old single-task
//                                hmmscan produced, so annotate_presence_matrix.py is
//                                unchanged downstream)

process SPLIT_CANDIDATES {
    label 'low_cpu'
    tag "split_candidates"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(candidates_fa)
    val(chunk_size)

    output:
    path("pfam_chunk_*.fa"), emit: chunks

    shell:
    '''
    split_fasta_chunks.py \
        --input !{candidates_fa} \
        --chunk-size !{chunk_size} \
        --outdir . \
        --prefix pfam_chunk_
    '''
}

process PFAM_HMMSCAN_CHUNK {
    label 'hmmsearch'
    tag { chunk.baseName }
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    // Cap concurrent chunks -- each stages its own ~2.5 GB copy of the pressed Pfam DB
    // (see scratch=true), so uncapped concurrency means N simultaneous large GPFS reads.
    maxForks params.pfam_chunk_max_forks

    input:
    path(pfam_files)     // Pfam-A.hmm.h3f/.h3i/.h3m/.h3p ONLY -- not the flat .hmm or .dat,
                          // hmmscan against a pressed DB never reads either (confirmed).
    val(pfam_basename)   // 'Pfam-A.hmm' -- the shared prefix hmmscan resolves .h3* files
                          // from; the flat file need not actually exist on disk.
    path(chunk)

    output:
    path("${chunk.baseName}.pfam.tblout"), emit: tblout

    script:
    """
    if [ -s ${chunk} ]; then
        hmmscan --cpu ${task.cpus} \
            --tblout ${chunk.baseName}.pfam.tblout \
            --noali \
            ${pfam_basename} \
            ${chunk} > /dev/null
    else
        touch ${chunk.baseName}.pfam.tblout
    fi
    """
}

process MERGE_PFAM_TBLOUT {
    label 'low_cpu'
    tag "merge_pfam_tblout"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(chunk_tblouts)
    val(output_prefix)   // '' (novelty) | 'loss_'

    output:
    path("${output_prefix}candidates.pfam.tblout"), emit: tblout

    shell:
    '''
    cat !{chunk_tblouts} > !{output_prefix}candidates.pfam.tblout
    '''
}

// pfam_hmm was unset (params.pfam_hmm not provided) -- skip the scatter-gather entirely and
// hand ANNOTATE_MATRIX an empty tblout, same as the old single-task "touch" fallback.
process EMPTY_PFAM_TBLOUT {
    label 'low_cpu'
    tag "empty_pfam_tblout"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    val(output_prefix)

    output:
    path("${output_prefix}candidates.pfam.tblout"), emit: tblout

    shell:
    '''
    touch !{output_prefix}candidates.pfam.tblout
    '''
}

workflow ANNOTATE_PFAM {
    take:
    candidates_fa    // path: candidates.fa
    pfam_hmm         // val: absolute path to Pfam-A.hmm (pressed), or ''
    output_prefix    // val: '' (novelty) | 'loss_'

    main:
    def chunk_size = params.pfam_chunk_size as int
    // Only the pressed index files hmmscan actually reads (.h3f/.h3i/.h3m/.h3p) -- not the
    // 2.1 GB flat Pfam-A.hmm text file or the .dat file, which would otherwise nearly double
    // every chunk's DB copy for no benefit.
    def pfam_files    = files("${pfam_hmm}.h3?")
    def pfam_basename = file(pfam_hmm).name

    SPLIT_CANDIDATES(candidates_fa, chunk_size)
    PFAM_HMMSCAN_CHUNK(pfam_files, pfam_basename, SPLIT_CANDIDATES.out.chunks.flatten())
    MERGE_PFAM_TBLOUT(PFAM_HMMSCAN_CHUNK.out.tblout.collect(), output_prefix)

    emit:
    tblout = MERGE_PFAM_TBLOUT.out.tblout
}
