// BUILD_FAMILY_PROFILES — turn seed-group gene families into a concatenated HMM database
// (ADR-0002 Q4/Q6), as a scatter-gather so the per-family famsa+hmmbuild work parallelises
// across the cluster instead of grinding serially in one task (issue #14). The seed group
// is the ingroup for novelty, the outgroup for the loss mirror (out_prefix disambiguates).
//
//   EXTRACT_FAMILY_SEQS  (once)  → per-family FASTAs + families.tsv
//        │ flatten + collate(params.family_chunk_size)
//   BUILD_CHUNK          (many)  → famsa + hmmbuild for a chunk of N families → partial HMM
//        │ collect
//   MERGE_PROFILES       (once)  → cat the partials → family_profiles.hmm
//
// The HMM name is the family representative id (hmmbuild -n <rep>), so a domtblout query
// name ties back to the cluster tsv. Singletons / sub-threshold families are dropped by
// extract_family_seqs.py. Each BUILD_CHUNK is short enough to fit the 2h `short` queue and
// carries the famsa AVX2 node constraint (see nextflow.config `.*BUILD_CHUNK`).

// Split the seed-group families into per-family FASTAs (one process, cheap).
process EXTRACT_FAMILY_SEQS {
    label 'low_cpu'
    tag "extract_families"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/${out_prefix}families" }, mode: 'copy',
        pattern: '{families.tsv,oversized_families.tsv}'

    input:
    path(cluster_tsv)
    path(seed_fa)      // concatenated seed-group proteomes (ingroup for novelty, outgroup for loss)
    val(out_prefix)    // '' (novelty) | 'loss_' — families/ vs loss_families/

    output:
    path("fam_*.faa"),              emit: family_fastas, optional: true
    path("families.tsv"),           emit: families
    path("oversized_families.tsv"), emit: oversized, optional: true

    shell:
    '''
    if [ ! -s !{cluster_tsv} ]; then
        printf 'family_index\\trepresentative_id\\tn_members\\n' > families.tsv
        exit 0
    fi
    extract_family_seqs.py \
        --cluster-tsv !{cluster_tsv} \
        --fasta !{seed_fa} \
        --min-members !{params.family_min_members} \
        --max-members !{params.family_max_members} \
        --oversized-report oversized_families.tsv \
        --outdir .
    '''
}

// Build profiles for one chunk of families (famsa + hmmbuild). Fans out; AVX2-constrained.
process BUILD_CHUNK {
    label 'high_cpu'
    tag { "chunk_${task.index}" }

    input:
    path(fam_files)       // a chunk: N per-family FASTAs (fam_*.faa)
    path(families_tsv)    // families.tsv, for the rep → HMM-name lookup

    output:
    path("chunk.*.hmm"), emit: partial_hmm

    shell:
    '''
    : > chunk.!{task.index}.hmm
    for fa in fam_*.faa; do
        [ -e "$fa" ] || continue
        base=$(basename "$fa" .faa)
        rep=$(awk -F'\\t' -v b="$base" '$1==b{print $2}' !{families_tsv})
        # Per-family safety net: --max-members should already exclude the pathological
        # multi-copy superfamilies that dominate famsa's runtime, but a rare outlier
        # (e.g. one unusually long/low-complexity sequence) could still stall. Bound each
        # family's alignment+build so one bad family can't cost the whole chunk a SLURM
        # walltime kill — skip it (logged) and keep going instead.
        if ! timeout !{params.family_align_timeout} famsa -t !{task.cpus} "$fa" "${fa}.aln" 2>/dev/null; then
            echo "WARN: famsa exceeded !{params.family_align_timeout}s or failed for $base ($rep) — skipping family" >&2
            continue
        fi
        if ! timeout !{params.family_align_timeout} hmmbuild --cpu !{task.cpus} -n "$rep" "${fa}.hmm" "${fa}.aln" > /dev/null; then
            echo "WARN: hmmbuild exceeded !{params.family_align_timeout}s or failed for $base ($rep) — skipping family" >&2
            continue
        fi
        cat "${fa}.hmm" >> chunk.!{task.index}.hmm
    done
    '''
}

// Concatenate the per-chunk partial HMMs into the family HMM database (one process, cheap).
process MERGE_PROFILES {
    label 'low_cpu'
    tag "merge_profiles"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/${out_prefix}families" }, mode: 'copy'

    input:
    path(partial_hmms)
    val(out_prefix)    // '' (novelty) | 'loss_' — families/ vs loss_families/

    output:
    path("family_profiles.hmm"), emit: profiles

    shell:
    '''
    : > family_profiles.hmm
    for f in chunk.*.hmm; do
        [ -e "$f" ] || continue
        cat "$f" >> family_profiles.hmm
    done
    '''
}

// Scatter-gather wrapper — same take/emit contract the single-process version had, so
// workflows/profile_search.nf and its BUILD_FAMILY_PROFILES.out.profiles/families callers
// are unchanged.
workflow BUILD_FAMILY_PROFILES {
    take:
    cluster_tsv
    seed_fa       // concatenated seed-group proteomes (ingroup for novelty, outgroup for loss)
    out_prefix    // '' (novelty) | 'loss_' — threaded to the publishing sub-processes

    main:
    EXTRACT_FAMILY_SEQS(cluster_tsv, seed_fa, out_prefix)

    // family_fastas emits every fam_*.faa (one list, or a lone file for a single family);
    // flatten to one item per family, then buffer into fixed-size chunks — each chunk is
    // one BUILD_CHUNK task. (buffer, not collate: collate is not a channel operator here.)
    // Coerce the size: a --family_chunk_size CLI override arrives as a String.
    def chunk_size = params.family_chunk_size as int
    chunks = EXTRACT_FAMILY_SEQS.out.family_fastas
                 .flatten()
                 .buffer( size: chunk_size, remainder: true )

    BUILD_CHUNK(chunks, EXTRACT_FAMILY_SEQS.out.families.first())

    // ifEmpty([]) so MERGE_PROFILES still emits an (empty) family_profiles.hmm when there
    // were no families to profile — FAMILY_HMMSEARCH depends on the file existing.
    MERGE_PROFILES(BUILD_CHUNK.out.partial_hmm.collect().ifEmpty([]), out_prefix)

    emit:
    profiles = MERGE_PROFILES.out.profiles
    families = EXTRACT_FAMILY_SEQS.out.families
}
