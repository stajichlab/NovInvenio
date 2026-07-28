nextflow.enable.dsl=2

// PROFILE_SEARCH — the family-profile alternative to the pairwise SEARCH workflow
// (ADR-0002). Clusters the SEED group into gene families, builds a family HMM per
// >=2-member family (famsa + hmmbuild), scans EVERY proteome (both groups) with the family
// HMM database, and turns the results into the same presence_matrix.tsv + candidates.txt
// contract SEARCH emits — so main.nf can branch on --cluster_tool and reuse the whole
// downstream chain unchanged.
//
// Direction-agnostic: the novelty direction seeds families from the ingroup (query-group
// IN); the loss mirror (issue #12) aliases this workflow as PROFILE_LOSS_SEARCH, seeding
// from the outgroup (query-group OUT) with distinct `loss_` output names so both directions
// coexist in the same publishDir.

include { MMSEQS_FAMILY_CLUSTER   } from '../modules/mmseqs_family_cluster'
include { BUILD_FAMILY_PROFILES   } from '../modules/build_family_profiles'
include { FAMILY_HMMSEARCH        } from '../modules/family_hmmsearch'
include { PROFILE_PRESENCE_MATRIX } from '../modules/profile_presence_matrix'

// Per-seed-proteome protein_id -> proteome Short map, merged into one TSV so
// profile_to_matrix.py can attribute each family member to its source proteome.
process SEED_PROTEIN_MAP {
    label 'low_cpu'
    tag "${meta.id}"

    input:
    tuple val(meta), path(proteome_fa)

    output:
    path("${meta.id}.protmap.tsv")

    script:
    """
    grep '^>' ${proteome_fa} \
        | sed 's/^>//; s/[[:space:]].*//' \
        | awk -v s='${meta.id}' '{print \$1"\\t"s}' > ${meta.id}.protmap.tsv
    """
}

workflow PROFILE_SEARCH {
    take:
    seed_ch         // [meta, protein_fa] — the group families are seeded from
                    //   (ingroup for novelty, outgroup for the loss mirror)
    other_ch        // [meta, protein_fa] — the other group (probed for presence too)
    config_csv      // path to analysis CSV
    query_group     // 'IN' (novelty) | 'OUT' (loss)
    query_min_frac  // presence threshold within the query group
    other_max_frac  // max fraction of the other group a candidate may survive in
    out_prefix      // '' (novelty) | 'loss_' — distinct output filenames/subdirs

    main:
    all_proteomes_ch = seed_ch.mix(other_ch)

    // Concatenate the seed group once (order irrelevant) and build the protein->proteome
    // map. .first() turns the single collected file into a value channel reusable by many
    // tasks. Names carry out_prefix so the loss alias's staged files don't collide.
    seed_concat = seed_ch.map { meta, fa -> fa }
                         .collectFile(name: "${out_prefix}seed_all.faa")
                         .first()

    protein_map = SEED_PROTEIN_MAP(seed_ch)
                      .collectFile(name: "${out_prefix}protein_to_proteome.tsv")
                      .first()

    // Cluster seed -> families; build a family HMM db (famsa + hmmbuild, >=N members).
    MMSEQS_FAMILY_CLUSTER(seed_concat, out_prefix)
    BUILD_FAMILY_PROFILES(MMSEQS_FAMILY_CLUSTER.out.cluster_tsv, seed_concat, out_prefix)

    // Scan every proteome (both groups) with the family HMM database.
    FAMILY_HMMSEARCH(all_proteomes_ch, BUILD_FAMILY_PROFILES.out.profiles.first(), out_prefix)

    // Presence matrix + candidates (same contract as BUILD_PRESENCE_MATRIX).
    PROFILE_PRESENCE_MATRIX(
        FAMILY_HMMSEARCH.out.domtblout.map { meta, dom -> dom }.collect(),
        MMSEQS_FAMILY_CLUSTER.out.cluster_tsv,
        BUILD_FAMILY_PROFILES.out.families.first(),
        protein_map,
        config_csv,
        query_group,
        query_min_frac,
        other_max_frac,
        out_prefix
    )

    emit:
    matrix          = PROFILE_PRESENCE_MATRIX.out.matrix
    candidates      = PROFILE_PRESENCE_MATRIX.out.candidates
    // Family clustering + reps + the concatenated seed group, for the family-as-cluster
    // path (ADR-0002 Q7): main.nf feeds these to PROFILE_CANDIDATE_CLUSTERS instead of
    // re-clustering the candidates.
    family_cluster_tsv = MMSEQS_FAMILY_CLUSTER.out.cluster_tsv
    family_reps        = MMSEQS_FAMILY_CLUSTER.out.representatives
    seed_concat        = seed_concat
}
