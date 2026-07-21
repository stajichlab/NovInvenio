nextflow.enable.dsl=2

// PROFILE_SEARCH — the family-profile alternative to the pairwise SEARCH workflow
// (ADR-0002). Clusters the ingroup into gene families, builds a family HMM per >=2-member
// family (famsa + hmmbuild), scans EVERY proteome (ingroup + outgroup) with the family HMM
// database, and turns the results into the same presence_matrix.tsv + candidates.txt
// contract SEARCH emits — so main.nf can branch on --cluster_tool and reuse the whole
// downstream chain unchanged.
//
// Novelty direction only in this first cut; the outgroup-seeded loss mirror is a follow-up
// (the ingroup-direction pairwise LOSS_SEARCH still runs meanwhile).

include { MMSEQS_FAMILY_CLUSTER   } from '../modules/mmseqs_family_cluster'
include { BUILD_FAMILY_PROFILES   } from '../modules/build_family_profiles'
include { FAMILY_HMMSEARCH        } from '../modules/family_hmmsearch'
include { PROFILE_PRESENCE_MATRIX } from '../modules/profile_presence_matrix'

// Per-ingroup-proteome protein_id -> proteome Short map, merged into one TSV so
// profile_to_matrix.py can attribute each family member to its source proteome.
process INGROUP_PROTEIN_MAP {
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
    ingroup_ch   // [meta, protein_fa]  — ingroup proteomes (families seeded from these)
    outgroup_ch  // [meta, protein_fa]  — outgroup proteomes (probed for presence too)
    config_csv   // path to analysis CSV

    main:
    all_proteomes_ch = ingroup_ch.mix(outgroup_ch)

    // Concatenate the ingroup once (order irrelevant) and build the protein->proteome map.
    // .first() turns the single collected file into a value channel reusable by many tasks.
    ingroup_concat = ingroup_ch.map { meta, fa -> fa }
                               .collectFile(name: 'ingroup_all.faa')
                               .first()

    protein_map = INGROUP_PROTEIN_MAP(ingroup_ch)
                      .collectFile(name: 'protein_to_proteome.tsv')
                      .first()

    // Cluster ingroup -> families; build a family HMM db (famsa + hmmbuild, >=N members).
    MMSEQS_FAMILY_CLUSTER(ingroup_concat)
    BUILD_FAMILY_PROFILES(MMSEQS_FAMILY_CLUSTER.out.cluster_tsv, ingroup_concat)

    // Scan every proteome (both groups) with the family HMM database.
    FAMILY_HMMSEARCH(all_proteomes_ch, BUILD_FAMILY_PROFILES.out.profiles.first())

    // Presence matrix + candidates (same contract as BUILD_PRESENCE_MATRIX).
    PROFILE_PRESENCE_MATRIX(
        FAMILY_HMMSEARCH.out.domtblout.map { meta, dom -> dom }.collect(),
        MMSEQS_FAMILY_CLUSTER.out.cluster_tsv,
        BUILD_FAMILY_PROFILES.out.families,
        protein_map,
        config_csv
    )

    emit:
    matrix     = PROFILE_PRESENCE_MATRIX.out.matrix
    candidates = PROFILE_PRESENCE_MATRIX.out.candidates
}
