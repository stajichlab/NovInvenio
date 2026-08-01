nextflow.enable.dsl=2

// NOVELTY_DISCOVERY — two-phase targeted novelty pipeline (Ticket #26).
//
// Takes the species roles defined in the config CSV (DISCOVERY_TARGET, DISCOVERY_OUT --
// old aliases TARGET/DISC_OUT still accepted, see lib/config_parser.py GROUP_ALIASES) and
// produces calibrated family HMMs for in-group specific families.
//
// Branching logic:
//   |DISCOVERY_TARGET| == 1  → pairwise search of the single target proteome against
//                    each DISCOVERY_OUT proteome.  No clustering, no HMMs.
//   |DISCOVERY_TARGET| >= 2  → mmseqs cluster the target proteomes into gene families;
//                    build family HMMs (famsa + hmmbuild) for multi-member
//                    families; pairwise search for singletons.  Merge both
//                    into one presence matrix.  Calibrate each surviving
//                    family HMM using the DISCOVERY_OUT panel as a negative
//                    control.  TBLASTN validation against DISCOVERY_OUT genomes.

include { DIAMOND_MAKEDB        } from '../modules/diamond'
include { BLAST_MAKEDB          } from '../modules/blast'
include { PARSE_HITS            } from '../modules/parse_hits'
include { MMSEQS_FAMILY_CLUSTER } from '../modules/mmseqs_family_cluster'
include { BUILD_FAMILY_PROFILES } from '../modules/build_family_profiles'
include { FAMILY_HMMSEARCH      } from '../modules/family_hmmsearch'
include { TBLASTN_MAKEDB        } from '../modules/tblastn'
include { TBLASTN               } from '../modules/tblastn'
include { SUMMARIZE_TBLASTN     } from '../modules/summarize_tblastn'
include { PHMMER_SELF           } from '../modules/self_search'
include { DIAMOND_SELF          } from '../modules/self_search'
include { BLAST_SELF            } from '../modules/self_search'
include { PARSE_SELF_HITS       } from '../modules/parse_self_hits'

// ---------------------------------------------------------------------------
// Per-seed-proteome protein_id → proteome Short map, merged into one TSV so
// novelty_presence_matrix.py can attribute each family member / singleton to
// its source proteome.
// ---------------------------------------------------------------------------
process SEED_PROTEIN_MAP {
    label 'low_cpu'
    tag "${meta.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

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

// ---------------------------------------------------------------------------
// Extract singleton sequences from the mmseqs cluster results.
// Singletons = clusters with exactly one member (rep == member).
// ---------------------------------------------------------------------------
process EXTRACT_SINGLETONS {
    label 'low_cpu'
    tag "extract_singletons"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(cluster_tsv)
    path(seed_fa)

    output:
    path("singletons.fa"), emit: singletons_fa

    script:
    """
    extract_singletons.py \
        --cluster-tsv ${cluster_tsv} \
        --fasta ${seed_fa} \
        --output singletons.fa
    """
}

// ---------------------------------------------------------------------------
// Extend the singleton query with each singleton's own within-genome paralog (issue
// #52), so the singleton search's paralog-competition filter has the paralog's own
// hits to compare against. See bin/extend_singleton_query.py.
// ---------------------------------------------------------------------------
process EXTEND_SINGLETON_QUERY {
    label 'low_cpu'
    tag "extend_singleton_query"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(singletons_fa)
    path(seed_fa)
    path(paralog_cutoffs)

    output:
    path("singleton_query.fa"), emit: query_fa

    script:
    """
    extend_singleton_query.py \
        --singletons-fa ${singletons_fa} \
        --seed-fasta ${seed_fa} \
        --paralog-cutoffs ${paralog_cutoffs} \
        --output singleton_query.fa
    """
}

// ---------------------------------------------------------------------------
// Singleton pairwise search ("hybrid approach", todo/novelty-discovery-screen.md).
// mmseqs clustering can miss real orthologs when annotated gene models differ
// substantially in length (e.g. a large N/C-terminal extension unique to one
// species' gene model) -- its bidirectional coverage requirement
// (-c 0.8 --cov-mode 0) structurally excludes such pairs regardless of identity
// in the aligned region. A protein that ends up an mmseqs "singleton" (no
// cluster partners) is otherwise invisible to the discovery presence matrix
// beyond its own source proteome. This searches every singleton, once, against
// every TARGET+DISCOVERY_OUT proteome -- presence in other TARGET genomes
// (mmseqs-clustering-independent) and absence in DISCOVERY_OUT.
//
// No storeDir, unlike modules/phmmer.nf|diamond.nf|blast.nf's per-species search
// processes: singletons.fa is DERIVED (its content changes whenever the TARGET
// set or clustering parameters change) rather than a stable per-species proteome,
// so filename-only caching would silently reuse a stale result for a different
// target set under the same --project name.
// ---------------------------------------------------------------------------
process SINGLETON_PHMMER_SEARCH {
    label 'med_cpu'
    tag "singletons_vs_${meta_t.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(singletons_fa)
    tuple val(meta_t), path(target_fa)

    output:
    tuple val(meta_pair), path("singletons_vs_${meta_t.id}.phmmer.tblout")

    script:
    meta_pair = [query_id: 'singletons', target_id: meta_t.id, tool: 'phmmer']
    """
    if [ -s ${singletons_fa} ]; then
        phmmer \
            --cpu ${task.cpus} \
            --tblout singletons_vs_${meta_t.id}.phmmer.tblout \
            --noali \
            ${singletons_fa} ${target_fa} \
            > /dev/null
    else
        touch singletons_vs_${meta_t.id}.phmmer.tblout
    fi
    """
}

process SINGLETON_DIAMOND_SEARCH {
    label 'high_cpu'
    tag "singletons_vs_${meta_t.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(singletons_fa)
    tuple val(meta_t), path(target_db)

    output:
    tuple val(meta_pair), path("singletons_vs_${meta_t.id}.diamond.tsv")

    script:
    meta_pair = [query_id: 'singletons', target_id: meta_t.id, tool: 'diamond']
    """
    if [ -s ${singletons_fa} ]; then
        diamond blastp \
            --query ${singletons_fa} \
            --db ${target_db.baseName} \
            --outfmt 6 qseqid sseqid evalue bitscore \
            --evalue ${params.parse_evalue} \
            --threads ${task.cpus} \
            --quiet \
            --out singletons_vs_${meta_t.id}.diamond.tsv
    else
        touch singletons_vs_${meta_t.id}.diamond.tsv
    fi
    """
}

process SINGLETON_BLAST_SEARCH {
    label 'high_cpu'
    tag "singletons_vs_${meta_t.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(singletons_fa)
    tuple val(meta_t), path(target_db)

    output:
    tuple val(meta_pair), path("singletons_vs_${meta_t.id}.blast.tsv")

    script:
    meta_pair = [query_id: 'singletons', target_id: meta_t.id, tool: 'blast']
    """
    if [ -s ${singletons_fa} ]; then
        blastp \
            -query ${singletons_fa} \
            -db ${meta_t.id}.blast_db \
            -outfmt "6 qseqid sseqid evalue bitscore" \
            -evalue ${params.parse_evalue} \
            -num_threads ${task.cpus} \
            -out singletons_vs_${meta_t.id}.blast.tsv
    else
        touch singletons_vs_${meta_t.id}.blast.tsv
    fi
    """
}

// ---------------------------------------------------------------------------
// Calibrate family HMM thresholds using the DISCOVERY_OUT panel as a negative
// control.  For each surviving family HMM, we search it against the DISCOVERY_OUT
// proteomes and record the highest-scoring (lowest E-value) hit.  The
// per-family threshold is set to that E-value, so that the family is only
// called "present" in a downstream proteome if the hit is at least as strong
// as the best false-positive.  If a family HMM has no hit in any DISCOVERY_OUT
// proteome, the threshold falls back to the global E-value parameter.
// ---------------------------------------------------------------------------
process CALIBRATE_FAMILY_HMMS {
    label 'low_cpu'
    tag "calibrate_hmms"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(domtblouts)
    path(families_tsv)
    val(default_evalue)

    output:
    path("family_thresholds.tsv"), emit: thresholds

    script:
    """
    calibrate_family_hmms.py \
        --domtblout ${domtblouts} \
        --families ${families_tsv} \
        --default-evalue ${default_evalue} \
        --output family_thresholds.tsv
    """
}

// ---------------------------------------------------------------------------
// Build the novelty discovery presence matrix from two evidence sources:
//   1. Family HMM search (multi-member families)
//   2. Singleton pairwise search
// Merge both into one presence matrix and candidate list.
// ---------------------------------------------------------------------------
process NOVELTY_PRESENCE_MATRIX {
    label 'low_cpu'
    tag "novelty_matrix"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(family_domtblouts)
    path(singleton_hits)
    path(cluster_tsv)
    path(protein_map)
    path(config_csv)
    path(family_thresholds)
    val(default_family_evalue)
    val(min_coverage)
    val(target_min_frac)
    val(disc_out_max_frac)
    val(singleton_evalue)
    path(paralog_cutoffs)
    val(paralog_competition_scope)

    output:
    path("presence_matrix.tsv"), emit: matrix
    path("candidates.txt"),      emit: candidates
    path("presence_matrix.evalues.tsv"), emit: evalues

    script:
    // singleton_hits/paralog_cutoffs are [] when there were no singletons / no paralogs
    // detected at all -- omit the flag rather than pass an empty Groovy list literal
    // ("[]") on the command line.
    def singleton_arg = singleton_hits ? "--singleton-hits ${singleton_hits}" : ''
    def paralog_arg = paralog_cutoffs ? "--paralog-cutoffs ${paralog_cutoffs}" : ''
    """
    novelty_presence_matrix.py \
        --family-domtblout ${family_domtblouts} \
        ${singleton_arg} \
        --cluster-tsv ${cluster_tsv} \
        --protein-map ${protein_map} \
        --config ${config_csv} \
        --family-thresholds ${family_thresholds} \
        --default-family-evalue ${default_family_evalue} \
        --min-coverage ${min_coverage} \
        --target-min-frac ${target_min_frac} \
        --disc-out-max-frac ${disc_out_max_frac} \
        --singleton-evalue ${singleton_evalue} \
        ${paralog_arg} \
        --paralog-competition-scope ${paralog_competition_scope} \
        --output-matrix presence_matrix.tsv \
        --output-candidates candidates.txt \
        --output-evalues presence_matrix.evalues.tsv
    """
}

// ---------------------------------------------------------------------------
// Main workflow
// ---------------------------------------------------------------------------
workflow NOVELTY_DISCOVERY {
    take:
    target_ch      // [meta, protein_fa] — DISCOVERY_TARGET genomes (2-3)
    disc_out_ch    // [meta, protein_fa] — DISCOVERY_OUT proteomes
    disc_out_dna   // [meta, dna_fa] — DISCOVERY_OUT genomes (for TBLASTN)
    config_csv     // path to analysis CSV

    main:
    all_proteomes_ch = target_ch.mix(disc_out_ch)

    // Build the protein→proteome map from target proteomes.
    protein_map = SEED_PROTEIN_MAP(target_ch)
                      .collectFile(name: "protein_to_proteome.tsv")
                      .first()

    // Concatenate target proteomes once.
    seed_concat = target_ch.map { meta, fa -> fa }
                         .collectFile(name: "seed_all.faa")
                         .first()

    // Cluster target proteomes into gene families.
    MMSEQS_FAMILY_CLUSTER(seed_concat, '')

    // Build family HMMs for multi-member families.
    BUILD_FAMILY_PROFILES(MMSEQS_FAMILY_CLUSTER.out.cluster_tsv, seed_concat, '')
    family_profiles = BUILD_FAMILY_PROFILES.out.profiles.first()

    // Search family HMMs against ALL proteomes (target + DISCOVERY_OUT).
    FAMILY_HMMSEARCH(all_proteomes_ch, family_profiles, '')

    // Extract singleton sequences from the cluster results.
    EXTRACT_SINGLETONS(MMSEQS_FAMILY_CLUSTER.out.cluster_tsv, seed_concat)
    singletons_fa = EXTRACT_SINGLETONS.out.singletons_fa.first()

    // Self-vs-self search on DISCOVERY_TARGET proteomes (issue #52) -- captures each
    // target protein's best within-proteome paralog, the same rank-2 calibration
    // workflows/search.nf does for the classic pairwise pathway. Needed so the singleton
    // search below can apply the same paralog-cutoff + paralog-competition filters
    // bin/build_presence_matrix.py uses, instead of a flat --singleton-evalue threshold
    // with no protection against cross-reactivity with a conserved paralog (the
    // NCU08332/HEX-1-vs-eIF5A false positive found validating job 26997324).
    if (params.run_tool == 'phmmer') {
        raw_self_ch = PHMMER_SELF(target_ch)
    }
    else if (params.run_tool == 'diamond') {
        raw_self_ch = DIAMOND_SELF(target_ch)
    }
    else if (params.run_tool == 'blast') {
        raw_self_ch = BLAST_SELF(target_ch)
    }
    else {
        error "Unknown --run_tool '${params.run_tool}': choose phmmer, diamond, or blast"
    }
    PARSE_SELF_HITS(raw_self_ch)
    paralog_cutoffs_ch = PARSE_SELF_HITS.out.tsv.map { meta, tsv -> tsv }.collect().ifEmpty([])

    // Extend the singleton query with each singleton's own paralog (issue #52) before
    // searching, so the paralog-competition filter in bin/novelty_presence_matrix.py has
    // the paralog's own hits to compare against.
    EXTEND_SINGLETON_QUERY(singletons_fa, seed_concat, paralog_cutoffs_ch)
    singleton_query_fa = EXTEND_SINGLETON_QUERY.out.query_fa.first()

    // Singleton pairwise search ("hybrid approach") against every TARGET+DISCOVERY_OUT
    // proteome -- see the SINGLETON_*_SEARCH process docs above for why mmseqs clustering
    // alone isn't enough.
    if (params.run_tool == 'phmmer') {
        SINGLETON_PHMMER_SEARCH(singleton_query_fa, all_proteomes_ch)
        singleton_raw_ch = SINGLETON_PHMMER_SEARCH.out
    }
    else if (params.run_tool == 'diamond') {
        singleton_db_ch = DIAMOND_MAKEDB(all_proteomes_ch)
        SINGLETON_DIAMOND_SEARCH(singleton_query_fa, singleton_db_ch)
        singleton_raw_ch = SINGLETON_DIAMOND_SEARCH.out
    }
    else if (params.run_tool == 'blast') {
        singleton_db_ch = BLAST_MAKEDB(all_proteomes_ch)
        SINGLETON_BLAST_SEARCH(singleton_query_fa, singleton_db_ch)
        singleton_raw_ch = SINGLETON_BLAST_SEARCH.out
    }
    else {
        error "Unknown --run_tool '${params.run_tool}': choose phmmer, diamond, or blast"
    }
    PARSE_HITS(singleton_raw_ch)
    singleton_hits_ch = PARSE_HITS.out.map { meta_pair, tsv -> tsv }.collect().ifEmpty([])

    // Calibrate family HMM thresholds using DISCOVERY_OUT as negative control. Must be
    // restricted to DISCOVERY_OUT proteomes' domtblouts only: FAMILY_HMMSEARCH returns
    // [meta, domtblout] for ALL proteomes (target + DISCOVERY_OUT), and a family's near-perfect
    // self-hit against its own DISCOVERY_TARGET source protein would otherwise poison the "negative
    // control" with an all-but-impossible-to-beat E-value, making every family's threshold
    // effectively unmatchable in phase 2 (found via issue #29's integration test).
    disc_out_domtblouts = FAMILY_HMMSEARCH.out.domtblout
        .filter { meta, dom -> meta.group == 'DISCOVERY_OUT' }
        .map { meta, dom -> dom }
        .collect()
        .ifEmpty([])
    CALIBRATE_FAMILY_HMMS(
        disc_out_domtblouts,
        BUILD_FAMILY_PROFILES.out.families.first(),
        params.hmm_presence_evalue
    )

    // Build the combined presence matrix.
    NOVELTY_PRESENCE_MATRIX(
        FAMILY_HMMSEARCH.out.domtblout.map { meta, dom -> dom }.collect(),
        singleton_hits_ch,
        MMSEQS_FAMILY_CLUSTER.out.cluster_tsv,
        protein_map,
        config_csv,
        CALIBRATE_FAMILY_HMMS.out.thresholds,
        params.hmm_presence_evalue,
        params.hmm_presence_cov,
        params.ingroup_min_frac,
        0.0,  // disc_out_max_frac: strictly absent from DISCOVERY_OUT
        params.evalue,  // singleton_evalue: fallback when a singleton has no detected paralog
        paralog_cutoffs_ch,  // issue #52: DISCOVERY_TARGET self-search paralog cutoffs
        params.paralog_competition_scope
    )

    // TBLASTN validation against DISCOVERY_OUT genomes — query is the family representative
    // protein sequences, not the HMM database.
    family_reps_fa = MMSEQS_FAMILY_CLUSTER.out.representatives.first()
    genome_db_ch = TBLASTN_MAKEDB(disc_out_dna)
    TBLASTN(genome_db_ch, family_reps_fa)

    // Summarize TBLASTN hits.
    SUMMARIZE_TBLASTN(
        TBLASTN.out.tsv.map { meta, tsv -> tsv }.collect(),
        MMSEQS_FAMILY_CLUSTER.out.cluster_tsv,
        'tblastn_summary.tsv'
    )

    emit:
    matrix             = NOVELTY_PRESENCE_MATRIX.out.matrix
    candidates         = NOVELTY_PRESENCE_MATRIX.out.candidates
    // Report-only hit e-value evidence (issue #44). Phase 1-scoped: DISCOVERY_TARGET and
    // DISCOVERY_OUT columns are populated; NEAR_INGROUP/BROAD_OUTGROUP columns stay empty
    // since novelty_screen.py's phase-2 HMM search only tracks presence, not e-value (a
    // follow-up if that evidence is wanted too). Keyed by protein_id/source_proteome, not
    // row order, so NOVELTY_SCREEN carrying rows forward doesn't need to touch this.
    evalues            = NOVELTY_PRESENCE_MATRIX.out.evalues
    summary            = SUMMARIZE_TBLASTN.out.tsv
    // Family clustering + reps + the concatenated seed group, for the family-as-cluster
    // path (ADR-0002 Q7): main.nf feeds these to PROFILE_CANDIDATE_CLUSTERS instead of
    // re-clustering the candidates (mirrors PROFILE_SEARCH's emit contract).
    family_cluster_tsv = MMSEQS_FAMILY_CLUSTER.out.cluster_tsv
    family_reps        = family_reps_fa
    seed_concat        = seed_concat
    calibrated_hmms    = family_profiles
    // Per-family E-value thresholds (rep_id -> threshold), calibrated against DISCOVERY_OUT as a
    // negative control -- novelty_screen (issue #27) reuses these same thresholds when
    // calling family presence against NEAR_INGROUP/BROAD_OUTGROUP, so a family's definition of
    // "present" stays consistent across both phases.
    family_thresholds  = CALIBRATE_FAMILY_HMMS.out.thresholds
    // issue #52: the extended singleton query (singletons + their own paralogs) and the
    // paralog cutoffs used to filter phase 1's singleton search -- NOVELTY_SCREEN reuses
    // both to search+filter the SAME singletons against NEAR_INGROUP/BROAD_OUTGROUP with
    // identical paralog-aware logic, closing the "singletons always default to
    // target_specific" gap.
    singleton_query_fa = singleton_query_fa
    paralog_cutoffs    = paralog_cutoffs_ch
}
