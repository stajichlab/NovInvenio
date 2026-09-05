nextflow.enable.dsl=2

// NOVELTY_SCREEN — phase 2 of the two-phase targeted novelty pipeline (issue #27, see
// todo/novelty-discovery-screen.md). Takes the calibrated family HMMs + phase-1 candidate
// list from NOVELTY_DISCOVERY and re-searches them against two broader proteome sets to
// reclassify each candidate:
//
//   target_specific -- not found in NEAR_INGROUP or BROAD_OUTGROUP
//   clade_specific  -- found in NEAR_INGROUP, not found in BROAD_OUTGROUP
//   false_novelty   -- found in BROAD_OUTGROUP (removed from the screened candidate list)
//
// Reuses FAMILY_HMMSEARCH (same module NOVELTY_DISCOVERY uses against DISCOVERY_OUT) -- called
// once per proteome set (aliased NEAR_INGROUP_HMMSEARCH / BROAD_OUTGROUP_HMMSEARCH so each call's
// .out is captured separately, matching main.nf's CLUSTER/LOSS_CLUSTER convention for
// invoking the same workflow twice in one scope) -- and TBLASTN/SUMMARIZE_TBLASTN (against
// BROAD_OUTGROUP genomes instead of DISCOVERY_OUT genomes). No new search machinery, just a new
// proteome/genome panel and a classification step (bin/novelty_screen.py).

include { FAMILY_HMMSEARCH as NEAR_INGROUP_HMMSEARCH   } from '../modules/family_hmmsearch'
include { FAMILY_HMMSEARCH as BROAD_OUTGROUP_HMMSEARCH } from '../modules/family_hmmsearch'
include { TBLASTN_MAKEDB    } from '../modules/tblastn'
include { TBLASTN           } from '../modules/tblastn'
include { SUMMARIZE_TBLASTN } from '../modules/summarize_tblastn'

// Singleton screening (issue #52): reuse NOVELTY_DISCOVERY's extended singleton query
// (singletons + their own paralogs, see bin/extend_singleton_query.py) and search it
// against NEAR_INGROUP/BROAD_OUTGROUP too, with the same paralog-aware filtering phase 1
// applies -- closes the gap where singleton candidates got no phase-2 scrutiny at all
// and always defaulted to target_specific. Aliased per proteome set (near/broad) since a
// process can only be invoked once per execution path within one workflow scope.
include { SINGLETON_PHMMER_SEARCH  as SCREEN_SINGLETON_PHMMER_SEARCH_NEAR   } from './novelty_discovery'
include { SINGLETON_PHMMER_SEARCH  as SCREEN_SINGLETON_PHMMER_SEARCH_BROAD  } from './novelty_discovery'
include { SINGLETON_DIAMOND_SEARCH as SCREEN_SINGLETON_DIAMOND_SEARCH_NEAR  } from './novelty_discovery'
include { SINGLETON_DIAMOND_SEARCH as SCREEN_SINGLETON_DIAMOND_SEARCH_BROAD } from './novelty_discovery'
include { SINGLETON_BLAST_SEARCH   as SCREEN_SINGLETON_BLAST_SEARCH_NEAR    } from './novelty_discovery'
include { SINGLETON_BLAST_SEARCH   as SCREEN_SINGLETON_BLAST_SEARCH_BROAD   } from './novelty_discovery'
include { DIAMOND_MAKEDB as SCREEN_NEAR_DIAMOND_MAKEDB  } from '../modules/diamond'
include { DIAMOND_MAKEDB as SCREEN_BROAD_DIAMOND_MAKEDB } from '../modules/diamond'
include { BLAST_MAKEDB   as SCREEN_NEAR_BLAST_MAKEDB    } from '../modules/blast'
include { BLAST_MAKEDB   as SCREEN_BROAD_BLAST_MAKEDB   } from '../modules/blast'
include { PARSE_HITS as SCREEN_PARSE_HITS_NEAR  } from '../modules/parse_hits'
include { PARSE_HITS as SCREEN_PARSE_HITS_BROAD } from '../modules/parse_hits'

// ---------------------------------------------------------------------------
// Reclassify phase-1 candidates using NEAR_INGROUP/BROAD_OUTGROUP hmmsearch evidence.
// ---------------------------------------------------------------------------
process NOVELTY_SCREEN_CLASSIFY {
    label 'low_cpu'
    tag "novelty_screen_classify"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(discovery_matrix)
    path(discovery_candidates)
    path(cluster_tsv)
    path(near_ingroup_domtblouts)
    path(broad_outgroup_domtblouts)
    path(family_thresholds)
    path(config_csv)
    val(default_family_evalue)
    val(min_coverage)
    val(min_covered_residues)
    path(near_singleton_hits)
    path(broad_singleton_hits)
    path(paralog_cutoffs)
    val(singleton_evalue)
    val(paralog_competition_scope)

    output:
    path("screened_presence_matrix.tsv"), emit: matrix
    path("screened_candidates.txt"),      emit: candidates

    script:
    def near_ingroup_arg = near_ingroup_domtblouts ? "--near-in-domtblout ${near_ingroup_domtblouts}" : ''
    def broad_outgroup_arg = broad_outgroup_domtblouts ? "--broad-out-domtblout ${broad_outgroup_domtblouts}" : ''
    def near_singleton_arg = near_singleton_hits ? "--near-in-singleton-hits ${near_singleton_hits}" : ''
    def broad_singleton_arg = broad_singleton_hits ? "--broad-out-singleton-hits ${broad_singleton_hits}" : ''
    def paralog_arg = paralog_cutoffs ? "--paralog-cutoffs ${paralog_cutoffs}" : ''
    """
    novelty_screen.py \
        --discovery-matrix ${discovery_matrix} \
        --discovery-candidates ${discovery_candidates} \
        --cluster-tsv ${cluster_tsv} \
        ${near_ingroup_arg} \
        ${broad_outgroup_arg} \
        --family-thresholds ${family_thresholds} \
        --default-family-evalue ${default_family_evalue} \
        --min-coverage ${min_coverage} \
        --min-covered-residues ${min_covered_residues} \
        ${near_singleton_arg} \
        ${broad_singleton_arg} \
        ${paralog_arg} \
        --singleton-evalue ${singleton_evalue} \
        --paralog-competition-scope ${paralog_competition_scope} \
        --config ${config_csv} \
        --output-matrix screened_presence_matrix.tsv \
        --output-candidates screened_candidates.txt
    """
}

// ---------------------------------------------------------------------------
// Main workflow
// ---------------------------------------------------------------------------
workflow NOVELTY_SCREEN {
    take:
    near_ingroup_ch           // [meta, protein_fa] — NEAR_INGROUP proteomes (20-50 species)
    broad_outgroup_ch         // [meta, protein_fa] — BROAD_OUTGROUP proteomes (20-100 species)
    broad_outgroup_dna        // [meta, dna_fa] — BROAD_OUTGROUP genomes (for TBLASTN)
    calibrated_hmms      // path: family_profiles.hmm from NOVELTY_DISCOVERY
    family_thresholds    // path: family_thresholds.tsv from NOVELTY_DISCOVERY -- passed
                          //   through for process interface stability only; presence is
                          //   gated by the flat --default-family-evalue, not this file (see
                          //   lib/family_presence.py's module docstring)
    family_reps          // path: families_rep_seq.fasta from NOVELTY_DISCOVERY
    cluster_tsv          // path: families_cluster.tsv from NOVELTY_DISCOVERY
    discovery_matrix      // path: presence_matrix.tsv from NOVELTY_DISCOVERY
    discovery_candidates  // path: candidates.txt from NOVELTY_DISCOVERY
    singleton_query_fa    // path: singleton_query.fa from NOVELTY_DISCOVERY (issue #52) --
                           //   singletons + their own paralogs, same query phase 1 searched
    paralog_cutoffs       // path(s): paralog_cutoffs.tsv from NOVELTY_DISCOVERY's
                           //   DISCOVERY_TARGET self-vs-self search
    config_csv            // path: analysis CSV

    main:
    // hmmsearch the calibrated family HMMs against NEAR_INGROUP and BROAD_OUTGROUP separately (rather
    // than mixed) so each proteome set's domtblouts can be told apart without a join. Empty
    // proteome channels (a config with no NEAR_INGROUP/BROAD_OUTGROUP rows) degrade gracefully: 0
    // proteomes to search means 0 hits found, so every phase-1 candidate stays
    // target_specific by default -- there is simply no broader-screen evidence to demote it.
    NEAR_INGROUP_HMMSEARCH(near_ingroup_ch, calibrated_hmms, 'screen_near_')
    near_ingroup_domtblouts = NEAR_INGROUP_HMMSEARCH.out.domtblout
        .map { meta, dom -> dom }.collect().ifEmpty([])

    BROAD_OUTGROUP_HMMSEARCH(broad_outgroup_ch, calibrated_hmms, 'screen_broad_')
    broad_outgroup_domtblouts = BROAD_OUTGROUP_HMMSEARCH.out.domtblout
        .map { meta, dom -> dom }.collect().ifEmpty([])

    // Singleton screening (issue #52): the SAME extended singleton query phase 1 searched
    // (singletons + their own paralogs), searched against NEAR_INGROUP and BROAD_OUTGROUP
    // separately so bin/novelty_screen.py can tell the two sets' hits apart. Degrades
    // gracefully like the family HMM search above: empty proteome channels mean empty hit
    // lists, so singletons simply stay target_specific with no phase-2 evidence to demote them.
    if (params.run_tool == 'phmmer') {
        SCREEN_SINGLETON_PHMMER_SEARCH_NEAR(singleton_query_fa, near_ingroup_ch)
        near_singleton_raw_ch = SCREEN_SINGLETON_PHMMER_SEARCH_NEAR.out
        SCREEN_SINGLETON_PHMMER_SEARCH_BROAD(singleton_query_fa, broad_outgroup_ch)
        broad_singleton_raw_ch = SCREEN_SINGLETON_PHMMER_SEARCH_BROAD.out
    }
    else if (params.run_tool == 'diamond') {
        near_singleton_db_ch = SCREEN_NEAR_DIAMOND_MAKEDB(near_ingroup_ch)
        SCREEN_SINGLETON_DIAMOND_SEARCH_NEAR(singleton_query_fa, near_singleton_db_ch)
        near_singleton_raw_ch = SCREEN_SINGLETON_DIAMOND_SEARCH_NEAR.out
        broad_singleton_db_ch = SCREEN_BROAD_DIAMOND_MAKEDB(broad_outgroup_ch)
        SCREEN_SINGLETON_DIAMOND_SEARCH_BROAD(singleton_query_fa, broad_singleton_db_ch)
        broad_singleton_raw_ch = SCREEN_SINGLETON_DIAMOND_SEARCH_BROAD.out
    }
    else if (params.run_tool == 'blast') {
        near_singleton_db_ch = SCREEN_NEAR_BLAST_MAKEDB(near_ingroup_ch)
        SCREEN_SINGLETON_BLAST_SEARCH_NEAR(singleton_query_fa, near_singleton_db_ch)
        near_singleton_raw_ch = SCREEN_SINGLETON_BLAST_SEARCH_NEAR.out
        broad_singleton_db_ch = SCREEN_BROAD_BLAST_MAKEDB(broad_outgroup_ch)
        SCREEN_SINGLETON_BLAST_SEARCH_BROAD(singleton_query_fa, broad_singleton_db_ch)
        broad_singleton_raw_ch = SCREEN_SINGLETON_BLAST_SEARCH_BROAD.out
    }
    else {
        error "Unknown --run_tool '${params.run_tool}': choose phmmer, diamond, or blast"
    }
    SCREEN_PARSE_HITS_NEAR(near_singleton_raw_ch)
    near_singleton_hits_ch = SCREEN_PARSE_HITS_NEAR.out.map { meta_pair, tsv -> tsv }.collect().ifEmpty([])
    SCREEN_PARSE_HITS_BROAD(broad_singleton_raw_ch)
    broad_singleton_hits_ch = SCREEN_PARSE_HITS_BROAD.out.map { meta_pair, tsv -> tsv }.collect().ifEmpty([])

    // TBLASTN genomic validation against BROAD_OUTGROUP genomes (reporting-only, mirroring
    // NOVELTY_DISCOVERY's own DISCOVERY_OUT validation -- see bin/make_novelties.py's
    // --skip_tblastn_filter rationale: an absence-of-hit in the protein search does not
    // prove absence in the genome).
    genome_db_ch = TBLASTN_MAKEDB(broad_outgroup_dna)
    TBLASTN(genome_db_ch, family_reps)

    SUMMARIZE_TBLASTN(
        TBLASTN.out.tsv.map { meta, tsv -> tsv }.collect().ifEmpty([]),
        cluster_tsv,
        'screen_tblastn_summary.tsv'
    )

    NOVELTY_SCREEN_CLASSIFY(
        discovery_matrix,
        discovery_candidates,
        cluster_tsv,
        near_ingroup_domtblouts,
        broad_outgroup_domtblouts,
        family_thresholds,
        config_csv,
        params.hmm_presence_evalue,
        params.hmm_presence_cov,
        params.hmm_presence_min_residues,
        near_singleton_hits_ch,
        broad_singleton_hits_ch,
        paralog_cutoffs,
        params.evalue,
        params.paralog_competition_scope
    )

    emit:
    matrix          = NOVELTY_SCREEN_CLASSIFY.out.matrix
    candidates      = NOVELTY_SCREEN_CLASSIFY.out.candidates
    tblastn_summary = SUMMARIZE_TBLASTN.out.tsv
}
