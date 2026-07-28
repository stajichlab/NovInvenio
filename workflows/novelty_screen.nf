nextflow.enable.dsl=2

// NOVELTY_SCREEN — phase 2 of the two-phase targeted novelty pipeline (issue #27, see
// todo/novelty-discovery-screen.md). Takes the calibrated family HMMs + phase-1 candidate
// list from NOVELTY_DISCOVERY and re-searches them against two broader proteome sets to
// reclassify each candidate:
//
//   target_specific -- not found in NEAR_IN or BROAD_OUT
//   clade_specific  -- found in NEAR_IN, not found in BROAD_OUT
//   false_novelty   -- found in BROAD_OUT (removed from the screened candidate list)
//
// Reuses FAMILY_HMMSEARCH (same module NOVELTY_DISCOVERY uses against DISC_OUT) -- called
// once per proteome set (aliased NEAR_IN_HMMSEARCH / BROAD_OUT_HMMSEARCH so each call's
// .out is captured separately, matching main.nf's CLUSTER/LOSS_CLUSTER convention for
// invoking the same workflow twice in one scope) -- and TBLASTN/SUMMARIZE_TBLASTN (against
// BROAD_OUT genomes instead of DISC_OUT genomes). No new search machinery, just a new
// proteome/genome panel and a classification step (bin/novelty_screen.py).

include { FAMILY_HMMSEARCH as NEAR_IN_HMMSEARCH   } from '../modules/family_hmmsearch'
include { FAMILY_HMMSEARCH as BROAD_OUT_HMMSEARCH } from '../modules/family_hmmsearch'
include { TBLASTN_MAKEDB    } from '../modules/tblastn'
include { TBLASTN           } from '../modules/tblastn'
include { SUMMARIZE_TBLASTN } from '../modules/summarize_tblastn'

// ---------------------------------------------------------------------------
// Reclassify phase-1 candidates using NEAR_IN/BROAD_OUT hmmsearch evidence.
// ---------------------------------------------------------------------------
process NOVELTY_SCREEN_CLASSIFY {
    label 'low_cpu'
    tag "novelty_screen_classify"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(discovery_matrix)
    path(discovery_candidates)
    path(cluster_tsv)
    path(near_in_domtblouts)
    path(broad_out_domtblouts)
    path(family_thresholds)
    path(config_csv)
    val(default_family_evalue)
    val(min_coverage)

    output:
    path("screened_presence_matrix.tsv"), emit: matrix
    path("screened_candidates.txt"),      emit: candidates

    script:
    def near_in_arg = near_in_domtblouts ? "--near-in-domtblout ${near_in_domtblouts}" : ''
    def broad_out_arg = broad_out_domtblouts ? "--broad-out-domtblout ${broad_out_domtblouts}" : ''
    """
    novelty_screen.py \
        --discovery-matrix ${discovery_matrix} \
        --discovery-candidates ${discovery_candidates} \
        --cluster-tsv ${cluster_tsv} \
        ${near_in_arg} \
        ${broad_out_arg} \
        --family-thresholds ${family_thresholds} \
        --default-family-evalue ${default_family_evalue} \
        --min-coverage ${min_coverage} \
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
    near_in_ch           // [meta, protein_fa] — NEAR_IN proteomes (20-50 species)
    broad_out_ch         // [meta, protein_fa] — BROAD_OUT proteomes (20-100 species)
    broad_out_dna        // [meta, dna_fa] — BROAD_OUT genomes (for TBLASTN)
    calibrated_hmms      // path: family_profiles.hmm from NOVELTY_DISCOVERY
    family_thresholds    // path: family_thresholds.tsv from NOVELTY_DISCOVERY
    family_reps          // path: families_rep_seq.fasta from NOVELTY_DISCOVERY
    cluster_tsv          // path: families_cluster.tsv from NOVELTY_DISCOVERY
    discovery_matrix      // path: presence_matrix.tsv from NOVELTY_DISCOVERY
    discovery_candidates  // path: candidates.txt from NOVELTY_DISCOVERY
    config_csv            // path: analysis CSV

    main:
    // hmmsearch the calibrated family HMMs against NEAR_IN and BROAD_OUT separately (rather
    // than mixed) so each proteome set's domtblouts can be told apart without a join. Empty
    // proteome channels (a config with no NEAR_IN/BROAD_OUT rows) degrade gracefully: 0
    // proteomes to search means 0 hits found, so every phase-1 candidate stays
    // target_specific by default -- there is simply no broader-screen evidence to demote it.
    NEAR_IN_HMMSEARCH(near_in_ch, calibrated_hmms, 'screen_near_')
    near_in_domtblouts = NEAR_IN_HMMSEARCH.out.domtblout
        .map { meta, dom -> dom }.collect().ifEmpty([])

    BROAD_OUT_HMMSEARCH(broad_out_ch, calibrated_hmms, 'screen_broad_')
    broad_out_domtblouts = BROAD_OUT_HMMSEARCH.out.domtblout
        .map { meta, dom -> dom }.collect().ifEmpty([])

    // TBLASTN genomic validation against BROAD_OUT genomes (reporting-only, mirroring
    // NOVELTY_DISCOVERY's own DISC_OUT validation -- see bin/make_novelties.py's
    // --skip_tblastn_filter rationale: an absence-of-hit in the protein search does not
    // prove absence in the genome).
    genome_db_ch = TBLASTN_MAKEDB(broad_out_dna)
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
        near_in_domtblouts,
        broad_out_domtblouts,
        family_thresholds,
        config_csv,
        params.hmm_presence_evalue,
        params.hmm_presence_cov
    )

    emit:
    matrix          = NOVELTY_SCREEN_CLASSIFY.out.matrix
    candidates      = NOVELTY_SCREEN_CLASSIFY.out.candidates
    tblastn_summary = SUMMARIZE_TBLASTN.out.tsv
}
