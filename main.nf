#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { SEARCH   } from './workflows/search'
include { PROFILE_SEARCH } from './workflows/profile_search'
include { PROFILE_SEARCH as PROFILE_LOSS_SEARCH } from './workflows/profile_search'
include { PROFILE_CANDIDATE_CLUSTERS } from './modules/profile_candidate_clusters'
include { PROFILE_CANDIDATE_CLUSTERS as LOSS_PROFILE_CANDIDATE_CLUSTERS } from './modules/profile_candidate_clusters'
include { LOSS_SEARCH } from './workflows/loss_search'
include { CLUSTER  } from './workflows/cluster'
include { CLUSTER  as LOSS_CLUSTER  } from './workflows/cluster'
include { VALIDATE } from './workflows/validate'
include { VALIDATE as LOSS_VALIDATE } from './workflows/validate'
include { ANNOTATE } from './workflows/annotate'
include { ANNOTATE as LOSS_ANNOTATE } from './workflows/annotate'
include { SUMMARIZE } from './workflows/summarize'
include { REPORT   } from './workflows/report'
include { NOVELTY_DISCOVERY } from './workflows/novelty_discovery'
include { NOVELTY_SCREEN } from './workflows/novelty_screen'
include { EMPTY_LOSS_STUB } from './modules/empty_loss_stub'

// Original novelty_discovery/novelty_screen GROUP labels (issues #24-#29), renamed for
// clarity (todo/rename-novelty-discovery-group-labels.md) -- still accepted in a config
// CSV's GROUP column and normalized to the canonical spelling here, mirroring
// lib/config_parser.py's GROUP_ALIASES so both the Python and Nextflow sides agree.
def normalizeGroup(String group) {
    def aliases = [
        'TARGET':    'DISCOVERY_TARGET',
        'DISC_OUT':  'DISCOVERY_OUT',
        'NEAR_IN':   'NEAR_INGROUP',
        'BROAD_OUT': 'BROAD_OUTGROUP',
    ]
    return aliases[group] ?: group
}

// Resolve a FASTA basename against data_dir, checking the flat layout first
// then the listed subdirectories (so configs that reference bare basenames
// still find files under data_dir/pep/ and data_dir/dna/).
def resolve_fa(String basename, List<String> subdirs) {
    if (!basename) return []
    def candidates = ([''] + subdirs).collect { sub ->
        file(sub ? "${params.data_dir}/${sub}/${basename}" : "${params.data_dir}/${basename}")
    }
    def hit = candidates.find { it.exists() }
    if (!hit) error "Cannot locate FASTA '${basename}' under ${params.data_dir} (also tried subdirs: ${subdirs.join(', ')})"
    return hit
}

workflow {
    if (!params.config)   error "ERROR: --config <analysis_csv> is required"
    if (!params.data_dir) error "ERROR: --data_dir <fasta_directory> is required"
    if (!file(params.config).exists())   error "ERROR: --config file not found: ${params.config}"
    if (!file(params.data_dir).isDirectory()) error "ERROR: --data_dir is not a directory: ${params.data_dir}"
    if (params.run_tool !in ['phmmer', 'diamond', 'blast']) error "ERROR: --run_tool must be phmmer, diamond, or blast (got: ${params.run_tool})"
    if (params.cluster_tool !in ['pairwise', 'mmseqs', 'novelty_discovery']) error "ERROR: --cluster_tool must be pairwise, mmseqs, or novelty_discovery (got: ${params.cluster_tool})"

    // Resolve DB paths to absolute at launch time and pass them as val inputs —
    // params mutations do not reliably propagate into process script closures.
    def pfam_abs  = params.pfam_hmm       ? file(params.pfam_hmm).toAbsolutePath().toString()       : ''
    def sprot_abs = params.swissprot_dmnd ? file(params.swissprot_dmnd).toAbsolutePath().toString() : ''
    def morgs_abs = params.modelorgs_config
        ? file(params.modelorgs_config).toAbsolutePath().toString()
        : ''

    // Parse the analysis description CSV into a channel of [meta, protein_fa, dna_fa]
    samples_ch = Channel
        .fromPath(params.config)
        .splitCsv(header: true)
        .map { row ->
            def group = normalizeGroup(row.GROUP?.trim())
            def meta = [
                id:      row.Short,
                group:   group,
                species: row.Species,
                strain:  row.Strain ?: '',
                taxon:   row.TaxonGroup
            ]
            def protein_fa = resolve_fa(row.Protein, ['pep', 'proteins'])
            def dna_fa     = resolve_fa(row.DNA,     ['dna', 'genome', 'scaffolds'])
            [ meta, protein_fa, dna_fa ]
        }

    ingroup_prot_ch   = samples_ch.filter { meta, prot, dna -> meta.group == 'IN' }
                                   .map    { meta, prot, dna -> [ meta, prot ] }
    outgroup_prot_ch  = samples_ch.filter { meta, prot, dna -> meta.group == 'OUT' }
                                   .map    { meta, prot, dna -> [ meta, prot ] }
    outgroup_dna_ch   = samples_ch.filter { meta, prot, dna -> meta.group == 'OUT' && dna }
                                   .map    { meta, prot, dna -> [ meta, dna ] }
    ingroup_dna_ch    = samples_ch.filter { meta, prot, dna -> meta.group == 'IN' && dna }
                                   .map    { meta, prot, dna -> [ meta, dna ] }

    // Channels for the two-phase novelty_discovery / novelty_screen workflow.
    target_prot_ch    = samples_ch.filter { meta, prot, dna -> meta.group == 'DISCOVERY_TARGET' }
                                   .map    { meta, prot, dna -> [ meta, prot ] }
    disc_out_prot_ch  = samples_ch.filter { meta, prot, dna -> meta.group == 'DISCOVERY_OUT' }
                                   .map    { meta, prot, dna -> [ meta, prot ] }
    disc_out_dna_ch   = samples_ch.filter { meta, prot, dna -> meta.group == 'DISCOVERY_OUT' && dna }
                                   .map    { meta, prot, dna -> [ meta, dna ] }
    near_in_prot_ch   = samples_ch.filter { meta, prot, dna -> meta.group == 'NEAR_INGROUP' }
                                   .map    { meta, prot, dna -> [ meta, prot ] }
    broad_out_prot_ch = samples_ch.filter { meta, prot, dna -> meta.group == 'BROAD_OUTGROUP' }
                                   .map    { meta, prot, dna -> [ meta, prot ] }
    broad_out_dna_ch  = samples_ch.filter { meta, prot, dna -> meta.group == 'BROAD_OUTGROUP' && dna }
                                   .map    { meta, prot, dna -> [ meta, dna ] }

    // Novelty-direction presence matrix + candidates. --cluster_tool selects the producer:
    //   pairwise (default) — the O(N^2) phmmer/diamond/blast SEARCH workflow.
    //   mmseqs             — the scalable family-profile pathway (ADR-0002). Both emit the
    //                        same matrix/candidates contract, so everything below is shared.
    if (params.cluster_tool == 'mmseqs') {
        // Novelty direction: seed families from the ingroup (query-group IN), absent from
        // every outgroup proteome (other_max_frac 0.0).
        PROFILE_SEARCH(ingroup_prot_ch, outgroup_prot_ch, file(params.config),
                       'IN', params.ingroup_min_frac, 0.0, '')
        novelty_matrix     = PROFILE_SEARCH.out.matrix
        novelty_candidates = PROFILE_SEARCH.out.candidates

        // Family-as-cluster (ADR-0002 Q7): reuse the profile pathway's gene families
        // (restricted to candidate-containing ones) instead of re-clustering candidates.
        PROFILE_CANDIDATE_CLUSTERS(
            novelty_candidates,
            PROFILE_SEARCH.out.family_cluster_tsv,
            PROFILE_SEARCH.out.family_reps,
            PROFILE_SEARCH.out.seed_concat,
            'candidates.fa',
            ''
        )
        cand_fa          = PROFILE_CANDIDATE_CLUSTERS.out.candidates_fa
        cand_reps        = PROFILE_CANDIDATE_CLUSTERS.out.representatives
        cand_cluster_tsv = PROFILE_CANDIDATE_CLUSTERS.out.cluster_tsv
    }
    else if (params.cluster_tool == 'novelty_discovery') {
        // Two-phase targeted novelty pipeline (see todo/novelty-discovery-screen.md).
        NOVELTY_DISCOVERY(target_prot_ch, disc_out_prot_ch, disc_out_dna_ch, file(params.config))

        // Phase 2 (issue #27): reclassify phase-1 candidates against NEAR_INGROUP
        // (clade-mates) and BROAD_OUTGROUP (distant lineages). Reuses the FULL calibrated
        // family HMM db and
        // family clustering (not narrowed to phase-1 candidates) — simpler wiring than
        // re-clustering twice, at the cost of a little extra hmmsearch/TBLASTN work on
        // families phase 1 already rejected.
        NOVELTY_SCREEN(
            near_in_prot_ch,
            broad_out_prot_ch,
            broad_out_dna_ch,
            NOVELTY_DISCOVERY.out.calibrated_hmms,
            NOVELTY_DISCOVERY.out.family_thresholds,
            NOVELTY_DISCOVERY.out.family_reps,
            NOVELTY_DISCOVERY.out.family_cluster_tsv,
            NOVELTY_DISCOVERY.out.matrix,
            NOVELTY_DISCOVERY.out.candidates,
            file(params.config)
        )
        novelty_matrix     = NOVELTY_SCREEN.out.matrix
        novelty_candidates = NOVELTY_SCREEN.out.candidates

        // Family-as-cluster (ADR-0002 Q7): reuse NOVELTY_DISCOVERY's own gene families
        // (restricted to candidate-containing ones) instead of re-clustering candidates.
        // Uses the *screened* candidate list (false_novelty already removed) so annotation
        // only runs on the surviving novelty candidates, per todo/novelty-discovery-screen.md
        // "Post-Screen: Annotation" — annotation is expensive, no reason to spend it on
        // families the screen phase already ruled out.
        PROFILE_CANDIDATE_CLUSTERS(
            novelty_candidates,
            NOVELTY_DISCOVERY.out.family_cluster_tsv,
            NOVELTY_DISCOVERY.out.family_reps,
            NOVELTY_DISCOVERY.out.seed_concat,
            'candidates.fa',
            ''
        )
        cand_fa          = PROFILE_CANDIDATE_CLUSTERS.out.candidates_fa
        cand_reps        = PROFILE_CANDIDATE_CLUSTERS.out.representatives
        cand_cluster_tsv = PROFILE_CANDIDATE_CLUSTERS.out.cluster_tsv

        // NOVELTY_DISCOVERY already ran its own TBLASTN vs DISCOVERY_OUT genomes and
        // summarized it (SUMMARIZE_TBLASTN) — the generic VALIDATE workflow below would be
        // redundant (and its outgroup_dna_ch is empty for DISCOVERY_TARGET/DISCOVERY_OUT
        // configs anyway), so this branch's TBLASTN summary is carried straight through to
        // REPORT/SUMMARIZE. NOVELTY_SCREEN.out.tblastn_summary (vs BROAD_OUTGROUP genomes)
        // is published separately
        // (screen_tblastn_summary.tsv) but not yet wired into the report -- which TBLASTN
        // evidence the final report surfaces is a report-rendering decision left to #28
        // alongside the novelty_category column.
        novelty_tblastn_summary = NOVELTY_DISCOVERY.out.summary
    }
    else {
        SEARCH(ingroup_prot_ch, outgroup_prot_ch, file(params.config))
        novelty_matrix     = SEARCH.out.matrix
        novelty_candidates = SEARCH.out.candidates

        CLUSTER(novelty_candidates, ingroup_prot_ch, file(params.config), 'candidates.fa', 'clusters')
        cand_fa          = CLUSTER.out.candidates_fa
        cand_reps        = CLUSTER.out.representatives
        cand_cluster_tsv = CLUSTER.out.cluster_tsv
    }

    // novelty_discovery already produced its own TBLASTN summary (see above); the other two
    // cluster_tool paths still need the generic VALIDATE (TBLASTN vs the OUT proteomes' DNA).
    if (params.cluster_tool != 'novelty_discovery') {
        VALIDATE(cand_reps, outgroup_dna_ch, cand_cluster_tsv, 'tblastn_summary.tsv')
        novelty_tblastn_summary = VALIDATE.out.summary
    }

    ANNOTATE(cand_fa, novelty_matrix, pfam_abs, sprot_abs, morgs_abs, '')

    SUMMARIZE(ANNOTATE.out.annotated_matrix, novelty_tblastn_summary, cand_cluster_tsv, file(params.config))

    // Loss direction — candidate lineage-specific gene losses (present in the outgroup,
    // absent from the ingroup). --cluster_tool selects the producer, mirroring the novelty
    // direction: pairwise LOSS_SEARCH, or the outgroup-seeded family-profile mirror
    // (ADR-0002 / issue #12). Both emit the same loss matrix/candidates contract.
    if (params.cluster_tool == 'mmseqs') {
        // Seed families from the OUTGROUP (query-group OUT): conserved in >= outgroup_min_frac
        // of the outgroup and present in <= loss_ingroup_max_frac of the ingroup.
        PROFILE_LOSS_SEARCH(outgroup_prot_ch, ingroup_prot_ch, file(params.config),
                            'OUT', params.outgroup_min_frac, params.loss_ingroup_max_frac, 'loss_')
        loss_matrix = PROFILE_LOSS_SEARCH.out.matrix

        LOSS_PROFILE_CANDIDATE_CLUSTERS(
            PROFILE_LOSS_SEARCH.out.candidates,
            PROFILE_LOSS_SEARCH.out.family_cluster_tsv,
            PROFILE_LOSS_SEARCH.out.family_reps,
            PROFILE_LOSS_SEARCH.out.seed_concat,
            'loss_candidates.fa',
            'loss_'
        )
        loss_cand_fa          = LOSS_PROFILE_CANDIDATE_CLUSTERS.out.candidates_fa
        loss_cand_reps        = LOSS_PROFILE_CANDIDATE_CLUSTERS.out.representatives
        loss_cand_cluster_tsv = LOSS_PROFILE_CANDIDATE_CLUSTERS.out.cluster_tsv
    }
    else if (params.cluster_tool == 'novelty_discovery') {
        // Loss analysis is an explicitly deferred future extension for the two-phase
        // novelty_discovery/novelty_screen plan (todo/novelty-discovery-screen.md); a
        // DISCOVERY_TARGET/DISCOVERY_OUT config has no IN/OUT rows, so LOSS_SEARCH would only ever see
        // empty channels. REPORT's COLLATE_REPORTS still needs a (zero-row) losses.html
        // to assemble view/<project>/, so stub the three loss artifacts instead.
        EMPTY_LOSS_STUB(ANNOTATE.out.annotated_matrix)
        loss_annotated_matrix   = EMPTY_LOSS_STUB.out.matrix
        loss_tblastn_summary    = EMPTY_LOSS_STUB.out.tblastn_summary
        loss_cand_cluster_tsv   = EMPTY_LOSS_STUB.out.cluster_tsv
    }
    else {
        // See workflows/loss_search.nf for why this needs its own search direction.
        LOSS_SEARCH(ingroup_prot_ch, outgroup_prot_ch, file(params.config))
        loss_matrix = LOSS_SEARCH.out.matrix

        LOSS_CLUSTER(LOSS_SEARCH.out.candidates, outgroup_prot_ch, file(params.config), 'loss_candidates.fa', 'loss_clusters')
        loss_cand_fa          = LOSS_CLUSTER.out.candidates_fa
        loss_cand_reps        = LOSS_CLUSTER.out.representatives
        loss_cand_cluster_tsv = LOSS_CLUSTER.out.cluster_tsv
    }

    if (params.cluster_tool != 'novelty_discovery') {
        LOSS_VALIDATE(loss_cand_reps, ingroup_dna_ch, loss_cand_cluster_tsv, 'loss_tblastn_summary.tsv')
        LOSS_ANNOTATE(loss_cand_fa, loss_matrix, pfam_abs, sprot_abs, morgs_abs, 'loss_')
        loss_annotated_matrix = LOSS_ANNOTATE.out.annotated_matrix
        loss_tblastn_summary  = LOSS_VALIDATE.out.summary
    }

    REPORT(
        ANNOTATE.out.annotated_matrix,
        novelty_tblastn_summary,
        SUMMARIZE.out.novelties,
        cand_fa,
        cand_cluster_tsv,
        loss_annotated_matrix,
        loss_tblastn_summary,
        loss_cand_cluster_tsv,
        file(params.config)
    )
}
