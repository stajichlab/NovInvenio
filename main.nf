#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { SEARCH   } from './workflows/search'
include { CONTEXT_SEARCH } from './workflows/context_search'
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
include { UNIPROT_LINK } from './modules/uniprot_link'
include { UNIPROT_INDEX_BUILD } from './workflows/uniprot_index'
include { SUMMARIZE } from './workflows/summarize'
include { REPORT   } from './workflows/report'
include { NOVELTY_DISCOVERY } from './workflows/novelty_discovery'
include { NOVELTY_SCREEN } from './workflows/novelty_screen'
include { EMPTY_LOSS_STUB } from './modules/empty_loss_stub'
include { EMPTY_EVALUES_STUB } from './modules/empty_evalues_stub'
// A process can only be invoked once per (mutually-exclusive if/else) branch, so the
// context-search stub (issue #48) needs its own aliased imports alongside the e-value
// stub (issue #44) -- both are the same trivial "touch an empty file" process.
include { EMPTY_EVALUES_STUB as EMPTY_CONTEXT_MATRIX_STUB  } from './modules/empty_evalues_stub'
include { EMPTY_EVALUES_STUB as EMPTY_CONTEXT_EVALUES_STUB } from './modules/empty_evalues_stub'

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

// True if at least one config row in GROUP `grp` (after alias normalization) has a
// non-empty DNA value -- used by the launch-time TBLASTN panel check (issue #166).
def panelHasDna(List rows, String grp) {
    return rows.any { row -> normalizeGroup(row.GROUP?.trim()) == grp && row.DNA?.trim() }
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

workflow NOVINVENIO {
    if (!params.config)   error "ERROR: --config <analysis_csv> is required"
    if (!params.data_dir) error "ERROR: --data_dir <fasta_directory> is required"
    if (!file(params.config).exists())   error "ERROR: --config file not found: ${params.config}"
    if (!file(params.data_dir).isDirectory()) error "ERROR: --data_dir is not a directory: ${params.data_dir}"
    if (params.run_tool !in ['phmmer', 'diamond', 'blast']) error "ERROR: --run_tool must be phmmer, diamond, or blast (got: ${params.run_tool})"
    if (params.diamond_sensitivity !in ['', 'sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive']) error "ERROR: --diamond_sensitivity must be empty (default), sensitive, more-sensitive, very-sensitive, or ultra-sensitive (got: ${params.diamond_sensitivity})"
    if (params.other_coverage_floor_qcov && params.run_tool == 'phmmer') error "ERROR: --other_coverage_floor_qcov needs alignment coverage (qcov), which phmmer --tblout does not report. Use --run_tool diamond or blast, or drop the floor (issue #158)."
    if (params.other_coverage_floor_qcov && params.cluster_tool == 'mmseqs') log.warn "--other_coverage_floor_qcov only judges pairwise hits (pairwise matrix, novelty_discovery singletons, context search); --cluster_tool mmseqs family presence uses --hmm_presence_cov instead, so the floor has no effect here."
    if (params.cluster_tool !in ['pairwise', 'mmseqs', 'novelty_discovery']) error "ERROR: --cluster_tool must be pairwise, mmseqs, or novelty_discovery (got: ${params.cluster_tool})"

    // TBLASTN validation needs at least one genome per searched panel. DNA is optional per
    // species, but if a whole panel has none, VALIDATE's / NOVELTY_DISCOVERY's
    // TBLASTN.out...collect() emits nothing, SUMMARIZE_TBLASTN never runs, and neither do
    // SUMMARIZE or REPORT -- the run still exits successfully with no reports (issue #166).
    // Fail at launch instead.
    def cfg_rows = file(params.config).splitCsv(header: true)
    def dna_panels = params.cluster_tool == 'novelty_discovery'
        ? ['DISCOVERY_OUT': 'phase-1 TBLASTN validation']
        : ['OUT': 'novelty-direction TBLASTN validation', 'IN': 'loss-direction TBLASTN validation']
    dna_panels.each { grp, use ->
        if (!panelHasDna(cfg_rows, grp)) error "ERROR: no ${grp} row in ${params.config} has a DNA file, but ${use} needs at least one genome. Add a DNA column value for at least one ${grp} species."
    }

    // Resolve DB paths to absolute at launch time and pass them as val inputs —
    // params mutations do not reliably propagate into process script closures.
    def pfam_abs  = params.pfam_hmm       ? file(params.pfam_hmm).toAbsolutePath().toString()       : ''
    def sprot_abs = params.swissprot_dmnd ? file(params.swissprot_dmnd).toAbsolutePath().toString() : ''
    def morgs_abs = params.modelorgs_config
        ? file(params.modelorgs_config).toAbsolutePath().toString()
        : ''
    // Report-only: resolves each species' optional GFF3 config-CSV column (chrom/start
    // columns in novelties.html/core.html/losses.html) -- not staged/channeled like the
    // Protein/DNA FASTAs, just passed through as an absolute path val for bin/make_*.py
    // to re-resolve GFF3 filenames against directly (see lib/gff3_genes.py).
    def data_dir_abs = file(params.data_dir).toAbsolutePath().toString()

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

    // UniProt linking (docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md):
    // every config proteome, when --uniprot_index is set. UniProtDatGz is a deprecated
    // alias kept for one release: its basename's UP... id restricts the own-species
    // proteome (bin/uniprot_link.py --restrict-proteome); the .dat.gz itself is not read.
    uniprot_link_ch = Channel
        .fromPath(params.config)
        .splitCsv(header: true)
        .map { row ->
            def dat = row.UniProtDatGz?.trim()
            def restrict = dat ? file(dat).name.tokenize('_')[0] : null
            if (dat) log.warn "UniProtDatGz is deprecated (${row.Short}); use --uniprot_index. " +
                              "Treated as --restrict-proteome ${restrict}"
            [ [id: row.Short, species: row.Species, taxid: row.NCBI_TaxID?.trim() ?: null,
               uniprot_restrict: restrict],
              resolve_fa(row.Protein, ['pep', 'proteins']) ]
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

    // UniProt linking -- once per config proteome, regardless of IN/OUT/direction; the
    // whole collected list goes to both ANNOTATE calls below (see workflows/annotate.nf's
    // take: comment). Skipped, with an empty list, when --uniprot_index is unset.
    if (params.uniprot_index) {
        UNIPROT_LINK(uniprot_link_ch, file(params.uniprot_index).toAbsolutePath().toString())
        uniprot_xref_files = UNIPROT_LINK.out.tsv.collect().ifEmpty([])
    } else {
        log.info "--uniprot_index not set: UniProt linking skipped"
        uniprot_xref_files = Channel.value([])
    }

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
        // mmseqs/PROFILE_SEARCH doesn't track hit e-values (or target IDs, or protein
        // descriptions) yet (issue #44 follow-up) -- one empty stub file covers all
        // three (read_evalues()/read_targets()/read_descriptions() all treat an empty
        // file as "no evidence available").
        EMPTY_EVALUES_STUB()
        novelty_evalues      = EMPTY_EVALUES_STUB.out.evalues
        novelty_targets      = EMPTY_EVALUES_STUB.out.evalues
        novelty_descriptions = EMPTY_EVALUES_STUB.out.evalues
        // Query-group low-coverage counts (issue #159) come from the pairwise matrix
        // builder only; family-HMM presence has no per-hit qcov.
        novelty_query_lowcov = EMPTY_EVALUES_STUB.out.evalues
        // NEAR_INGROUP/BROAD_OUTGROUP context search (issue #48) is pairwise-only for now
        // -- mmseqs/PROFILE_SEARCH has no self-vs-self paralog cutoffs to filter against.
        EMPTY_CONTEXT_MATRIX_STUB()
        context_matrix     = EMPTY_CONTEXT_MATRIX_STUB.out.evalues
        EMPTY_CONTEXT_EVALUES_STUB()
        context_evalues    = EMPTY_CONTEXT_EVALUES_STUB.out.evalues

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
            NOVELTY_DISCOVERY.out.singleton_query_fa,  // issue #52
            NOVELTY_DISCOVERY.out.paralog_cutoffs,     // issue #52
            file(params.config)
        )
        novelty_matrix     = NOVELTY_SCREEN.out.matrix
        novelty_candidates = NOVELTY_SCREEN.out.candidates
        // Phase-1-scoped e-value evidence (issue #44) -- see NOVELTY_DISCOVERY's emit: block.
        novelty_evalues    = NOVELTY_DISCOVERY.out.evalues
        // novelty_discovery doesn't track target IDs or protein descriptions yet either --
        // same empty-stub convention as the mmseqs branch above.
        EMPTY_EVALUES_STUB()
        novelty_targets      = EMPTY_EVALUES_STUB.out.evalues
        novelty_descriptions = EMPTY_EVALUES_STUB.out.evalues
        novelty_query_lowcov = EMPTY_EVALUES_STUB.out.evalues   // pairwise-only (issue #159)
        // novelty_discovery already has its own NEAR_INGROUP/BROAD_OUTGROUP screen
        // (NOVELTY_SCREEN) -- the pairwise-only context search (issue #48) doesn't apply.
        EMPTY_CONTEXT_MATRIX_STUB()
        context_matrix     = EMPTY_CONTEXT_MATRIX_STUB.out.evalues
        EMPTY_CONTEXT_EVALUES_STUB()
        context_evalues    = EMPTY_CONTEXT_EVALUES_STUB.out.evalues

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
        novelty_matrix       = SEARCH.out.matrix
        novelty_candidates   = SEARCH.out.candidates
        novelty_evalues      = SEARCH.out.evalues
        novelty_targets      = SEARCH.out.targets
        novelty_descriptions = SEARCH.out.descriptions
        novelty_query_lowcov = SEARCH.out.query_lowcov

        CLUSTER(novelty_candidates, ingroup_prot_ch, file(params.config), 'candidates.fa', 'clusters')
        cand_fa          = CLUSTER.out.candidates_fa
        cand_reps        = CLUSTER.out.representatives
        cand_cluster_tsv = CLUSTER.out.cluster_tsv

        // NEAR_INGROUP/BROAD_OUTGROUP context (issue #48): report-only presence/e-value
        // evidence for the already-fixed candidate list, never affecting novelty calling.
        // Runs automatically -- when the config has no NEAR_INGROUP/BROAD_OUTGROUP rows,
        // near_in_prot_ch/broad_out_prot_ch are simply empty and CONTEXT_SEARCH is a no-op
        // (empty context columns in the report, not skipped entirely).
        CONTEXT_SEARCH(
            novelty_candidates,
            SEARCH.out.self_hits,
            ingroup_prot_ch,
            near_in_prot_ch,
            broad_out_prot_ch,
            file(params.config)
        )
        context_matrix  = CONTEXT_SEARCH.out.matrix
        context_evalues = CONTEXT_SEARCH.out.evalues
    }

    // novelty_discovery already produced its own TBLASTN summary (see above); the other two
    // cluster_tool paths still need the generic VALIDATE (TBLASTN vs the OUT proteomes' DNA).
    if (params.cluster_tool != 'novelty_discovery') {
        VALIDATE(cand_reps, outgroup_dna_ch, cand_cluster_tsv, 'tblastn_summary.tsv',
                 novelty_candidates, 'alignments', novelty_descriptions)
        novelty_tblastn_summary = VALIDATE.out.summary
    }

    ANNOTATE(cand_fa, novelty_matrix, pfam_abs, sprot_abs, morgs_abs, '', uniprot_xref_files)

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
        loss_matrix     = PROFILE_LOSS_SEARCH.out.matrix
        loss_candidates = PROFILE_LOSS_SEARCH.out.candidates

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
        // to assemble docs/<project>/, so stub the three loss artifacts instead.
        EMPTY_LOSS_STUB(ANNOTATE.out.annotated_matrix)
        loss_annotated_matrix   = EMPTY_LOSS_STUB.out.matrix
        loss_tblastn_summary    = EMPTY_LOSS_STUB.out.tblastn_summary
        loss_cand_cluster_tsv   = EMPTY_LOSS_STUB.out.cluster_tsv
    }
    else {
        // See workflows/loss_search.nf for why this needs its own search direction.
        LOSS_SEARCH(ingroup_prot_ch, outgroup_prot_ch, file(params.config))
        loss_matrix     = LOSS_SEARCH.out.matrix
        loss_candidates = LOSS_SEARCH.out.candidates

        LOSS_CLUSTER(LOSS_SEARCH.out.candidates, outgroup_prot_ch, file(params.config), 'loss_candidates.fa', 'loss_clusters')
        loss_cand_fa          = LOSS_CLUSTER.out.candidates_fa
        loss_cand_reps        = LOSS_CLUSTER.out.representatives
        loss_cand_cluster_tsv = LOSS_CLUSTER.out.cluster_tsv
    }

    if (params.cluster_tool != 'novelty_discovery') {
        LOSS_VALIDATE(loss_cand_reps, ingroup_dna_ch, loss_cand_cluster_tsv, 'loss_tblastn_summary.tsv',
                      loss_candidates, 'loss_alignments', novelty_descriptions)
        LOSS_ANNOTATE(loss_cand_fa, loss_matrix, pfam_abs, sprot_abs, morgs_abs, 'loss_', uniprot_xref_files)
        loss_annotated_matrix = LOSS_ANNOTATE.out.annotated_matrix
        loss_tblastn_summary  = LOSS_VALIDATE.out.summary
    }

    REPORT(
        ANNOTATE.out.annotated_matrix,
        novelty_tblastn_summary,
        SUMMARIZE.out.novelties,
        cand_fa,
        cand_cluster_tsv,
        novelty_evalues,
        novelty_targets,
        novelty_descriptions,
        context_matrix,
        context_evalues,
        novelty_query_lowcov,
        loss_annotated_matrix,
        loss_tblastn_summary,
        loss_cand_cluster_tsv,
        file(params.config),
        data_dir_abs
    )
}

// One-time UniProt library index build (docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md):
//   nextflow run main.nf --build_uniprot_index --uniprot_library <dir> \
//       --uniprot_library_csv <csv name> --uniprot_index <out dir>
// (Nextflow's strict parser does not support -entry, so a param selects it.)
workflow UNIPROT_INDEX {
    main:
    if (!params.uniprot_library || !params.uniprot_library_csv || !params.uniprot_index) {
        error "--build_uniprot_index needs --uniprot_library, --uniprot_library_csv and --uniprot_index"
    }
    if (file("${params.uniprot_index}/manifest.json").exists()) {
        // The index is complete (manifest.json is written last): skip re-parsing the
        // whole library, which storeDir on the final step alone would not prevent.
        log.info "UniProt index already built: ${params.uniprot_index}/manifest.json exists; nothing to do"
    } else {
        UNIPROT_INDEX_BUILD(Channel.value(file(params.uniprot_library).toAbsolutePath().toString()))
    }
}

workflow {
    if (Helpers.asBool(params.build_uniprot_index)) {
        UNIPROT_INDEX()
    } else {
        NOVINVENIO()
    }
}
