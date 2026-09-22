nextflow.enable.dsl=2

// PANGENOME_PROFILE — reusable pangenome cluster-profiling subworkflow.
//
// Ported from a one-off study script chain (NovInvenio_Investigations'
// Afumigatus_pangenome study, `bin/*.py` + `run_*.sh` SLURM launchers) into a
// composable DSL2 subworkflow any future study can invoke against its own
// species set. See NEXTFLOW_MIGRATION_NOTES.md (in that study's directory)
// for the full design record this implements.
//
// Pipeline (mirrors that study's PANGENOME_CLUSTER_PROFILE_NOTES.md "Run
// order" section):
//   PREFIX_PROTEOME/PREFIX_GENOME (per-strain) -> CONCAT_PROTEOMES -> CLUSTER_TIER1
//     -> PRESENCE_MATRIX -> GENE_POSITIONS
//     -> [rescue branch] EXTRACT_ABSENT_QUERIES -> per-strain scatter
//        (MAKE_STRAIN_GENOME_DB + TBLASTN_PER_STRAIN) -> RESCUE_PASS (issue
//        #133 structural filter: gene_positions + cluster_tsv + rep_fasta)
//        -> EXTRACT_RESCUE_POSITIONS
//     -> [clade branch] MASH_SKETCH (x2: all-strains for dedup, ingroup-only
//        for clades) -> DEREPLICATE / ASSIGN_CLADES -> FILL_TAXON_GROUP
//     -> FREQUENCY_BINS -> COOCCURRENCE
//     -> [positions branch] FAMILY_POSITIONS (gene_positions + rescue_positions)
//     -> [captain branch, optional] HMMFETCH_CAPTAIN? -> CAPTAIN_HMMSEARCH
//     -> PAIR_CLASSIFICATION
//
// GENE_POSITIONS runs right after PRESENCE_MATRIX (moved ahead of the rescue
// branch when issue #133 added the structural rescue filter) rather than
// alongside FAMILY_POSITIONS where it originally lived -- it only ever
// depended on the samplesheet/GFF3s, never on clustering or rescue, so
// nothing about what it computes changed, only when.
//
// Two documented departures from a literal 1:1 port of the originating
// study, both explained further at their point of use below:
//   - MASH_SKETCH runs TWICE (once over the full strain set for DEREPLICATE,
//     once ingroup-only for ASSIGN_CLADES), not once for both consumers --
//     see modules/pangenome/mash.nf's module docstring for why a single
//     combined sketch would reintroduce a real outgroup-domination artifact.
//   - The samplesheet's TaxonGroup fill (manual/study-specific in the
//     originating study) is automated here as a generic "fill only what's
//     blank from the Mash clade labels" policy (FILL_TAXON_GROUP) -- a study
//     wanting a richer priority scheme (e.g. cross-referencing a published
//     population table) should pre-fill TaxonGroup in its own samplesheet
//     before calling this subworkflow; this subworkflow never overwrites an
//     existing value.

include { PREFIX_PROTEOME; PREFIX_GENOME; CONCAT_PROTEOMES; CLUSTER_TIER1 } from '../modules/pangenome/prefix_and_cluster'
include { PRESENCE_MATRIX }                                                 from '../modules/pangenome/presence_matrix'
include { EXTRACT_ABSENT_QUERIES; MAKE_STRAIN_GENOME_DB; TBLASTN_PER_STRAIN; RESCUE_PASS } from '../modules/pangenome/rescue'
include { MASH_SKETCH as MASH_SKETCH_ALL }      from '../modules/pangenome/mash'
include { MASH_SKETCH as MASH_SKETCH_INGROUP }  from '../modules/pangenome/mash'
include { DEREPLICATE; ASSIGN_CLADES; FILL_TAXON_GROUP }                    from '../modules/pangenome/mash'
include { FREQUENCY_BINS; COOCCURRENCE }                                    from '../modules/pangenome/frequency_cooccurrence'
include { GENE_POSITIONS; EXTRACT_RESCUE_POSITIONS; FAMILY_POSITIONS }      from '../modules/pangenome/positions'
include { HMMFETCH_CAPTAIN; CAPTAIN_HMMSEARCH }                             from '../modules/pangenome/captain'
include { PAIR_CLASSIFICATION }                                             from '../modules/pangenome/pair_classification'
include { BUILD_ISLANDS; MARKER_HMMSEARCH }                                from '../modules/pangenome/islands'
include { SELECT_BACKGROUND_REPS; HMMPRESS_PFAM; FAMILY_PFAM_SCAN;
          MERGE_PFAM_DOMTBLOUT; DOMAIN_ENRICHMENT }                 from '../modules/pangenome/pfam_enrichment'
include { REPORT_TABLES; REPORT_RENDER; DIAGNOSTICS }                      from '../modules/pangenome/report'
include { ASSEMBLY_QUALITY_QC }                                            from '../modules/pangenome/assembly_quality_qc'
include { ISLAND_SYNTENY }                                                 from '../modules/pangenome/island_synteny'
include { LEIDEN_MODULES; MODULE_DOMAINS }                                 from '../modules/pangenome/trans_modules'
include { PFAM2GO } from '../modules/pangenome/pfam2go'
include { EMPTY_EVALUES_STUB as EMPTY_RESCUE_POSITIONS_STUB } from '../modules/empty_evalues_stub'
include { EMPTY_EVALUES_STUB as EMPTY_CAPTAIN_STUB }          from '../modules/empty_evalues_stub'
include { EMPTY_EVALUES_STUB as EMPTY_INVENTORY_STUB }        from '../modules/empty_evalues_stub'
include { EMPTY_EVALUES_STUB as EMPTY_RESCUE_TSV_STUB }       from '../modules/empty_evalues_stub'
include { EMPTY_EVALUES_STUB as EMPTY_RESCUE_FUNNEL_STUB }    from '../modules/empty_evalues_stub'

workflow PANGENOME_PROFILE {
    take:
    samplesheet      // path: GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup CSV
    data_dir_abs     // val:  absolute path to the directory holding pep/, dna/, gff3/
    gff3_dir_abs     // val:  absolute path to the GFF3 directory (usually "${data_dir_abs}/gff3")
    samples_ch       // channel: [meta, protein_fa, dna_fa] — one entry per strain in
                      //   {ingroup_label, outgroup_label} (the clustering-input strain set)
    ingroup_dna_ch   // channel: [meta, dna_fa] — ingroup strains only (for clade assignment)

    main:
    // --- 1. Prefix + concat + tier-1 cluster ------------------------------
    PREFIX_PROTEOME(samples_ch.map { meta, prot, dna -> [meta, prot] })
    PREFIX_GENOME(samples_ch.map { meta, prot, dna -> [meta, dna] })

    CONCAT_PROTEOMES(PREFIX_PROTEOME.out.fasta.map { meta, fa -> fa }.collect())
    CLUSTER_TIER1(CONCAT_PROTEOMES.out.fasta)

    // --- 2. Presence matrix ------------------------------------------------
    PRESENCE_MATRIX(CLUSTER_TIER1.out.cluster_tsv, samplesheet)

    // --- 2b. Gene positions --------------------------------------------------
    // protein_dir_abs is always "<data_dir_abs>/pep" -- the fixed layout
    // build_study_config.py always produces, same fixed-subdir convention
    // gff3_dir_abs's own default already uses in pangenome.nf. Needed so
    // GENE_POSITIONS can cross-check resolved GFF3 IDs against each strain's
    // real protein FASTA headers, rather than trust GFF3 attribute presence
    // alone. See this file's header comment for why this now runs here
    // rather than alongside FAMILY_POSITIONS.
    protein_dir_abs = "${data_dir_abs}/pep"
    GENE_POSITIONS(samplesheet, gff3_dir_abs, protein_dir_abs)

    // --- 3. Rescue pass (optional; per-strain scatter, see modules/pangenome/rescue.nf) ---
    if (params.pangenome_rescue_enable) {
        EXTRACT_ABSENT_QUERIES(PRESENCE_MATRIX.out.matrix, CLUSTER_TIER1.out.rep_fasta)

        // Join each strain's absent-family query FASTA back to that same
        // strain's prefixed genome FASTA by Short ID (parsed from the query
        // FASTA's own filename, "<Short>.absent.fa" -- see
        // pangenome_extract_absent_family_queries.py). Strains with zero
        // absent families never got a query file and are naturally excluded
        // here (no tblastn work needed for them).
        query_ch = EXTRACT_ABSENT_QUERIES.out.query_fastas
            .ifEmpty([])
            .flatten()
            .map { fa -> tuple(fa.name.replaceFirst(/\.absent\.fa$/, ''), fa) }
        genome_ch = PREFIX_GENOME.out.fasta.map { meta, fa -> tuple(meta.id, meta, fa) }
        joined = query_ch.join(genome_ch)
            .map { strain_id, query_fa, meta, genome_fa -> tuple(meta, genome_fa, query_fa) }

        MAKE_STRAIN_GENOME_DB(joined.map { meta, genome_fa, query_fa -> tuple(meta, genome_fa) })

        db_ch    = MAKE_STRAIN_GENOME_DB.out.map { meta, db -> tuple(meta.id, meta, db) }
        query_only_ch = joined.map { meta, genome_fa, query_fa -> tuple(meta.id, query_fa) }
        tblastn_in = db_ch.join(query_only_ch)
            .map { id, meta, db, query_fa -> tuple(meta, db, query_fa) }

        TBLASTN_PER_STRAIN(tblastn_in)
        tblastn_tsv_files = TBLASTN_PER_STRAIN.out.tsv.map { meta, tsv -> tsv }.collect().ifEmpty([])

        RESCUE_PASS(
            PRESENCE_MATRIX.out.matrix, tblastn_tsv_files,
            GENE_POSITIONS.out.positions, CLUSTER_TIER1.out.cluster_tsv, CLUSTER_TIER1.out.rep_fasta,
        )
        rescued_matrix = RESCUE_PASS.out.matrix
        rescue_funnel = RESCUE_PASS.out.funnel

        EXTRACT_RESCUE_POSITIONS(rescued_matrix, tblastn_tsv_files)
        rescue_positions = EXTRACT_RESCUE_POSITIONS.out.positions
    }
    else {
        rescued_matrix = PRESENCE_MATRIX.out.matrix
        EMPTY_RESCUE_POSITIONS_STUB()
        rescue_positions = EMPTY_RESCUE_POSITIONS_STUB.out.evalues
        EMPTY_RESCUE_TSV_STUB()
        tblastn_tsv_files = EMPTY_RESCUE_TSV_STUB.out.evalues
        // issue #134: no rescue funnel was ever computed -- DIAGNOSTICS
        // reports rescue_redundancy as not_computed rather than erroring.
        EMPTY_RESCUE_FUNNEL_STUB()
        rescue_funnel = EMPTY_RESCUE_FUNNEL_STUB.out.evalues
    }

    // --- 4. Strain dedup + clade assignment --------------------------------
    // MASH_SKETCH runs twice (see this file's header comment / mash.nf's
    // module docstring): the full strain set for dedup, ingroup-only for
    // clade assignment -- a single shared sketch would let an outgroup
    // species dominate the ingroup clade signal.
    // MASH_SKETCH_ALL's only consumer is DEREPLICATE -- skip it entirely
    // when pangenome_dereplicate is disabled instead of sketching for
    // nothing (see review item 7).
    if (params.pangenome_dereplicate) {
        all_dna_ch = samples_ch.map { meta, prot, dna -> dna }.collect()
        MASH_SKETCH_ALL(all_dna_ch, 'all_strains')
        DEREPLICATE(samplesheet, data_dir_abs, MASH_SKETCH_ALL.out.dist_tsv)
        strain_inventory = DEREPLICATE.out.inventory
    }
    else {
        EMPTY_INVENTORY_STUB()
        strain_inventory = EMPTY_INVENTORY_STUB.out.evalues
    }

    if (params.pangenome_assign_clades) {
        ingroup_dna_files = ingroup_dna_ch.map { meta, dna -> dna }.collect()
        MASH_SKETCH_INGROUP(ingroup_dna_files, 'ingroup')
        ASSIGN_CLADES(samplesheet, data_dir_abs, MASH_SKETCH_INGROUP.out.dist_tsv)
        FILL_TAXON_GROUP(samplesheet, ASSIGN_CLADES.out.clades)
        effective_samplesheet = FILL_TAXON_GROUP.out.samplesheet
    }
    else {
        effective_samplesheet = samplesheet
    }

    // --- 5. Frequency binning + co-occurrence ------------------------------
    FREQUENCY_BINS(rescued_matrix, effective_samplesheet, strain_inventory)
    COOCCURRENCE(rescued_matrix, FREQUENCY_BINS.out.table, effective_samplesheet, strain_inventory)

    // --- 5b. Assembly-quality vs pangenome-content QC (issue #130) --------
    // Runs unconditionally (cheap; no gating param) -- see
    // modules/pangenome/assembly_quality_qc.nf's module docstring. Diagnostic
    // only: never excludes a strain or alters a presence call.
    ASSEMBLY_QUALITY_QC(
        effective_samplesheet, data_dir_abs, rescued_matrix,
        FREQUENCY_BINS.out.table, GENE_POSITIONS.out.positions, CLUSTER_TIER1.out.cluster_tsv,
    )

    // --- 6. Family positions -------------------------------------------------
    // GENE_POSITIONS itself now runs at step 2b, above -- see this file's
    // header comment.
    FAMILY_POSITIONS(GENE_POSITIONS.out.positions, CLUSTER_TIER1.out.cluster_tsv, rescue_positions)

    // --- 7. Optional captain/mobile-element marker gene evidence -----------
    def captain_requested = params.pangenome_captain_hmm || (params.pangenome_captain_hmm_name && params.pangenome_pfam_hmm)
    if (captain_requested) {
        if (params.pangenome_captain_hmm) {
            captain_hmm_ch = Channel.value(file(params.pangenome_captain_hmm))
        }
        else {
            HMMFETCH_CAPTAIN(Channel.value(file(params.pangenome_pfam_hmm)), params.pangenome_captain_hmm_name)
            captain_hmm_ch = HMMFETCH_CAPTAIN.out.hmm
        }
        CAPTAIN_HMMSEARCH(captain_hmm_ch, CONCAT_PROTEOMES.out.fasta)
        captain_tblout = CAPTAIN_HMMSEARCH.out.tblout
    }
    else {
        EMPTY_CAPTAIN_STUB()
        captain_tblout = EMPTY_CAPTAIN_STUB.out.evalues
    }

    // --- 8. Pair classification (physical linkage / mobile-element mechanism) ---
    PAIR_CLASSIFICATION(
        COOCCURRENCE.out.pairs,
        FAMILY_POSITIONS.out.positions,
        CLUSTER_TIER1.out.cluster_tsv,
        captain_tblout,
    )

    // --- 8b. Leiden trans-module detection (chained co-occurring families) ---
    // Runs unconditionally (cheap; degrades gracefully to empty output for a
    // study with zero `trans` pairs, e.g. a single-species ingroup -- see
    // bin/pangenome_detect_trans_modules.py's module docstring).
    LEIDEN_MODULES(PAIR_CLASSIFICATION.out.classification)

    // --- 9. Accessory islands + Pfam functional enrichment (optional) --------
    // Named marker searches (0+): ONE MARKER_HMMSEARCH invocation over a
    // Channel.fromList of [name, hmm_path] tuples runs one task per marker via
    // Nextflow's normal channel-based fan-out -- not a loop calling the process
    // repeatedly (DSL2 forbids invoking a process more than once per workflow).
    // CONCAT_PROTEOMES.out.fasta is a single-emission channel; `.first()` turns
    // it into a value reused for every marker task rather than being consumed
    // by only the first one.
    if (params.pangenome_island_pfam_hmm) {
        def marker_names_list = params.pangenome_marker_names ? params.pangenome_marker_names.split(',')*.trim() as List : []
        def marker_hmm_paths_list = params.pangenome_marker_hmm_paths ? params.pangenome_marker_hmm_paths.split(',') as List : []
        if (marker_names_list.size() != marker_hmm_paths_list.size()) {
            error "ERROR: --pangenome_marker_names and --pangenome_marker_hmm_paths must have " +
                  "the same number of comma-separated entries (got ${marker_names_list.size()} names, " +
                  "${marker_hmm_paths_list.size()} paths)"
        }
        // Fail fast at parse time (not deep inside a later process, e.g. a
        // broken output filename or a malformed --marker_tblout CLI arg): a
        // natural-but-wrong list like 'captain, sm_backbone' (space after
        // the comma) would otherwise produce a marker literally named
        // " sm_backbone". Names are trimmed above; validate what remains is
        // a well-formed identifier.
        def invalid_marker_names = marker_names_list.findAll { !(it ==~ /^\w+$/) }
        if (invalid_marker_names) {
            error "ERROR: --pangenome_marker_names entries must be non-empty and contain only " +
                  "word characters (letters, digits, underscore) -- invalid: ${invalid_marker_names}"
        }

        if (marker_names_list) {
            marker_input_ch = Channel.fromList(
                [marker_names_list, marker_hmm_paths_list.collect { file(it) }].transpose()
            )
            MARKER_HMMSEARCH(marker_input_ch, CONCAT_PROTEOMES.out.fasta.first())
            // .toList() (NOT .collect(), which flattens [[n,p],[n,p]] to
            // [n,p,n,p] by default) keeps each [name, tblout] pair intact as one
            // list-of-pairs value; .collect{it[0]}/.collect{it[1]} then split
            // that into the two PARALLEL lists BUILD_ISLANDS's val+path inputs
            // expect (see Task 4).
            marker_pairs_ch = MARKER_HMMSEARCH.out.result.toList()
            marker_names_ch = marker_pairs_ch.map { pairs -> pairs.collect { it[0] } }
            marker_tblout_files_ch = marker_pairs_ch.map { pairs -> pairs.collect { it[1] } }
        } else {
            marker_names_ch = Channel.value([])
            marker_tblout_files_ch = Channel.value([])
        }

        BUILD_ISLANDS(
            FAMILY_POSITIONS.out.positions,
            FREQUENCY_BINS.out.table,
            PAIR_CLASSIFICATION.out.classification,
            CLUSTER_TIER1.out.cluster_tsv,
            marker_names_ch,
            marker_tblout_files_ch,
        )

        SELECT_BACKGROUND_REPS(CLUSTER_TIER1.out.rep_fasta, FREQUENCY_BINS.out.table)

        // Scatter the Pfam scan across chunks of the background set (issue
        // #112): Nextflow submits one independent SLURM job per chunk, so a
        // 2-4 h monolithic hmmscan that twice hit a 2 h wall-clock cap
        // becomes N short jobs, each individually retryable. The database is
        // hmmpress'd once and shared. MERGE_PFAM_DOMTBLOUT reassembles the
        // single pfam.domtblout every downstream consumer already expects,
        // so nothing below this point changes.
        // `as int` is required, not cosmetic: a value supplied on the command
        // line (--pangenome_pfam_chunk_size 500) arrives as a String, and
        // splitFasta's `by:` rejects it ("Value don't match: class
        // java.lang.Integer"). The nextflow.config default is already an
        // Integer, so without this the failure appears ONLY when a user
        // overrides the default.
        HMMPRESS_PFAM(file(params.pangenome_island_pfam_hmm))
        pfam_chunks_ch = SELECT_BACKGROUND_REPS.out.fasta
            .splitFasta(by: params.pangenome_pfam_chunk_size as int, file: true)
        FAMILY_PFAM_SCAN(pfam_chunks_ch, HMMPRESS_PFAM.out.db.collect())
        MERGE_PFAM_DOMTBLOUT(FAMILY_PFAM_SCAN.out.domtblout.collect())
        pfam_domtblout = MERGE_PFAM_DOMTBLOUT.out.domtblout

        DOMAIN_ENRICHMENT(BUILD_ISLANDS.out.islands, pfam_domtblout, FREQUENCY_BINS.out.table)
        MODULE_DOMAINS(LEIDEN_MODULES.out.family_modules, pfam_domtblout)

        if (params.pangenome_pfam2go) {
            PFAM2GO(DOMAIN_ENRICHMENT.out.enrichment, file(params.pangenome_pfam2go))
            enrichment_for_report = PFAM2GO.out.annotated
        } else {
            enrichment_for_report = DOMAIN_ENRICHMENT.out.enrichment
        }

        REPORT_TABLES(
            BUILD_ISLANDS.out.islands,
            enrichment_for_report,
            PAIR_CLASSIFICATION.out.classification,
            rescued_matrix,
            FREQUENCY_BINS.out.table,
            pfam_domtblout,
            CLUSTER_TIER1.out.cluster_tsv,
            GENE_POSITIONS.out.positions,
        )

        // Issue #134: pipeline diagnostics (rescue_redundancy today; other
        // issue #134 table rows declared not_computed until their own
        // statistic exists as a real pipeline output -- see
        // bin/pangenome_diagnostics.py's module docstring). Feeds the
        // Markdown/HTML banners REPORT_RENDER/ISLAND_SYNTENY prepend.
        DIAGNOSTICS(rescue_funnel)

        REPORT_RENDER(
            FREQUENCY_BINS.out.table,
            rescued_matrix,
            REPORT_TABLES.out.islands_with_domains,
            REPORT_TABLES.out.size_distribution,
            REPORT_TABLES.out.classification_counts,
            enrichment_for_report,
            REPORT_TABLES.out.marker_summary,
            REPORT_TABLES.out.per_strain_summary,
            DIAGNOSTICS.out.banner_md,
        )

        ISLAND_SYNTENY(
            REPORT_TABLES.out.islands_with_domains,
            rescued_matrix,
            FAMILY_POSITIONS.out.positions,
            FAMILY_PFAM_SCAN.out.domtblout,
            samplesheet,
            DIAGNOSTICS.out.banner_html,
        )
    }

    emit:
    cluster_tsv         = CLUSTER_TIER1.out.cluster_tsv
    rep_fasta           = CLUSTER_TIER1.out.rep_fasta
    presence_matrix     = PRESENCE_MATRIX.out.matrix
    rescued_matrix       = rescued_matrix
    strain_inventory     = strain_inventory
    effective_samplesheet = effective_samplesheet
    frequency_table      = FREQUENCY_BINS.out.table
    cooccurring_pairs    = COOCCURRENCE.out.pairs
    family_positions     = FAMILY_POSITIONS.out.positions
    pair_classification  = PAIR_CLASSIFICATION.out.classification
    family_modules       = LEIDEN_MODULES.out.family_modules
    module_summary       = LEIDEN_MODULES.out.module_summary
    assembly_quality_table         = ASSEMBLY_QUALITY_QC.out.table
    assembly_quality_correlations  = ASSEMBLY_QUALITY_QC.out.correlations
    assembly_quality_report        = ASSEMBLY_QUALITY_QC.out.report
}
