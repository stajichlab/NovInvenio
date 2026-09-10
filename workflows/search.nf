nextflow.enable.dsl=2

include { PHMMER_SEARCH         } from '../modules/phmmer'
include { DIAMOND_MAKEDB        } from '../modules/diamond'
include { DIAMOND_SEARCH        } from '../modules/diamond'
include { BLAST_MAKEDB          } from '../modules/blast'
include { BLAST_SEARCH          } from '../modules/blast'
include { PARSE_HITS            } from '../modules/parse_hits'
include { BUILD_PRESENCE_MATRIX } from '../modules/build_presence_matrix'
include { PHMMER_SELF           } from '../modules/self_search'
include { DIAMOND_SELF          } from '../modules/self_search'
include { BLAST_SELF            } from '../modules/self_search'
include { PARSE_SELF_HITS       } from '../modules/parse_self_hits'
include { EXTRACT_PROTEIN_DESCRIPTIONS } from '../modules/protein_descriptions'

workflow SEARCH {
    take:
    ingroup_ch   // [meta, protein_fa]  — ingroup proteomes
    outgroup_ch  // [meta, protein_fa]  — outgroup proteomes
    config_csv   // path to analysis CSV

    main:
    all_proteomes_ch = ingroup_ch.mix(outgroup_ch)

    // Run the selected search tool. diamond/blast build a per-proteome database
    // once (storeDir-cached) and reuse it across every query, rather than
    // rebuilding the target DB for each of the |IN|×(N-1) pairwise jobs.
    // phmmer searches the target FASTA directly, so no DB step is needed.
    // Self-vs-self searches (rank-2 = best within-proteome paralog) feed the
    // paralog-competition filter in build_presence_matrix.py.
    if (params.run_tool == 'phmmer') {
        pairs_ch = ingroup_ch
            .combine(all_proteomes_ch)
            .filter { meta_q, fa_q, meta_t, fa_t -> meta_q.id != meta_t.id }
        raw_hits_ch = PHMMER_SEARCH(pairs_ch)
        raw_self_ch = PHMMER_SELF(ingroup_ch)
    }
    else if (params.run_tool == 'diamond') {
        db_ch = DIAMOND_MAKEDB(all_proteomes_ch)
        pairs_ch = ingroup_ch
            .combine(db_ch)
            .filter { meta_q, fa_q, meta_t, db_t -> meta_q.id != meta_t.id }
        raw_hits_ch = DIAMOND_SEARCH(pairs_ch)
        raw_self_ch = DIAMOND_SELF(ingroup_ch)
    }
    else if (params.run_tool == 'blast') {
        db_ch = BLAST_MAKEDB(all_proteomes_ch)
        pairs_ch = ingroup_ch
            .combine(db_ch)
            .filter { meta_q, fa_q, meta_t, db_t -> meta_q.id != meta_t.id }
        raw_hits_ch = BLAST_SEARCH(pairs_ch)
        raw_self_ch = BLAST_SELF(ingroup_ch)
    }
    else {
        error "Unknown --run_tool '${params.run_tool}': choose phmmer, diamond, or blast"
    }

    PARSE_HITS(raw_hits_ch)
    PARSE_SELF_HITS(raw_self_ch)

    // Collect all parsed TSVs and per-species paralog cutoffs; build presence/absence matrix.
    // Matrix building waits for both pairwise hits and all paralog cutoff files.
    BUILD_PRESENCE_MATRIX(
        PARSE_HITS.out.map { meta_pair, tsv -> tsv }.collect(),
        PARSE_SELF_HITS.out.tsv.map { meta, tsv -> tsv }.collect(),
        config_csv,
        'IN',
        params.ingroup_min_frac,
        0.0,
        'presence_matrix.tsv',
        'candidates.txt',
        'presence_matrix.evalues.tsv',
        'presence_matrix.targets.tsv'
    )

    // Shared protein_id -> gene_name/description lookup (issue: TBLASTN popup query
    // name + pairwise-hit target name) -- one pass over every proteome's own FASTA
    // headers, reused by both BUILD_ALIGNMENT_SHARDS (validate.nf) and REPORT
    // (report.nf) so it's built exactly once per run, not once per consumer.
    EXTRACT_PROTEIN_DESCRIPTIONS(all_proteomes_ch.map { meta, fa -> fa }.collect())

    emit:
    hits          = PARSE_HITS.out             // [meta_pair, parsed_tsv]
    matrix        = BUILD_PRESENCE_MATRIX.out.matrix
    candidates    = BUILD_PRESENCE_MATRIX.out.candidates
    evalues       = BUILD_PRESENCE_MATRIX.out.evalues
    targets       = BUILD_PRESENCE_MATRIX.out.targets
    self_hits     = PARSE_SELF_HITS.out.tsv    // [meta, self_hits_tsv] — one per ingroup proteome
    descriptions  = EXTRACT_PROTEIN_DESCRIPTIONS.out.tsv
}
