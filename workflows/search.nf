nextflow.enable.dsl=2

include { PHMMER_SEARCH         } from '../modules/phmmer'
include { DIAMOND_SEARCH        } from '../modules/diamond'
include { BLAST_SEARCH          } from '../modules/blast'
include { PARSE_HITS            } from '../modules/parse_hits'
include { BUILD_PRESENCE_MATRIX } from '../modules/build_presence_matrix'
include { PHMMER_SELF           } from '../modules/self_search'
include { DIAMOND_SELF          } from '../modules/self_search'
include { BLAST_SELF            } from '../modules/self_search'
include { PARSE_SELF_HITS       } from '../modules/parse_self_hits'

workflow SEARCH {
    take:
    ingroup_ch   // [meta, protein_fa]  — ingroup proteomes
    outgroup_ch  // [meta, protein_fa]  — outgroup proteomes
    config_csv   // path to analysis CSV

    main:
    all_proteomes_ch = ingroup_ch.mix(outgroup_ch)

    // Every ingroup proteome searched against all other proteomes (no self-hits)
    pairs_ch = ingroup_ch
        .combine(all_proteomes_ch)
        .filter { meta_q, fa_q, meta_t, fa_t -> meta_q.id != meta_t.id }

    // Run the selected search tool
    raw_hits_ch = params.run_tool == 'phmmer'  ? PHMMER_SEARCH(pairs_ch)  :
                  params.run_tool == 'diamond' ? DIAMOND_SEARCH(pairs_ch) :
                  params.run_tool == 'blast'   ? BLAST_SEARCH(pairs_ch)   :
                  { error "Unknown --run_tool '${params.run_tool}': choose phmmer, diamond, or blast" }()

    PARSE_HITS(raw_hits_ch)

    // Self-vs-self searches: each ingroup proteome against itself (top 2 hits, relaxed e-value)
    // Rank-1 = self-hit; rank-2 = best within-proteome paralog — used to calibrate score cutoffs
    raw_self_ch = params.run_tool == 'phmmer'  ? PHMMER_SELF(ingroup_ch)  :
                  params.run_tool == 'diamond' ? DIAMOND_SELF(ingroup_ch) :
                  params.run_tool == 'blast'   ? BLAST_SELF(ingroup_ch)   :
                  { error "Unknown --run_tool '${params.run_tool}'" }()

    PARSE_SELF_HITS(raw_self_ch)

    // Collect all parsed TSVs and per-species paralog cutoffs; build presence/absence matrix.
    // Matrix building waits for both pairwise hits and all paralog cutoff files.
    BUILD_PRESENCE_MATRIX(
        PARSE_HITS.out.map { meta_pair, tsv -> tsv }.collect(),
        PARSE_SELF_HITS.out.tsv.map { meta, tsv -> tsv }.collect(),
        config_csv
    )

    emit:
    hits        = PARSE_HITS.out             // [meta_pair, parsed_tsv]
    matrix      = BUILD_PRESENCE_MATRIX.out.matrix
    candidates  = BUILD_PRESENCE_MATRIX.out.candidates
    self_hits   = PARSE_SELF_HITS.out.tsv    // [meta, self_hits_tsv] — one per ingroup proteome
}
