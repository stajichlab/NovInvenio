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

// Mirrors workflows/search.nf with the ingroup/outgroup roles swapped:
// outgroup proteomes are the query, looking for genes conserved in the
// outgroup but absent from the ingroup — candidate lineage-specific losses.
//
// This needs its own search direction because SEARCH only ever runs with
// ingroup proteomes as the query (see search.nf) — presence_matrix.tsv has
// no rows sourced from an outgroup protein, so "is this outgroup gene
// present in the ingroup?" cannot be answered from data SEARCH already
// produced. DIAMOND_MAKEDB/BLAST_MAKEDB are storeDir-cached per proteome, so
// re-running them here for the same proteome universe is a cache hit, not
// duplicated work.
workflow LOSS_SEARCH {
    take:
    ingroup_ch    // [meta, protein_fa]  — ingroup proteomes (search target)
    outgroup_ch   // [meta, protein_fa]  — outgroup proteomes (search query)
    config_csv    // path to analysis CSV

    main:
    all_proteomes_ch = ingroup_ch.mix(outgroup_ch)

    if (params.run_tool == 'phmmer') {
        pairs_ch = outgroup_ch
            .combine(all_proteomes_ch)
            .filter { meta_q, fa_q, meta_t, fa_t -> meta_q.id != meta_t.id }
        raw_hits_ch = PHMMER_SEARCH(pairs_ch)
        raw_self_ch = PHMMER_SELF(outgroup_ch)
    }
    else if (params.run_tool == 'diamond') {
        db_ch = DIAMOND_MAKEDB(all_proteomes_ch)
        pairs_ch = outgroup_ch
            .combine(db_ch)
            .filter { meta_q, fa_q, meta_t, db_t -> meta_q.id != meta_t.id }
        raw_hits_ch = DIAMOND_SEARCH(pairs_ch)
        raw_self_ch = DIAMOND_SELF(outgroup_ch)
    }
    else if (params.run_tool == 'blast') {
        db_ch = BLAST_MAKEDB(all_proteomes_ch)
        pairs_ch = outgroup_ch
            .combine(db_ch)
            .filter { meta_q, fa_q, meta_t, db_t -> meta_q.id != meta_t.id }
        raw_hits_ch = BLAST_SEARCH(pairs_ch)
        raw_self_ch = BLAST_SELF(outgroup_ch)
    }
    else {
        error "Unknown --run_tool '${params.run_tool}': choose phmmer, diamond, or blast"
    }

    PARSE_HITS(raw_hits_ch)
    PARSE_SELF_HITS(raw_self_ch)

    // Same paralog-aware filters as SEARCH, run with the query/target groups
    // swapped: a loss candidate is present in >= outgroup_min_frac of the
    // outgroup and absent from every ingroup proteome.
    BUILD_PRESENCE_MATRIX(
        PARSE_HITS.out.map { meta_pair, tsv -> tsv }.collect(),
        PARSE_SELF_HITS.out.tsv.map { meta, tsv -> tsv }.collect(),
        config_csv,
        'OUT',
        params.outgroup_min_frac,
        params.loss_ingroup_max_frac,
        'loss_presence_matrix.tsv',
        'loss_candidates.txt'
    )

    emit:
    hits        = PARSE_HITS.out             // [meta_pair, parsed_tsv]
    matrix      = BUILD_PRESENCE_MATRIX.out.matrix
    candidates  = BUILD_PRESENCE_MATRIX.out.candidates
    self_hits   = PARSE_SELF_HITS.out.tsv    // [meta, self_hits_tsv] — one per outgroup proteome
}
