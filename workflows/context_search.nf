nextflow.enable.dsl=2

// CONTEXT_SEARCH (issue #48) — NEAR_INGROUP/BROAD_OUTGROUP as report-only context for
// the --cluster_tool pairwise pathway.
//
// Novelty filtering is unchanged: a candidate is fixed by the strict IN/OUT search
// (bin/build_presence_matrix.py) before this workflow ever runs. This answers a
// separate, exploratory question for the report -- "does this novelty candidate also
// turn up in a close relative, or a distant lineage?" -- without paying the cost of a
// full ingroup-proteome-vs-NEAR_INGROUP/BROAD_OUTGROUP-proteome pairwise search (which
// would be the same O(N_ingroup_proteins x N_context_proteomes) blowup SEARCH already
// pays for IN vs OUT, just multiplied by however many NEAR_INGROUP/BROAD_OUTGROUP
// species the config lists -- 20-100 per todo/novelty-discovery-screen.md's stated
// range). Searching only the (typically two-orders-of-magnitude smaller) candidate set
// keeps this cheap even as NEAR_INGROUP/BROAD_OUTGROUP scale up.
//
// bin/build_context_query.py extends candidates.txt with each candidate's own
// within-genome paralog (from the ingroup self-vs-self search) so the same
// paralog-competition filter build_presence_matrix.py uses can be applied here too
// (bin/context_presence.py) -- otherwise a candidate's hit against a context proteome
// can't be told apart from a hit better explained by its conserved paralog (the
// HEX-1/eIF5A pattern from the NCU08332 investigation).

include { EXTRACT_CANDIDATES } from './cluster'
include { PARSE_HITS         } from '../modules/parse_hits'
include { DIAMOND_MAKEDB     } from '../modules/diamond'
include { BLAST_MAKEDB       } from '../modules/blast'

process CONTEXT_PHMMER_SEARCH {
    label 'med_cpu'
    tag "context_vs_${meta_t.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(context_query_fa)
    tuple val(meta_t), path(target_fa)

    output:
    tuple val(meta_pair), path("context_vs_${meta_t.id}.phmmer.tblout")

    script:
    meta_pair = [query_id: 'context', target_id: meta_t.id, tool: 'phmmer']
    """
    if [ -s ${context_query_fa} ]; then
        phmmer \
            --cpu ${task.cpus} \
            --tblout context_vs_${meta_t.id}.phmmer.tblout \
            --noali \
            ${context_query_fa} ${target_fa} \
            > /dev/null
    else
        touch context_vs_${meta_t.id}.phmmer.tblout
    fi
    """
}

process CONTEXT_DIAMOND_SEARCH {
    label 'high_cpu'
    tag "context_vs_${meta_t.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(context_query_fa)
    tuple val(meta_t), path(target_db)

    output:
    tuple val(meta_pair), path("context_vs_${meta_t.id}.diamond.tsv")

    script:
    meta_pair = [query_id: 'context', target_id: meta_t.id, tool: 'diamond']
    """
    if [ -s ${context_query_fa} ]; then
        diamond blastp \
            --query ${context_query_fa} \
            --db ${target_db.baseName} \
            --outfmt 6 qseqid sseqid evalue bitscore \
            --evalue ${params.parse_evalue} \
            --threads ${task.cpus} \
            --quiet \
            --out context_vs_${meta_t.id}.diamond.tsv
    else
        touch context_vs_${meta_t.id}.diamond.tsv
    fi
    """
}

process CONTEXT_BLAST_SEARCH {
    label 'high_cpu'
    tag "context_vs_${meta_t.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(context_query_fa)
    tuple val(meta_t), path(target_db)

    output:
    tuple val(meta_pair), path("context_vs_${meta_t.id}.blast.tsv")

    script:
    meta_pair = [query_id: 'context', target_id: meta_t.id, tool: 'blast']
    """
    if [ -s ${context_query_fa} ]; then
        blastp \
            -query ${context_query_fa} \
            -db ${meta_t.id}.blast_db \
            -outfmt "6 qseqid sseqid evalue bitscore" \
            -evalue ${params.parse_evalue} \
            -num_threads ${task.cpus} \
            -out context_vs_${meta_t.id}.blast.tsv
    else
        touch context_vs_${meta_t.id}.blast.tsv
    fi
    """
}

process BUILD_CONTEXT_QUERY {
    label 'low_cpu'
    tag "build_context_query"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(candidates_txt)
    path(paralog_cutoffs)

    output:
    path("context_query.txt"), emit: query_list

    script:
    """
    build_context_query.py \
        --candidates ${candidates_txt} \
        --paralog-cutoffs ${paralog_cutoffs} \
        --output context_query.txt
    """
}

process CONTEXT_PRESENCE {
    label 'low_cpu'
    tag "context_presence"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(hit_tsvs)
    path(candidates_txt)
    path(paralog_cutoffs)
    path(config_csv)
    val(paralog_competition_scope)
    val(default_evalue)

    output:
    path("context_presence.tsv"),         emit: matrix
    path("context_presence.evalues.tsv"), emit: evalues

    script:
    """
    context_presence.py \
        --hits ${hit_tsvs} \
        --candidates ${candidates_txt} \
        --paralog-cutoffs ${paralog_cutoffs} \
        --config ${config_csv} \
        --paralog-competition-scope ${paralog_competition_scope} \
        --default-evalue ${default_evalue} \
        --output-matrix context_presence.tsv \
        --output-evalues context_presence.evalues.tsv
    """
}

workflow CONTEXT_SEARCH {
    take:
    candidates_txt  // path: candidates.txt (source_proteome::protein_id)
    self_hits       // [meta, paralog_cutoffs.tsv] — one per ingroup proteome (SEARCH.out.self_hits)
    ingroup_ch      // [meta, protein_fa] — to extract candidate + paralog sequences from
    near_in_ch      // [meta, protein_fa] — NEAR_INGROUP proteomes (may be empty)
    broad_out_ch    // [meta, protein_fa] — BROAD_OUTGROUP proteomes (may be empty)
    config_csv      // path to analysis CSV

    main:
    paralog_cutoffs_ch = self_hits.map { meta, tsv -> tsv }.collect().ifEmpty([])

    BUILD_CONTEXT_QUERY(candidates_txt, paralog_cutoffs_ch)

    EXTRACT_CANDIDATES(
        BUILD_CONTEXT_QUERY.out.query_list,
        ingroup_ch.map { m, fa -> fa }.collect(),
        config_csv,
        'context_query.fa'
    )
    context_query_fa = EXTRACT_CANDIDATES.out.candidates_fa.first()

    context_proteomes_ch = near_in_ch.mix(broad_out_ch)

    if (params.run_tool == 'phmmer') {
        CONTEXT_PHMMER_SEARCH(context_query_fa, context_proteomes_ch)
        raw_ch = CONTEXT_PHMMER_SEARCH.out
    }
    else if (params.run_tool == 'diamond') {
        db_ch = DIAMOND_MAKEDB(context_proteomes_ch)
        CONTEXT_DIAMOND_SEARCH(context_query_fa, db_ch)
        raw_ch = CONTEXT_DIAMOND_SEARCH.out
    }
    else if (params.run_tool == 'blast') {
        db_ch = BLAST_MAKEDB(context_proteomes_ch)
        CONTEXT_BLAST_SEARCH(context_query_fa, db_ch)
        raw_ch = CONTEXT_BLAST_SEARCH.out
    }
    else {
        error "Unknown --run_tool '${params.run_tool}': choose phmmer, diamond, or blast"
    }
    PARSE_HITS(raw_ch)
    hits_ch = PARSE_HITS.out.map { meta_pair, tsv -> tsv }.collect().ifEmpty([])

    CONTEXT_PRESENCE(
        hits_ch,
        candidates_txt,
        paralog_cutoffs_ch,
        config_csv,
        params.paralog_competition_scope,
        params.evalue
    )

    emit:
    matrix   = CONTEXT_PRESENCE.out.matrix
    evalues  = CONTEXT_PRESENCE.out.evalues
}
