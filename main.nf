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
    if (params.cluster_tool !in ['pairwise', 'mmseqs']) error "ERROR: --cluster_tool must be pairwise or mmseqs (got: ${params.cluster_tool})"

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
            def meta = [
                id:      row.Short,
                group:   row.GROUP,
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
    else {
        SEARCH(ingroup_prot_ch, outgroup_prot_ch, file(params.config))
        novelty_matrix     = SEARCH.out.matrix
        novelty_candidates = SEARCH.out.candidates

        CLUSTER(novelty_candidates, ingroup_prot_ch, file(params.config), 'candidates.fa', 'clusters')
        cand_fa          = CLUSTER.out.candidates_fa
        cand_reps        = CLUSTER.out.representatives
        cand_cluster_tsv = CLUSTER.out.cluster_tsv
    }

    VALIDATE(cand_reps, outgroup_dna_ch, cand_cluster_tsv, 'tblastn_summary.tsv')

    ANNOTATE(cand_fa, novelty_matrix, pfam_abs, sprot_abs, morgs_abs, '')

    SUMMARIZE(ANNOTATE.out.annotated_matrix, VALIDATE.out.summary, cand_cluster_tsv, file(params.config))

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
    else {
        // See workflows/loss_search.nf for why this needs its own search direction.
        LOSS_SEARCH(ingroup_prot_ch, outgroup_prot_ch, file(params.config))
        loss_matrix = LOSS_SEARCH.out.matrix

        LOSS_CLUSTER(LOSS_SEARCH.out.candidates, outgroup_prot_ch, file(params.config), 'loss_candidates.fa', 'loss_clusters')
        loss_cand_fa          = LOSS_CLUSTER.out.candidates_fa
        loss_cand_reps        = LOSS_CLUSTER.out.representatives
        loss_cand_cluster_tsv = LOSS_CLUSTER.out.cluster_tsv
    }

    LOSS_VALIDATE(loss_cand_reps, ingroup_dna_ch, loss_cand_cluster_tsv, 'loss_tblastn_summary.tsv')

    LOSS_ANNOTATE(loss_cand_fa, loss_matrix, pfam_abs, sprot_abs, morgs_abs, 'loss_')

    REPORT(
        ANNOTATE.out.annotated_matrix,
        VALIDATE.out.summary,
        SUMMARIZE.out.novelties,
        cand_fa,
        cand_cluster_tsv,
        LOSS_ANNOTATE.out.annotated_matrix,
        LOSS_VALIDATE.out.summary,
        loss_cand_cluster_tsv,
        file(params.config)
    )
}
