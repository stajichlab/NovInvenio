#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { SEARCH   } from './workflows/search'
include { CLUSTER  } from './workflows/cluster'
include { VALIDATE } from './workflows/validate'
include { ANNOTATE } from './workflows/annotate'
include { SUMMARIZE } from './workflows/summarize'

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

    // Resolve DB and config paths to absolute so they work inside Nextflow work directories
    if (params.pfam_hmm)         params.pfam_hmm         = file(params.pfam_hmm).toAbsolutePath().toString()
    if (params.swissprot_dmnd)   params.swissprot_dmnd   = file(params.swissprot_dmnd).toAbsolutePath().toString()
    if (params.modelorgs_config) params.modelorgs_config = file(params.modelorgs_config).toAbsolutePath().toString()

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

    SEARCH(ingroup_prot_ch, outgroup_prot_ch, file(params.config))

    CLUSTER(SEARCH.out.candidates, ingroup_prot_ch, file(params.config))

    VALIDATE(CLUSTER.out.representatives, outgroup_dna_ch, CLUSTER.out.cluster_tsv)

    ANNOTATE(CLUSTER.out.candidates_fa, SEARCH.out.matrix)

    SUMMARIZE(ANNOTATE.out.annotated_matrix, VALIDATE.out.summary, file(params.config))
}
