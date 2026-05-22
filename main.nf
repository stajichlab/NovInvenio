#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { SEARCH   } from './workflows/search'
include { CLUSTER  } from './workflows/cluster'
include { VALIDATE } from './workflows/validate'

workflow {
    // Derive project name from config filename when not explicitly set
    def project = params.project ?: file(params.config).baseName
    params.project = project

    if (!params.config)   error "ERROR: --config <analysis_csv> is required"
    if (!params.data_dir) error "ERROR: --data_dir <directory> is required"

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
            def protein_fa = file("${params.data_dir}/${row.Protein}")
            def dna_fa     = row.DNA ? file("${params.data_dir}/${row.DNA}") : []
            [ meta, protein_fa, dna_fa ]
        }

    ingroup_prot_ch   = samples_ch.filter { meta, prot, dna -> meta.group == 'IN' }
                                   .map    { meta, prot, dna -> [ meta, prot ] }
    outgroup_prot_ch  = samples_ch.filter { meta, prot, dna -> meta.group == 'OUT' }
                                   .map    { meta, prot, dna -> [ meta, prot ] }
    outgroup_dna_ch   = samples_ch.filter { meta, prot, dna -> meta.group == 'OUT' && dna }
                                   .map    { meta, prot, dna -> [ meta, dna ] }

    SEARCH(ingroup_prot_ch, outgroup_prot_ch, file(params.config))

    CLUSTER(SEARCH.out.candidates, ingroup_prot_ch)

    VALIDATE(CLUSTER.out.representatives, outgroup_dna_ch)
}
