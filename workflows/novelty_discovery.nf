nextflow.enable.dsl=2

// NOVELTY_DISCOVERY — stub for the two-phase targeted novelty pipeline.
// Ticket #26 will implement the real workflow. This stub exists so that
// main.nf can branch on --cluster_tool novelty_discovery and validate the
// CLI plumbing in Ticket #25.
workflow NOVELTY_DISCOVERY {
    take:
    target_ch      // [meta, protein_fa] — TARGET genomes
    disc_out_ch    // [meta, protein_fa] — discovery outgroup proteomes
    disc_out_dna   // [meta, dna_fa] — discovery outgroup genomes
    config_csv     // path to analysis CSV

    main:
    error "NOVELTY_DISCOVERY is not yet implemented (see ticket #26)"

    emit:
    matrix         = []
    candidates     = []
    cand_fa        = []
    cand_reps      = []
    cand_cluster_tsv = []
    calibrated_hmms  = []
}
