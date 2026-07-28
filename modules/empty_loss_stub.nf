// EMPTY_LOSS_STUB — placeholder loss-direction artifacts for --cluster_tool
// novelty_discovery. Loss analysis (present in outgroup, absent from ingroup) is an
// explicitly deferred future extension for the two-phase novelty_discovery/novelty_screen
// pathway (see todo/novelty-discovery-screen.md "Loss Direction") — TARGET/DISC_OUT configs
// have no IN/OUT rows, so the generic LOSS_SEARCH/LOSS_CLUSTER/LOSS_VALIDATE/LOSS_ANNOTATE
// mirror would only ever see empty channels and never fire.
//
// REPORT's COLLATE_REPORTS (which writes view/<project>/report.html) requires all three
// reports, including losses.html, so MAKE_LOSSES_REPORT must still run — with zero rows.
// build_losses_payload() requires at least one proteome column from the config to be
// present in the matrix header (otherwise it raises), so the stub matrix's header is
// copied from the real novelty-direction annotated matrix (same TARGET/DISC_OUT columns),
// just with zero data rows.
process EMPTY_LOSS_STUB {
    label 'low_cpu'
    tag "empty_loss_stub"

    input:
    path(header_source_matrix)   // annotated novelty-direction matrix, for its column header

    output:
    path("loss_presence_matrix.function.tsv"), emit: matrix
    path("loss_tblastn_summary.tsv"),          emit: tblastn_summary
    path("empty_loss_cluster.tsv"),            emit: cluster_tsv

    shell:
    '''
    head -n1 !{header_source_matrix} > loss_presence_matrix.function.tsv
    printf 'protein_id\n' > loss_tblastn_summary.tsv
    : > empty_loss_cluster.tsv
    '''
}
