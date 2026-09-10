nextflow.enable.dsl=2

// Render the interactive HTML reports:
//   - novelties.html — ingroup-specific candidates (annotated matrix + TBLASTN
//     summary + per-species novelties tables)
//   - core.html — near-universally conserved genes (annotated matrix +
//     cluster_tsv only, no new search/annotation step)
//   - losses.html — candidate lineage-specific gene losses (loss-search
//     annotated matrix + loss TBLASTN summary + loss cluster_tsv, all produced
//     by the LOSS_SEARCH/LOSS_CLUSTER/LOSS_VALIDATE/LOSS_ANNOTATE mirror
//     pipeline in main.nf)
// All three pages are self-contained — no network access is needed to open them.

workflow REPORT {
    take:
    annotated_matrix       // path: presence_matrix.function.tsv
    tblastn_summary        // path: tblastn_summary.tsv
    novelties               // path(s): novelties.<SHORT>.tsv
    candidates_fa           // path: candidates.fa (protein sequences)
    cluster_tsv              // path: mmseqs *_cluster.tsv (rep -> member) for gene-family grouping
    evalues                  // path: presence_matrix.evalues.tsv sidecar (issue #44) --
                              //   report-only search-hit e-value evidence, same shape as the
                              //   matrix. May be an empty stub (EMPTY_EVALUES_STUB) when the
                              //   producer pathway doesn't track e-values yet.
    targets                  // path: presence_matrix.targets.tsv sidecar -- target protein ID
                              //   per qualifying hit, aligned with `evalues`. --cluster_tool
                              //   pairwise only for now; empty stub otherwise (same
                              //   EMPTY_EVALUES_STUB reused -- read_targets()==read_evalues(),
                              //   same "empty means no evidence" contract).
    descriptions              // path: bin/extract_protein_descriptions.py output TSV, used to
                              //   resolve `targets`' IDs into a display name. Same empty-stub
                              //   convention when unavailable.
    context_matrix           // path: context_presence.tsv (issue #48) -- NEAR_INGROUP/
                              //   BROAD_OUTGROUP presence for the candidate list only,
                              //   report-only, never used for novelty calling. Empty stub
                              //   outside --cluster_tool pairwise or when the config has no
                              //   NEAR_INGROUP/BROAD_OUTGROUP rows.
    context_evalues          // path: context_presence.evalues.tsv sidecar for context_matrix
    loss_annotated_matrix    // path: loss_presence_matrix.function.tsv
    loss_tblastn_summary     // path: loss_tblastn_summary.tsv
    loss_cluster_tsv         // path: loss mmseqs *_cluster.tsv
    config_csv               // path: analysis CSV
    data_dir                  // val: absolute path to --data_dir -- report-only, used to
                               //   resolve each species' optional GFF3 config-CSV column for
                               //   the reports' chrom/start columns (not staged/channeled like
                               //   the FASTAs; the standalone bin/make_*.py scripts re-resolve
                               //   it themselves, see lib/gff3_genes.resolve_gff3_paths)

    main:
    // results/ copy: offline, file://-safe, never carries the TBLASTN
    // alignment popup (issue #74's ALIGNMENT_POPUP_JS).
    MAKE_REPORT(annotated_matrix, tblastn_summary, novelties, candidates_fa, cluster_tsv,
               evalues, targets, descriptions, context_matrix, context_evalues, config_csv, data_dir)
    MAKE_CORE_REPORT(annotated_matrix, cluster_tsv, config_csv, data_dir)
    MAKE_LOSSES_REPORT(loss_annotated_matrix, loss_tblastn_summary, loss_cluster_tsv, config_csv, data_dir)

    // docs/ copy: fetch-capable (GitHub Pages), regenerated (not copied) with
    // --online so it includes the alignment popup wired to its own
    // alignments/loss_alignments shards (issue #75). core.html has no TBLASTN
    // column at all (MAKE_CORE_REPORT never takes a tblastn_summary input), so
    // it never diverges between the two copies and is reused as-is.
    MAKE_REPORT_ONLINE(annotated_matrix, tblastn_summary, novelties, candidates_fa, cluster_tsv,
                       evalues, targets, descriptions, context_matrix, context_evalues, config_csv, data_dir)
    MAKE_LOSSES_REPORT_ONLINE(loss_annotated_matrix, loss_tblastn_summary, loss_cluster_tsv, config_csv, data_dir)

    // Publication-quality PDF summary (static figures) alongside the interactive HTML.
    MAKE_PDF_REPORT(annotated_matrix, tblastn_summary, cluster_tsv,
                    loss_annotated_matrix, loss_tblastn_summary, loss_cluster_tsv, config_csv)

    // Final step: gather the three (docs/-flavored) reports under docs/<project>/
    // with a report.html landing page describing the run (ingroup/outgroup, tool,
    // thresholds).
    COLLATE_REPORTS(
        MAKE_REPORT_ONLINE.out.report,
        MAKE_CORE_REPORT.out.report,
        MAKE_LOSSES_REPORT_ONLINE.out.report,
        config_csv,
    )

    emit:
    report        = MAKE_REPORT.out.report
    core_report   = MAKE_CORE_REPORT.out.report
    losses_report = MAKE_LOSSES_REPORT.out.report
    pdf_report    = MAKE_PDF_REPORT.out.pdf
    index_report  = COLLATE_REPORTS.out.index
}

process MAKE_PDF_REPORT {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${Helpers.docsDir(params, workflow.launchDir)}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(annotated_matrix)
    path(tblastn_summary)
    path(cluster_tsv)
    path(loss_annotated_matrix)
    path(loss_tblastn_summary)
    path(loss_cluster_tsv)
    path(config_csv)

    output:
    path("summary.pdf"), emit: pdf

    // Gate: skip the PDF step when --pdf_report false (e.g. no matplotlib available).
    when:
    params.pdf_report != false

    script:
    """
    make_pdf_report.py \
        --matrix ${annotated_matrix} \
        --config ${config_csv} \
        --tblastn_summary ${tblastn_summary} \
        --cluster_tsv ${cluster_tsv} \
        --project ${Helpers.projectName(params)} \
        --ingroup_min_frac ${params.ingroup_min_frac} \
        --loss_matrix ${loss_annotated_matrix} \
        --loss_tblastn_summary ${loss_tblastn_summary} \
        --loss_cluster_tsv ${loss_cluster_tsv} \
        --outgroup_min_frac ${params.outgroup_min_frac} \
        --loss_ingroup_max_frac ${params.loss_ingroup_max_frac} \
        --output summary.pdf
    """
}

process MAKE_REPORT {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(annotated_matrix)
    path(tblastn_summary)
    path(novelties)
    path(candidates_fa)
    path(cluster_tsv)
    // evalues/targets/context_matrix/context_evalues can all be the same stub process's
    // literal "empty_evalues.tsv" output when unused (EMPTY_EVALUES_STUB /
    // EMPTY_CONTEXT_MATRIX_STUB / EMPTY_CONTEXT_EVALUES_STUB) -- stageAs disambiguates
    // so Nextflow doesn't reject them as an input file name collision. descriptions can
    // likewise be the same empty stub file -- read_descriptions() treats an empty file
    // as "no descriptions" the same way read_evalues() does.
    path(evalues, stageAs: 'evalues.tsv')
    path(targets, stageAs: 'targets.tsv')
    path(descriptions, stageAs: 'descriptions.tsv')
    path(context_matrix, stageAs: 'context_matrix.tsv')
    path(context_evalues, stageAs: 'context_evalues.tsv')
    path(config_csv)
    val(data_dir)

    output:
    path("novelties.html"), emit: report

    script:
    // evalues/targets/descriptions/context_matrix/context_evalues may be empty stub
    // files -- make_report.py's read_evalues()/read_targets()/read_descriptions()/
    // read_context() treat a missing/empty/header-only file as "no evidence available"
    // (issue #44/#48).
    """
    make_report.py \
        --matrix ${annotated_matrix} \
        --config ${config_csv} \
        --tblastn_summary ${tblastn_summary} \
        --novelties ${novelties} \
        --candidates_fa ${candidates_fa} \
        --cluster_tsv ${cluster_tsv} \
        --evalues ${evalues} \
        --targets ${targets} \
        --descriptions ${descriptions} \
        --context_matrix ${context_matrix} \
        --context_evalues ${context_evalues} \
        --project ${Helpers.projectName(params)} \
        --ingroup_min_frac ${params.ingroup_min_frac} \
        --sequences ${params.report_sequences} \
        --data_dir ${data_dir} \
        --output novelties.html
    """
}

// docs/-only twin of MAKE_REPORT: identical inputs, --online flag included.
// No publishDir of its own -- COLLATE_REPORTS is what actually places this
// into docs/<project>/, exactly like MAKE_REPORT's output only reaches
// results/<project>/ through its own publishDir. Publishing this here too
// would just double-write the same bytes to the same docs/ path.
process MAKE_REPORT_ONLINE {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(annotated_matrix)
    path(tblastn_summary)
    path(novelties)
    path(candidates_fa)
    path(cluster_tsv)
    path(evalues, stageAs: 'evalues.tsv')
    path(targets, stageAs: 'targets.tsv')
    path(descriptions, stageAs: 'descriptions.tsv')
    path(context_matrix, stageAs: 'context_matrix.tsv')
    path(context_evalues, stageAs: 'context_evalues.tsv')
    path(config_csv)
    val(data_dir)

    output:
    path("novelties.html"), emit: report

    script:
    """
    make_report.py \
        --matrix ${annotated_matrix} \
        --config ${config_csv} \
        --tblastn_summary ${tblastn_summary} \
        --novelties ${novelties} \
        --candidates_fa ${candidates_fa} \
        --cluster_tsv ${cluster_tsv} \
        --evalues ${evalues} \
        --targets ${targets} \
        --descriptions ${descriptions} \
        --context_matrix ${context_matrix} \
        --context_evalues ${context_evalues} \
        --project ${Helpers.projectName(params)} \
        --ingroup_min_frac ${params.ingroup_min_frac} \
        --sequences ${params.report_sequences} \
        --data_dir ${data_dir} \
        --online \
        --output novelties.html
    """
}

process MAKE_CORE_REPORT {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(annotated_matrix)
    path(cluster_tsv)
    path(config_csv)
    val(data_dir)

    output:
    path("core.html"), emit: report

    script:
    """
    make_core_report.py \
        --matrix ${annotated_matrix} \
        --config ${config_csv} \
        --cluster_tsv ${cluster_tsv} \
        --project ${Helpers.projectName(params)} \
        --core_min_frac ${params.core_min_frac} \
        --data_dir ${data_dir} \
        --output core.html
    """
}

process MAKE_LOSSES_REPORT {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(loss_annotated_matrix)
    path(loss_tblastn_summary)
    path(loss_cluster_tsv)
    path(config_csv)
    val(data_dir)

    output:
    path("losses.html"), emit: report

    script:
    """
    make_losses_report.py \
        --matrix ${loss_annotated_matrix} \
        --config ${config_csv} \
        --tblastn_summary ${loss_tblastn_summary} \
        --cluster_tsv ${loss_cluster_tsv} \
        --data_dir ${data_dir} \
        --project ${Helpers.projectName(params)} \
        --outgroup_min_frac ${params.outgroup_min_frac} \
        --loss_ingroup_max_frac ${params.loss_ingroup_max_frac} \
        --output losses.html
    """
}

// docs/-only twin of MAKE_LOSSES_REPORT -- see MAKE_REPORT_ONLINE's comment.
process MAKE_LOSSES_REPORT_ONLINE {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(loss_annotated_matrix)
    path(loss_tblastn_summary)
    path(loss_cluster_tsv)
    path(config_csv)
    val(data_dir)

    output:
    path("losses.html"), emit: report

    script:
    """
    make_losses_report.py \
        --matrix ${loss_annotated_matrix} \
        --config ${config_csv} \
        --tblastn_summary ${loss_tblastn_summary} \
        --cluster_tsv ${loss_cluster_tsv} \
        --data_dir ${data_dir} \
        --project ${Helpers.projectName(params)} \
        --outgroup_min_frac ${params.outgroup_min_frac} \
        --loss_ingroup_max_frac ${params.loss_ingroup_max_frac} \
        --online \
        --output losses.html
    """
}

// Collate the three reports into docs/<project>/ and generate a report.html
// landing page. The three HTML files are re-published here (copied) so the whole
// result set lives in one shareable folder alongside its index.
process COLLATE_REPORTS {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${Helpers.docsDir(params, workflow.launchDir)}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(novelties_html)
    path(core_html)
    path(losses_html)
    path(config_csv)

    output:
    path("report.html"),     emit: index
    path("alignment.html"),  emit: alignment
    path(novelties_html),    emit: novelties
    path(core_html),         emit: core
    path(losses_html),       emit: losses

    script:
    """
    make_index_report.py \
        --config ${config_csv} \
        --project ${Helpers.projectName(params)} \
        --reports_dir . \
        --run_tool ${params.run_tool} \
        --ingroup_min_frac ${params.ingroup_min_frac} \
        --outgroup_min_frac ${params.outgroup_min_frac} \
        --loss_ingroup_max_frac ${params.loss_ingroup_max_frac} \
        --core_min_frac ${params.core_min_frac} \
        --output report.html
    """
}
