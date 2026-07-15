nextflow.enable.dsl=2

// Produce per-species novelty candidate files (novelties.<SHORT>.tsv).
// Combines the annotated presence matrix with TBLASTN summary to apply the
// final filter: present in ingroup, absent from all outgroup proteomes AND
// absent from all outgroup genomes.

workflow SUMMARIZE {
    take:
    annotated_matrix   // path: presence_matrix.function.tsv
    tblastn_summary    // path: tblastn_summary.tsv (protein × genome hit matrix)
    config_csv         // path: analysis CSV

    main:
    MAKE_NOVELTIES(annotated_matrix, tblastn_summary, config_csv)

    emit:
    novelties = MAKE_NOVELTIES.out.novelties
}

process MAKE_NOVELTIES {
    label 'low_cpu'
    publishDir { "${params.outdir}/${Helpers.projectName(params)}" }, mode: 'copy'

    input:
    path(annotated_matrix)
    path(tblastn_summary)
    path(config_csv)

    output:
    path("novelties.*.tsv"), emit: novelties

    script:
    """
    make_novelties.py \
        --matrix ${annotated_matrix} \
        --tblastn_summary ${tblastn_summary} \
        --config ${config_csv} \
        --ingroup_min ${params.ingroup_min_frac} \
        --output_dir . \
	--skip_tblastn_filter
    """
}
