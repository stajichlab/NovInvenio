// FAMILY_HMMSEARCH — scan the family HMM database against one proteome (ADR-0002 Q1/Q5).
// Run for EVERY proteome (ingroup + outgroup) so presence is one uniform, alignment-based
// criterion in both groups. `--domtblout` supplies per-domain HMM coordinates for the
// profile-coverage filter; a fixed `-Z` (params.hmm_Z) makes E-values comparable across
// proteomes of different size. The reporting threshold (-E) is intentionally loose — the
// real presence rule (E<params.hmm_presence_evalue AND coverage) is applied downstream by
// profile_to_matrix.py, so nothing faint is discarded before the coverage guard runs.
//
// No storeDir here: the output depends on the whole family HMM db (all ingroup + params),
// which the filename cannot encode, so caching is left to Nextflow -resume.
process FAMILY_HMMSEARCH {
    label 'hmmsearch'
    tag "${meta.id}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/${out_prefix}family_hmmsearch" }, mode: 'copy'

    input:
    tuple val(meta), path(proteome_fa)
    path(profiles_hmm)
    val(out_prefix)    // '' (novelty) | 'loss_' — family_hmmsearch/ vs loss_family_hmmsearch/

    output:
    tuple val(meta), path("${meta.id}.family.domtblout"), emit: domtblout

    script:
    def zflag = params.hmm_Z ? "-Z ${params.hmm_Z}" : ""
    """
    if [ -s ${profiles_hmm} ]; then
        hmmsearch --cpu ${task.cpus} \
            --domtblout ${meta.id}.family.domtblout \
            --noali \
            -E ${params.hmm_report_evalue} ${zflag} \
            ${profiles_hmm} ${proteome_fa} > /dev/null
    else
        touch ${meta.id}.family.domtblout
    fi
    """
}
