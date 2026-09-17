// PREFIX_PROTEOME / PREFIX_GENOME — Short-prefix every strain's protein and
// genome FASTA headers (><Short><id_sep><original_id>) before clustering,
// per the pangenome-profiling subworkflow's ID convention (see
// bin/pangenome_build_presence_matrix.py's module docstring). Was a manual
// awk one-liner in the originating study; pangenome_prefix_fasta.py gives it
// a real, tested process body (NEXTFLOW_MIGRATION_NOTES.md section A item 1).
process PREFIX_PROTEOME {
    label 'low_cpu'
    tag "${meta.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    tuple val(meta), path(protein_fa)

    output:
    tuple val(meta), path("${meta.id}.prefixed.pep.fa"), emit: fasta

    script:
    """
    pangenome_prefix_fasta.py --short ${meta.id} --input ${protein_fa} \
        --id_sep '${params.pangenome_id_sep}' --output ${meta.id}.prefixed.pep.fa
    """
}

// Same prefixing, applied to genome/DNA FASTAs -- consumed by the rescue-pass
// per-strain BLAST database build (modules/pangenome/rescue.nf), so a tblastn
// subject ID carries the strain prefix the same way a cluster member ID does.
process PREFIX_GENOME {
    label 'low_cpu'
    tag "${meta.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    tuple val(meta), path(genome_fa)

    output:
    tuple val(meta), path("${meta.id}.prefixed.dna.fa"), emit: fasta

    script:
    """
    pangenome_prefix_fasta.py --short ${meta.id} --input ${genome_fa} \
        --id_sep '${params.pangenome_id_sep}' --output ${meta.id}.prefixed.dna.fa
    """
}

// Concatenate every strain's prefixed proteome into one clustering input FASTA.
process CONCAT_PROTEOMES {
    label 'low_cpu'
    tag "concat_proteomes"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(prefixed_fastas)

    output:
    path("all_strains.pep.fa"), emit: fasta

    script:
    """
    cat ${prefixed_fastas} > all_strains.pep.fa
    """
}

// Tier-1 clustering (the allele/ortholog unit every downstream step relies
// on). Dispatches on params.pangenome_cluster_backend ('mmseqs' | 'diamond')
// -- pangenome_cluster_backend.py itself is the tested two-tier implementation
// ported from the originating study; only tier 1 is wired into this
// subworkflow today (tier 2 -- a loose superfamily label, never used for
// presence/frequency calls -- is left as a documented future extension, see
// this repo's CHANGES.md/README for how to add a CLUSTER_TIER2 process
// alongside this one using the same script's `*-tier2` subcommands).
process CLUSTER_TIER1 {
    label 'high_cpu'
    tag "cluster_tier1"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome/cluster" }, mode: 'copy'

    input:
    path(all_strains_fa)

    output:
    path("tier1_cluster.tsv"),   emit: cluster_tsv
    path("tier1_rep_seq.fasta"), emit: rep_fasta

    script:
    def backend = params.pangenome_cluster_backend
    if (backend == 'mmseqs')
        """
        pangenome_cluster_backend.py mmseqs-tier1 --fasta ${all_strains_fa} --out_prefix tier1 \
            --min_seq_id ${params.pangenome_tier1_min_id} --cov ${params.pangenome_tier1_cov} \
            --threads ${task.cpus}
        """
    else if (backend == 'diamond')
        // Unlike mmseqs (which needs bin/restore_mmseqs_cluster_ids.py to
        // undo its own header-collapsing -- see that script's docstring and
        // lib/fasta.py::mmseqs_id()), diamond cluster is expected to keep
        // the whole first-token header verbatim in its own *_cluster.tsv --
        // no collapsing, so no restoration needed. That expectation is
        // checked here, not assumed: verify_diamond_cluster_ids.py fails
        // loud if any id in the raw cluster.tsv doesn't match a real header
        // in ${all_strains_fa} verbatim, before anything downstream reads
        // it. See docs/adr/0003-diamond-tier1-clustering-backend.md -- this
        // closes the gap that previously kept
        // --pangenome_cluster_backend diamond hard-blocked in pangenome.nf.
        // Still not validated against a real multi-strain run (no diamond
        // binary in the environment used to build this check); the
        // mmseqs-vs-diamond concordance benchmark called for in that ADR is
        // a separate, real-data follow-up.
        """
        pangenome_cluster_backend.py diamond-tier1 --fasta ${all_strains_fa} --out_prefix tier1 \
            --approx_id ${(params.pangenome_tier1_min_id as double) * 100} \
            --member_cover ${(params.pangenome_tier1_cov as double) * 100} \
            --threads ${task.cpus}
        verify_diamond_cluster_ids.py --input-fasta ${all_strains_fa} --cluster-tsv tier1_cluster.tsv
        mv tier1_cluster.tsv tier1_cluster.tsv.raw
        # diamond cluster's TSV has no separate *_rep_seq.fasta; extract representatives
        # (col 1, unique) from the input FASTA to give both backends the same output contract.
        cut -f1 tier1_cluster.tsv.raw | sort -u > tier1_reps.txt
        awk 'BEGIN{while((getline l < "tier1_reps.txt")>0) want[l]=1}
             /^>/{id=substr(\$1,2); keep=(id in want)} keep' ${all_strains_fa} > tier1_rep_seq.fasta
        mv tier1_cluster.tsv.raw tier1_cluster.tsv
        """
    else
        error "pangenome_cluster_backend must be 'mmseqs' or 'diamond' (got: ${backend})"
}
