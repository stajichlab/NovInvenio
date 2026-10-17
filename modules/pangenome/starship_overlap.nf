// STARSHIP_OVERLAP / SM_OVERLAP -- optional overlap of precomputed Starship element
// calls and secondary-metabolite (SM) cluster calls with the pangenome islands.
//
// Neither process calls elements or clusters itself: the tables come from outside the
// pipeline (e.g. starfish element calls, a published element table, antiSMASH regions)
// and are matched to this run through gene_positions (coordinates) and tier-1 family
// membership. Schemas: lib/starship_tables.py; methods: paper/methods_pangenome.md
// section 3.25. Both are off unless their param is set (workflows/pangenome_profile.nf,
// step 9a); each runs without the other. A run without either param has no new process.
//
// Outputs (file names and columns follow NovInvenio_Investigations' analysis scripts, so
// NII's bin/sync_pangenome_report.py builds its starship_overlap.html page from them):
//   <pangenome>/starship_overlap/      STARSHIP_OVERLAP (bin/pangenome_starship_overlap.py)
//   <pangenome>/starship_overlap/sm/   SM_OVERLAP       (bin/pangenome_sm_overlap.py)
//
// Outputs are declared as files (starship_overlap/*, sm/*), not as directories: with
// both params set, publishing the starship_overlap/ directory into a target that already
// held SM_OVERLAP's starship_overlap/sm/ was skipped without an error in a local test
// run (Nextflow 26.04.6, 2026-10-17). File outputs publish one by one and do not clash.
//
// Inputs that may be EMPTY_EVALUES_STUB outputs (0-byte empty_evalues.tsv: rescue
// positions when the rescue pass is off, the captain tblout when no captain HMM is set)
// stage under unique names; the scripts read a 0-byte file as "not given".
process STARSHIP_OVERLAP {
    label 'low_cpu'
    tag "starship_overlap"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(starship_table, stageAs: 'inputs/starship/*')
    path(frequency_table)
    path(significant_islands)
    path(cluster_tsv)
    path(gene_positions)
    path(family_positions)
    path(rescue_positions, stageAs: 'inputs/rescue_positions.tsv')
    path(captain_tblout, stageAs: 'inputs/captain.tblout')
    path(samplesheet, stageAs: 'inputs/samplesheet.csv')

    output:
    path("starship_overlap/*"), emit: files

    script:
    """
    pangenome_starship_overlap.py \
        --frequency_table ${frequency_table} \
        --significant_islands ${significant_islands} \
        --cluster_tsv ${cluster_tsv} \
        --gene_positions ${gene_positions} \
        --family_positions ${family_positions} \
        --rescue_positions ${rescue_positions} \
        --starship_annotation ${starship_table} \
        --captain_tblout ${captain_tblout} \
        --samplesheet ${samplesheet} \
        --id_sep '${params.pangenome_id_sep}' \
        --n_perm ${params.pangenome_starship_permutations} \
        --seed ${params.pangenome_starship_seed} \
        --run_label '${Helpers.projectName(params)}' \
        --out_dir starship_overlap
    """
}

process SM_OVERLAP {
    label 'low_cpu'
    tag "sm_overlap"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome/starship_overlap" }, mode: 'copy'

    input:
    path(sm_table, stageAs: 'inputs/sm/*')
    path(starship_table, stageAs: 'inputs/starship/*')   // [] when --pangenome_starship_annotation is unset
    path(frequency_table)
    path(significant_islands)
    path(cluster_tsv)
    path(gene_positions)
    path(captain_tblout, stageAs: 'inputs/captain.tblout')
    path(samplesheet, stageAs: 'inputs/samplesheet.csv')

    output:
    path("sm/*"), emit: files

    script:
    def starship_arg = starship_table ? "--starship_annotation ${starship_table}" : ''
    """
    pangenome_sm_overlap.py \
        --sm_clusters ${sm_table} \
        ${starship_arg} \
        --frequency_table ${frequency_table} \
        --significant_islands ${significant_islands} \
        --cluster_tsv ${cluster_tsv} \
        --gene_positions ${gene_positions} \
        --captain_tblout ${captain_tblout} \
        --samplesheet ${samplesheet} \
        --id_sep '${params.pangenome_id_sep}' \
        --n_perm ${params.pangenome_starship_permutations} \
        --seed ${params.pangenome_starship_seed} \
        --out_dir sm
    """
}
