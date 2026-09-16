// CAPTAIN_HMMSEARCH — optional "mobile-element captain gene" hmmsearch used
// by pair classification's physical-linkage evidence (e.g. the fungal
// Starship-transposon DUF3435 captain gene in the originating study).
// Deliberately OPTIONAL and fungal/Starship-agnostic in its mechanics: any
// hmmsearch --tblout of ANY marker HMM against the all-strains proteome
// works the same way downstream (pangenome_pair_classification.py just
// checks "is a hit from this tblout near this family pair"). A study with no
// such marker (most non-fungal, non-Starship-relevant species) skips this
// process entirely -- see workflows/pangenome_profile.nf, which emits an
// empty stub file (EMPTY_EVALUES_STUB, reused) instead when
// params.pangenome_captain_hmm/--pangenome_captain_hmm_name are both unset.
//
// Two ways to supply the marker HMM:
//   --pangenome_captain_hmm <path.hmm>          a single-profile HMM file, used as-is.
//   --pangenome_captain_hmm_name <NAME> \
//     --pangenome_pfam_hmm <Pfam-A.hmm>         hmmfetch <NAME> out of a larger Pfam db
//                                                first (NAME, not accession -- e.g. Pfam's
//                                                PF11001 entry is fetched by its NAME field,
//                                                not the PF##### accession).
//
// This is a deliberately coarse proxy, not element-level annotation: it only
// asks "is a captain-family hmmsearch hit within k genes of this co-occurring
// pair, in any strain" (lib/pangenome_synteny.py's linkage_fraction, reused).
// For fungal Starship biology specifically, `starfish`
// (https://github.com/egluckthaler/starfish; Gluck-Thaler & Vogan 2024,
// NAR) is the purpose-built tool for real element-level work this
// subworkflow does not attempt: it re-annotates captains against an 11-family
// YR HMM profile set (metaeuk + hmmsearch), calls actual element BOUNDARIES
// (direct-repeat/TIR-supported flank/insert/extend confidence tiers, not a
// fixed gene-count window), clusters captains into orthology-based
// navis/haplotype groups (Orthofinder or mmseqs + mcl), and builds a
// genome-level element presence/absence matrix via orthologous-flank
// dereplication (searching up to 600kb around each boundary) that calls
// "empty"/"fragmented" haplotypes explicitly -- a materially different,
// finer-grained analysis than this subworkflow's family-level tier-1
// clustering + a rank-window adjacency heuristic. A study specifically
// investigating Starships should run starfish separately on the same genome
// set and cross-reference its element/navis/haplotype calls against this
// subworkflow's `unexplained_physical`-classified pairs (candidates for an
// undetected/unannotated element) rather than treat this module's captain
// hmmsearch as a substitute for it.
process HMMFETCH_CAPTAIN {
    label 'low_cpu'
    tag "hmmfetch_captain"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(pfam_hmm)
    val(hmm_name)

    output:
    path("captain.hmm"), emit: hmm

    script:
    """
    hmmfetch ${pfam_hmm} ${hmm_name} > captain.hmm
    """
}

process CAPTAIN_HMMSEARCH {
    label 'med_cpu'
    tag "captain_hmmsearch"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(captain_hmm)
    path(all_strains_fa)

    output:
    path("captain_vs_study.tblout"), emit: tblout

    script:
    """
    hmmsearch --tblout captain_vs_study.tblout -E ${params.pangenome_captain_hmm_evalue} \
        --cpu ${task.cpus} ${captain_hmm} ${all_strains_fa} > /dev/null
    """
}
