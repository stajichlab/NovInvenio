#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Entry point for the pangenome cluster-profiling subworkflow
// (workflows/pangenome_profile.nf). Kept as its own top-level script,
// separate from main.nf's IN/OUT novelty/loss-search workflow, since
// pangenome profiling is a structurally different analysis (whole-species-
// set clustering + presence/frequency/co-occurrence, not a two-group
// novelty/loss screen) that shares this repo's conventions (samplesheet
// schema, conf/ profiles, container image) but not its per-run channel
// wiring. See NEXTFLOW_MIGRATION_NOTES.md (in NovInvenio_Investigations'
// Afumigatus_pangenome study directory) for the design record this
// implements, and workflows/pangenome_profile.nf's own header comment for
// the pipeline shape.
//
// Usage:
//   nextflow run pangenome.nf -profile slurm -c conf/ucr_hpcc_slurm.config \
//       --pangenome_samplesheet /path/to/config.csv \
//       --pangenome_data_dir /path/to/data_dir \
//       [--pangenome_captain_hmm /path/to/captain.hmm | \
//        --pangenome_captain_hmm_name NAME --pangenome_pfam_hmm /path/to/Pfam-A.hmm]
//
// Samplesheet schema (same shape as main.nf's --config, plus GFF3 is
// REQUIRED here -- pair classification's physical-linkage test has no
// genomic-position source without it):
//   GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup
// GROUP values are compared against --pangenome_ingroup_label/
// --pangenome_outgroup_label (default IN/OUT) -- every other GROUP value in
// the samplesheet is ignored by this entry point, so a samplesheet shared
// with main.nf's novelty/loss run (which may also carry
// DISCOVERY_TARGET/NEAR_INGROUP/etc. rows) works unchanged.

include { PANGENOME_PROFILE } from './workflows/pangenome_profile'

// Resolve a FASTA basename against data_dir, checking the flat layout first
// then the listed subdirectories -- same convention as main.nf's resolve_fa.
def resolve_fa(String data_dir_abs, String basename, List<String> subdirs) {
    if (!basename) return []
    def candidates = ([''] + subdirs).collect { sub ->
        file(sub ? "${data_dir_abs}/${sub}/${basename}" : "${data_dir_abs}/${basename}")
    }
    def hit = candidates.find { it.exists() }
    if (!hit) error "Cannot locate FASTA '${basename}' under ${data_dir_abs} (also tried subdirs: ${subdirs.join(', ')})"
    return hit
}

def print_help() {
    log.info """
    Pangenome cluster-profiling subworkflow
    ========================================
    Usage:
      nextflow run pangenome.nf -profile slurm -c conf/ucr_hpcc_slurm.config \\
          --pangenome_samplesheet /path/to/config.csv \\
          --pangenome_data_dir /path/to/data_dir \\
          [--pangenome_captain_hmm /path/to/captain.hmm | \\
           --pangenome_captain_hmm_name NAME --pangenome_pfam_hmm /path/to/Pfam-A.hmm]

    Required arguments:
      --pangenome_samplesheet   Path to the samplesheet CSV (same schema as
                                 main.nf's --config, GFF3 column required).
      --pangenome_data_dir      Directory containing the Protein/DNA/GFF3 files
                                 referenced by the samplesheet.

    Optional arguments:
      --pangenome_project             Name used for the outdir/storeDir namespace
                                       (defaults to the samplesheet basename --
                                       see the "outdir collision" note below).
      --pangenome_gff3_dir             Override for the GFF3 subdirectory
                                       (default: <pangenome_data_dir>/gff3).
      --pangenome_ingroup_label        GROUP value treated as ingroup (default: IN).
      --pangenome_outgroup_label       GROUP value treated as outgroup (default: OUT).
      --pangenome_cluster_backend      mmseqs (default, validated) or diamond
                                       (opt-in/experimental -- ID-fidelity checked but
                                       not yet benchmarked on real data, see docs/adr/0003).
      --pangenome_captain_hmm          Path to a pre-built captain-gene HMM.
      --pangenome_captain_hmm_name     Named captain-gene model + --pangenome_pfam_hmm
                                       to build one from Pfam-A.hmm.
      --pangenome_dereplicate          Dereplicate near-identical strains via mash
                                       (default: true).
      --pangenome_assign_clades        Assign clades via mash + scipy clustering
                                       (default: true).
      --pangenome_island_pfam_hmm      Path to Pfam-A.hmm -- enables the accessory-island
                                       + Pfam functional-enrichment step (off by default).
                                       Distinct from --pangenome_pfam_hmm, which is only
                                       used by the captain-by-name branch above.
      --pangenome_island_min_size      Minimum island size to report (default: 2).
      --pangenome_pfam_domain_evalue   hmmscan domain-level E-value cutoff (default: 1e-3).
      --pangenome_marker_names         Comma list of named marker searches (e.g.
                                       'captain,sm_backbone'), run via MARKER_HMMSEARCH.
      --pangenome_marker_hmm_paths     Parallel comma list of HMM paths for each named marker.
      --pangenome_marker_evalue        hmmsearch E-value cutoff for marker searches (default: 1e-5).
      --pangenome_accumulation_permutations  Random strain-order permutations for the
                                       rarefaction/accumulation curve (default: 20).
      --pangenome_accumulation_seed    RNG seed for the accumulation curve (default: 0).
      --help                           Show this message and exit.

    Note: --pangenome_project (or a derivable default) is required so that two
    different species' runs never collide on outdir/storeDir -- see DESIGN
    notes / CLAUDE.md for the provenance rules that motivate this.
    """.stripIndent()
}

workflow {
    if (params.help) {
        print_help()
        exit 0
    }
    if (!params.pangenome_samplesheet) error "ERROR: --pangenome_samplesheet <config_csv> is required"
    if (!params.pangenome_data_dir)    error "ERROR: --pangenome_data_dir <data_dir> is required"
    if (!file(params.pangenome_samplesheet).exists())    error "ERROR: --pangenome_samplesheet file not found: ${params.pangenome_samplesheet}"
    if (!file(params.pangenome_data_dir).isDirectory())  error "ERROR: --pangenome_data_dir is not a directory: ${params.pangenome_data_dir}"
    if (params.pangenome_cluster_backend !in ['mmseqs', 'diamond'])
        error "ERROR: --pangenome_cluster_backend must be mmseqs or diamond (got: ${params.pangenome_cluster_backend})"
    // Previously hard-blocked here: diamond's tier-1 branch had no equivalent
    // of mmseqs' restore_mmseqs_cluster_ids.py safety net, and family IDs are
    // load-bearing for every downstream table (bin/pangenome_build_presence_matrix.py's
    // ID contract), so a silent ID mismatch would have corrupted every
    // downstream table without erroring. CLUSTER_TIER1 (modules/pangenome/prefix_and_cluster.nf)
    // now runs bin/verify_diamond_cluster_ids.py against the raw diamond cluster.tsv
    // before anything else reads it, and fails loud on any id mismatch --
    // see docs/adr/0003-diamond-tier1-clustering-backend.md. That check has
    // been unit-tested against adversarial headers (multi-pipe, Short-prefixed
    // UniProt-style ids) but not yet against a real diamond binary run at
    // pangenome scale -- treat --pangenome_cluster_backend diamond as
    // opt-in/experimental until that real-data concordance check (the ADR's
    // open item) has been done.

    // Helpers.projectName(params) (lib/Helpers.groovy) falls back to the bare
    // literal 'output' when neither params.project nor params.config is set.
    // pangenome.nf never sets params.config, so every run would otherwise
    // collide on the same outdir/storeDir regardless of species -- derive a
    // default from the samplesheet basename (mirrors main.nf's --config
    // convention) unless the caller already passed --pangenome_project or
    // --project explicitly.
    if (!params.project) {
        if (params.pangenome_project) {
            params.project = params.pangenome_project
        } else {
            params.project = new File(params.pangenome_samplesheet.toString()).name.replaceFirst(/\.[^.]+$/, '')
            log.warn "No --pangenome_project (or --project) given -- defaulting params.project to " +
                      "'${params.project}' (derived from --pangenome_samplesheet basename). Pass " +
                      "--pangenome_project explicitly to control the outdir/storeDir namespace."
        }
    }

    def data_dir_abs = file(params.pangenome_data_dir).toAbsolutePath().toString()
    def gff3_dir_abs = params.pangenome_gff3_dir ?: "${data_dir_abs}/gff3"
    def in_label  = params.pangenome_ingroup_label
    def out_label = params.pangenome_outgroup_label

    samples_ch = Channel
        .fromPath(params.pangenome_samplesheet)
        .splitCsv(header: true)
        .filter { row -> (row.GROUP?.trim()) in [in_label, out_label] }
        .map { row ->
            def meta = [id: row.Short, group: row.GROUP.trim(), taxon: row.TaxonGroup]
            def protein_fa = resolve_fa(data_dir_abs, row.Protein, ['pep', 'proteins'])
            def dna_fa     = resolve_fa(data_dir_abs, row.DNA,     ['dna', 'genome', 'scaffolds'])
            // Validate the GFF3 exists NOW, at samplesheet-parsing time, not
            // deep inside a later GENE_POSITIONS process failure after
            // clustering has already run (review item 9). GFF3 is REQUIRED
            // here (see this file's header comment) -- unlike Protein/DNA,
            // there is no subdir search convention for it, just
            // gff3_dir_abs/<basename>, matching GENE_POSITIONS' own lookup.
            if (!row.GFF3?.trim()) {
                error "ERROR: samplesheet row for Short='${row.Short}' has no GFF3 value -- " +
                      "GFF3 is required for pangenome.nf (pair classification's physical-linkage " +
                      "test has no genomic-position source without it)"
            }
            def gff3_file = file("${gff3_dir_abs}/${row.GFF3.trim()}")
            if (!gff3_file.exists()) {
                error "ERROR: GFF3 file for Short='${row.Short}' not found: ${gff3_file} " +
                      "(looked under --pangenome_gff3_dir=${gff3_dir_abs})"
            }
            [meta, protein_fa, dna_fa]
        }

    ingroup_dna_ch = samples_ch
        .filter { meta, prot, dna -> meta.group == in_label }
        .map    { meta, prot, dna -> [meta, dna] }

    PANGENOME_PROFILE(
        file(params.pangenome_samplesheet),
        data_dir_abs,
        gff3_dir_abs,
        samples_ch,
        ingroup_dna_ch,
    )
}
