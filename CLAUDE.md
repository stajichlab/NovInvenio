# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NovInvenio identifies lineage-specific genes: proteins present in a defined ingroup (≥N% of members) but absent from all outgroup proteomes. It uses pairwise protein searches (phmmer/diamond/blast), self-vs-self paralog lookup, mmseqs2 clustering, TBLASTN validation against outgroup genomes, functional annotation (Pfam + SwissProt), and model-organism gene-name lookup to produce per-species novelty candidate tables.

## Commands

```bash
# Set up environment
pixi install

# Run the full pipeline
nextflow run main.nf \
    --config configs/pezizo4_asco.csv \
    --data_dir /path/to/fastas \
    --run_tool phmmer

# Run tests
pytest tests/ -v

# Run a single test file
pytest tests/test_parse_hits.py -v

# Lint
pixi run lint
```

## Repository Layout

```
NovInvenio/
├── main.nf                        # Nextflow entry point; imports and calls all workflows
├── nextflow.config                # params{} defaults and executor/resource settings
├── pixi.toml                      # [workspace] table — tool dependencies
├── modules/
│   ├── phmmer.nf                  # PHMMER_SEARCH — phmmer pairwise search, storeDir
│   ├── diamond.nf                 # DIAMOND_SEARCH — diamond blastp pairwise, storeDir
│   ├── blast.nf                   # BLAST_SEARCH — blastp pairwise, storeDir
│   ├── self_search.nf             # PHMMER_SELF / DIAMOND_SELF / BLAST_SELF — self-vs-self, storeDir
│   ├── parse_hits.nf              # PARSE_HITS — normalises raw hits → common TSV
│   ├── parse_self_hits.nf         # PARSE_SELF_HITS — extracts per-gene paralog cutoffs
│   ├── build_presence_matrix.nf   # BUILD_PRESENCE_MATRIX — presence/absence matrix + candidates.txt
│   ├── mmseqs_cluster.nf          # MMSEQS_CLUSTER — mmseqs2 easy-cluster of candidates
│   ├── tblastn.nf                 # TBLASTN — translated search vs outgroup genomes
│   ├── hmmbuild.nf                # HMMBUILD (imported but not yet wired into cluster workflow)
│   └── hmmsearch.nf               # HMMSEARCH (imported but not yet wired into cluster workflow)
├── workflows/
│   ├── search.nf                  # SEARCH — pairwise search + self-hits + presence matrix (ingroup query)
│   ├── loss_search.nf             # LOSS_SEARCH — same, outgroup query (loss-search direction)
│   ├── cluster.nf                 # CLUSTER — extract candidates + mmseqs2 clustering (reused as LOSS_CLUSTER)
│   ├── validate.nf                # VALIDATE — TBLASTN + SUMMARIZE_TBLASTN (reused as LOSS_VALIDATE)
│   ├── annotate.nf                # ANNOTATE — Pfam hmmscan + SwissProt diamond + matrix merge (reused as LOSS_ANNOTATE)
│   ├── summarize.nf               # SUMMARIZE — MAKE_NOVELTIES per ingroup species
│   └── report.nf                  # REPORT — MAKE_REPORT (novelties.html), MAKE_CORE_REPORT (core.html),
│                                   #   MAKE_LOSSES_REPORT (losses.html), COLLATE_REPORTS (view/<project>/report.html + copies)
├── bin/                           # Executable Python helpers; NF adds bin/ to PATH automatically
│   ├── parse_hits.py              # Normalise phmmer/diamond/blast → query_id, target_id, evalue, bitscore, query_proteome, target_proteome TSV
│   ├── parse_self_hits.py         # Self-hit rank-2 → per-gene paralog_cutoffs.tsv
│   ├── build_presence_matrix.py   # Paralog-aware matrix construction + candidates.txt; --query-group IN|OUT
│   ├── extract_candidates.py      # Pull candidate sequences from the given proteome FASTAs (ingroup or outgroup)
│   ├── summarize_tblastn.py       # Aggregate per-genome TBLASTN TSVs → protein × genome matrix
│   ├── annotate_presence_matrix.py# Add gene_name / Pfam / SwissProt columns
│   ├── make_novelties.py          # Per-species novelties.<SHORT>.tsv (with --skip_tblastn_filter option)
│   ├── make_report.py             # Self-contained interactive novelties.html
│   ├── make_core_report.py        # Self-contained interactive core.html (near-universal genes)
│   ├── make_losses_report.py      # Self-contained interactive losses.html (candidate gene losses)
│   ├── make_index_report.py       # Self-contained view/<project>/report.html landing page linking the three reports
│   ├── make_config.py             # Generate a run config CSV from configs/samples.csv by taxon/lineage matching
│   ├── convert_bfd_samples.py     # Convert a BFD Fungi_BFD samples.csv into a make_config.py-compatible sample pool (config_support/)
│   ├── build_master_pool.py       # Build config_support/master_pool.csv from Fungi_BFD samples.csv + repr_assignments.tsv (targeted config builder)
│   ├── build_targeted_configs.py  # Render configs/<focal_short>_<batch>.csv per study in a batches/*.yaml spec (targeted config builder); --link-dir for a main.nf-ready per-batch symlink dir
│   ├── filter_candidates.py       # (standalone helper — not yet wired into pipeline)
│   ├── split_fasta.py             # (standalone helper)
│   └── summarize_clusters.py      # (standalone helper)
├── lib/                           # Shared Python package (lib/__init__.py present)
│   ├── hits.py                    # PARSERS dict + open_input + Hit namedtuple — used by parse_hits.py and parse_self_hits.py
│   ├── config_parser.py           # parse_config() → list[Sample]; short_to_group()
│   ├── fasta.py                   # FASTA utilities
│   ├── model_organisms.py         # ModelOrgAnnotator — YAML-driven gene name lookup
│   ├── clusters.py                # build_families() + FamilyIndex — mmseqs cluster -> gene-family grouping
│   ├── report_data.py             # build_payload() / build_core_payload() / build_losses_payload()
│   ├── report_template.py         # HTML_TEMPLATE — the novelties.html page (canvas heatmap, HTML/CSS/JS, no deps)
│   ├── skins.py                   # SKINS registry — every report page's CSS custom properties,
│   │                               #   plus the picker markup/JS and WCAG contrast helpers
│   ├── report_common.py           # Shared CSS/JS fragments for all report templates (page chrome,
│   │                               #   el()/linkout helpers, external-resource builder)
│   ├── index_page.py              # render_project_page() / render_gallery_page() — the
│   │                               #   view/<project>/report.html and view/index.html designs
│   ├── core_report_template.py    # CORE_HTML_TEMPLATE — the core.html page (single table, no deps)
│   ├── losses_report_template.py  # LOSSES_HTML_TEMPLATE — the losses.html page (single table, no deps)
│   ├── lineage.py                 # RANK_NAMES + lineage_match() — deepest shared lineage-rank name between two species
│   ├── master_pool.py             # MasterSample/load_master_pool(); make_short()/assign_shorts() (canonical Short-assignment rule, imported by bin/convert_bfd_samples.py); load_representative_picks() (repr_assignments.tsv join)
│   ├── trait_data.py              # load_trait_definitions()/load_traits() — config_support/traits/ YAML+CSV loader with undeclared-(trait,value) and none-coexistence hard errors
│   ├── targeted_selection.py      # select_nearest()/select_trait() — lineage-proximity-ranked ingroup-companion picking for bin/build_targeted_configs.py
│   └── Helpers.groovy             # Helpers.projectName(params) — derives results subdirectory name
├── tests/
│   ├── conftest.py
│   ├── test_config_parser.py
│   ├── test_hits.py
│   └── data/                      # test.csv + small FASTAs for -profile test
├── configs/
│   ├── samples.csv                 # Master species pool (Species,Strain,Protein,DNA,Short,Lineage) — source for bin/make_config.py
│   ├── pezizo4_asco.csv            # Main analysis config (Pezizomycotina ingroup, Asco outgroup)
│   ├── modelorgs.yaml             # Model organism YAML for gene name lookups
│   └── batches/                   # YAML batch specs for bin/build_targeted_configs.py — one focal
│       │                          #   species + companion-picker + shared outgroup pool per study
│       └── mucoromycota_focal_v1.yaml  # Worked example: Mucor/Phycomyces/Basidiobolus/Batrachochytrium vs a shared Dikarya reference pool
├── config_support/                # Intermediate sample pools, not run configs — e.g. bin/convert_bfd_samples.py
│                                   #   output, passed back into make_config.py's --samples. Never passed
│                                   #   directly as --config to main.nf. Also holds master_pool.csv
│                                   #   (bin/build_master_pool.py output) and traits/ (config_support/traits/,
│                                   #   see bin/build_targeted_configs.py).
└── results/
    └── <config_basename>/
        ├── search_cache/          # storeDir — pairwise + self-hit raw outputs, never re-run
        ├── self_hits/             # Per-species paralog_cutoffs.tsv files
        ├── clusters/              # mmseqs2 cluster FASTA and TSV
        ├── presence_matrix.tsv
        ├── presence_matrix.function.tsv
        ├── candidates.fa
        ├── candidates.txt
        ├── tblastn_summary.tsv
        ├── novelties.<SHORT>.tsv  # One per ingroup species
        ├── loss_presence_matrix.tsv
        ├── loss_presence_matrix.function.tsv
        ├── loss_candidates.fa
        ├── loss_candidates.txt
        ├── loss_tblastn_summary.tsv
        ├── novelties.html         # Self-contained interactive novelty-candidate report
        ├── core.html              # Self-contained interactive core-genes report
        └── losses.html            # Self-contained interactive candidate gene-loss report

view/                              # sibling of results/ — one shareable folder per project
└── <config_basename>/
    ├── report.html               # Landing page (run summary + links to the three reports)
    ├── novelties.html            # copies of the three results/ reports
    ├── core.html
    └── losses.html
```

## Architecture: Data Flow

1. **Config CSV** is parsed in `main.nf` using `Channel.fromPath(...).splitCsv(header:true)`. FASTA paths are resolved via `resolve_fa()`, which searches flat layout and subdirs (`pep/`, `dna/`, `genome/`, `scaffolds/`). Four channels are emitted: `ingroup_prot_ch`, `outgroup_prot_ch`, `outgroup_dna_ch`, `ingroup_dna_ch` (the last used only by the loss-search direction's TBLASTN validation).

2. **SEARCH workflow** (`workflows/search.nf`):
   - Cross-joins every ingroup proteome against all other proteomes (ingroup ∪ outgroup, excluding self-pairs).
   - Runs the selected tool (`PHMMER_SEARCH` / `DIAMOND_SEARCH` / `BLAST_SEARCH`), output cached in `search_cache/` via `storeDir`.
   - `PARSE_HITS` normalises each raw result into a TSV with columns `query_id, target_id, evalue, bitscore, query_proteome, target_proteome`.
   - Also runs self-vs-self searches (`PHMMER_SELF` / `DIAMOND_SELF` / `BLAST_SELF`, `-E 100 --max-target-seqs 2`) to capture the rank-2 (best paralog) hit per protein. Outputs also `storeDir`-cached.
   - `PARSE_SELF_HITS` → `<Short>.paralog_cutoffs.tsv` (columns: `protein_ID, paralog_protein_ID, bitscore, evalue`). Proteins with no within-proteome paralog are omitted. Only the `paralog_protein_ID` column is consumed downstream (by the paralog-competition filter); the `evalue` column is not used for gating (see the SEARCH workflow notes above).
   - `BUILD_PRESENCE_MATRIX` applies two filters:
     1. **Significance filter**: hit evalue must be < `--default-evalue` (flat, 1e-5), applied to every hit. (2026-09-03: this was previously a per-query "paralog-cutoff" — hit evalue < the query's own within-proteome paralog evalue — dropped as unsound: it swung both too tight, an unreachable bar whenever the query happened to have a close in-genome paralog regardless of relevance, and too loose, since the self-search that produces the paralog evalue is reported down to E=100 rather than a real significance floor, so a protein with no genuine paralog still got a noise-derived "cutoff" looser than the intended default. See `lib/singleton_presence.py`'s module docstring for the empirical writeup.)
     2. **Paralog-competition filter**: if the query's paralog hits the same target proteome with a better evalue, the hit is disqualified — this is the actual paralogy test (a direct head-to-head comparison), unaffected by the change above.
   - The `candidates.txt` keep filter is `query_frac >= min_frac AND other_frac <= --other-max-frac` (`other_frac` = fraction of the *other* group the protein is present in). `--other-max-frac` defaults to `0.0` — absent from every other-group proteome, the strict novelty/loss rule. The loss direction passes `params.loss_ingroup_max_frac` here to allow "nearly missing" candidates; the novelty direction always passes `0.0`. Note: this filter shapes only `candidates.txt`, not the matrix — the matrix always holds every scored row.
   - Produces `presence_matrix.tsv` (protein × proteome 0/1 matrix with `protein_id` and `source_proteome` columns) and `candidates.txt` (lines of `source_proteome::protein_id`).

3. **CLUSTER workflow** (`workflows/cluster.nf`):
   - `EXTRACT_CANDIDATES` reads `candidates.txt` and pulls matching sequences from all ingroup proteome FASTAs into `candidates.fa`.
   - `MMSEQS_CLUSTER` runs `mmseqs2 easy-cluster` (`--min-seq-id 0.3 -c 0.8 --cov-mode 0`) → `clusters_rep_seq.fasta`, `clusters_all_seqs.fasta`, `clusters_cluster.tsv`.
   - Emits: `candidates_fa`, `representatives`, `cluster_tsv`.

4. **VALIDATE workflow** (`workflows/validate.nf`):
   - `TBLASTN` searches all cluster representatives against each outgroup genome (one job per genome).
   - `SUMMARIZE_TBLASTN` merges per-genome TBLASTN TSVs, expands rep-level hits to cluster members, produces `tblastn_summary.tsv` (protein × outgroup-genome 0/1 matrix).

5. **ANNOTATE workflow** (`workflows/annotate.nf`):
   - Single `ANNOTATE_MATRIX` process runs (conditionally, based on whether paths are set):
     - `hmmscan` vs Pfam-A → `candidates.pfam.tblout` (MPI-capable via `--hmm_mpi`). Query is
       candidates.fa, target is the hmmpress'd Pfam-A database — this orientation (not
       hmmsearch) is what lets multithreading actually help: HMMER parallelizes over the
       target, and candidates.fa is far smaller than Pfam-A's ~20,800 profiles, so
       hmmsearch's threading barely engaged regardless of `--cpu` (observed ~2.8 real
       cores no matter how many were requested).
     - `diamond blastp` vs SwissProt → `candidates.swissprot.tsv` (best hit, outfmt 6 with stitle).
   - `annotate_presence_matrix.py` merges hits into the matrix with annotation priority:
     1. Model organism gene names (via `--modelorgs_config` YAML).
     2. Pfam domain names (all unique domains per protein).
     3. SwissProt description (best hit, `sp|ACCN|ID` prefix stripped).
   - Adds columns: `gene_name`, `product_description`, `function_source`, `Best_Swissprot`, `Pfam_Names`.
   - Produces `presence_matrix.function.tsv`.

6. **SUMMARIZE workflow** (`workflows/summarize.nf`):
   - `MAKE_NOVELTIES` calls `make_novelties.py` to produce one `novelties.<SHORT>.tsv` per ingroup species.
   - A protein is a novelty candidate if:
     - `source_proteome == short` (originates from that species).
     - Present in ≥ `ingroup_min_frac` of all ingroup species.
     - Absent from every outgroup proteome.
     - No TBLASTN hit in any outgroup genome (unless `--skip_tblastn_filter` is passed — hits still reported in `tblastn_outgroup_hits` column).
   - Currently wired with `--skip_tblastn_filter` in `workflows/summarize.nf`.

7. **Loss direction — candidate lineage-specific gene losses** (present in the outgroup,
   absent from the ingroup), run from `main.nf` right after the ingroup-direction stages.
   Like the novelty direction, `--cluster_tool` selects the producer:
   - **`pairwise`**: `LOSS_SEARCH → LOSS_CLUSTER → LOSS_VALIDATE → LOSS_ANNOTATE` (described
     below) — a mirror of steps 2–5 with ingroup/outgroup roles swapped.
   - **`mmseqs`**: the outgroup-seeded family-profile mirror (`PROFILE_LOSS_SEARCH`, an alias
     of `PROFILE_SEARCH` seeded from the outgroup with `--query-group OUT`, `out_prefix
     'loss_'`), then `LOSS_PROFILE_CANDIDATE_CLUSTERS` (family-as-cluster) feeding
     `LOSS_VALIDATE`/`LOSS_ANNOTATE`. Its family intermediates publish to `loss_families/`,
     `loss_family_hmmsearch/` so they don't collide with the novelty direction's. Emits the
     same `loss_presence_matrix.tsv` / `loss_candidates.txt` contract as the pairwise path.
   Both paths keep families/candidates present in ≥ `outgroup_min_frac` of the outgroup and
   ≤ `loss_ingroup_max_frac` of the ingroup, TBLASTN-validate against ingroup genomes, and
   render `losses.html`.

   The pairwise path in detail:
   - `workflows/loss_search.nf`'s `LOSS_SEARCH` runs outgroup proteomes as the *query*
     (`workflows/search.nf` only ever queries with the ingroup — see that file's module
     docstring for why this needs its own search direction rather than reusing SEARCH's
     output). Same tool selection, same self-search-based paralog lookup, same
     `BUILD_PRESENCE_MATRIX` process — called with `--query-group OUT` and
     `outgroup_min_frac` instead of `IN`/`ingroup_min_frac`. Produces
     `loss_presence_matrix.tsv` / `loss_candidates.txt`: a row exists only when that
     outgroup protein is present in ≥ `outgroup_min_frac` of the outgroup **and** absent
     from every ingroup proteome.
   - `LOSS_CLUSTER` is `workflows/cluster.nf`'s `CLUSTER` workflow imported under an
     alias (`include { CLUSTER as LOSS_CLUSTER } from './workflows/cluster'`) and called
     with `outgroup_prot_ch` (so `EXTRACT_CANDIDATES` pulls sequences from outgroup
     FASTAs) and distinct output names (`loss_candidates.fa`, `loss_clusters` mmseqs
     prefix) so it does not overwrite the ingroup-direction `CLUSTER` call's outputs in
     the same `publishDir`.
   - `LOSS_VALIDATE` is `VALIDATE` aliased likewise, called with `ingroup_dna_ch` as the
     TBLASTN target genomes (`TBLASTN_MAKEDB`/`TBLASTN` are already keyed by
     `meta_genome.id` in `storeDir`/output filenames, so no collision even without an
     alias there) → `loss_tblastn_summary.tsv`.
   - `LOSS_ANNOTATE` is `ANNOTATE` aliased likewise, called with an `output_prefix` of
     `'loss_'` → `loss_presence_matrix.function.tsv`.
   - Every reused process (`BUILD_PRESENCE_MATRIX`, `EXTRACT_CANDIDATES`,
     `MMSEQS_CLUSTER`, `SUMMARIZE_TBLASTN`, `ANNOTATE_MATRIX`) took a small,
     backward-compatible parameterization for this — an explicit query-group flag or
     output filename/prefix `val` input, always defaulted so the original
     ingroup-direction call sites in `main.nf` are unaffected. See
     `bin/build_presence_matrix.py`'s `--query-group` and `bin/extract_candidates.py`'s
     module docstring for the underlying logic.

8. **REPORT workflow** (`workflows/report.nf`):
   - `MAKE_REPORT` calls `make_report.py`, which merges the annotated matrix, the
     TBLASTN summary and the novelties tables into `novelties.html`.
   - `lib/report_data.py` builds the payload; `lib/report_template.py` holds the page.
   - The payload is embedded as JSON in a `<script type="application/json">` block, so
     the page opens from `file://` with no network access — reports get copied off the
     cluster. Any new asset must stay inlined.
   - Annotations come from the **matrix**; the novelties tables supply novelty status
     and protein sequences. When `--novelties` is omitted, novelty status is recomputed
     from the matrix (mirroring `make_novelties.py` with the TBLASTN filter skipped).
   - Sequences dominate file size, so `--sequences` (`params.report_sequences`) defaults
     to `novelties` — roughly 5 MB for a 21k-protein, 4.8k-novelty run.
   - `MAKE_CORE_REPORT` calls `make_core_report.py`, which asks the opposite question of
     the novelty report from the **same** annotated matrix — no new search or annotation
     step. A row (always ingroup-sourced, since that is the matrix's own scope) counts
     as core when its presence fraction across every proteome column (ingroup + outgroup)
     is `>= core_min_frac`. Sequences are never embedded in `core.html`; linkouts (Pfam →
     InterPro, SwissProt → UniProt/AlphaFold, model-organism → FungiDB, otherwise a
     generic NCBI Protein search) stand in for the BLASTP-with-sequence link the novelty
     report offers.
   - `MAKE_LOSSES_REPORT` calls `make_losses_report.py` against the **loss-direction**
     annotated matrix/TBLASTN summary/cluster_tsv (all from step 7, not step 2–6) —
     candidate lineage-specific gene losses, sourced from outgroup proteins conserved in
     the outgroup but (nearly) absent from the ingroup. The annotated matrix is the *full*
     presence table (`LOSS_ANNOTATE` runs on `LOSS_SEARCH.out.matrix`, not the filtered
     candidate list), so `build_losses_payload()` **re-applies the candidate predicate
     itself** — keep a row iff `outgroup_frac >= outgroup_min_frac AND ingroup_frac <=
     loss_ingroup_max_frac` — exactly mirroring `build_core_payload()`/`derive_novelties()`.
     Pass the same `--outgroup_min_frac`/`--loss_ingroup_max_frac` the loss search used so
     the reported rows match `loss_candidates.txt` and hence the mmseqs families built from
     it. Each row also carries two **family-level** aggregates — `out_breadth` (# distinct
     outgroup species carrying any family member) and `in_retained` (# ingroup species that
     still retain one; `0` = clean loss) — since the biological unit of a loss is the gene
     family, not one outgroup protein. Default priority sort: fewest `in_retained`, then
     widest `out_breadth`, then no ingroup-genome TBLASTN hit (reporting-only, not
     filtering — same `--skip_tblastn_filter` rationale as `MAKE_NOVELTIES`). No sequences
     embedded, same linkout fallback chain as `core.html`, resolved against the *outgroup*
     protein (that's where the gene is).
   - `COLLATE_REPORTS` (final `REPORT` step) copies `novelties.html`, `core.html` and
     `losses.html` into `view/<project>/` and generates `report.html` there via
     `bin/make_index_report.py` — a self-contained landing page linking the three reports
     with a run summary (ingroup/outgroup proteomes, search tool, thresholds). Published to
     `view/<project>/` (not `results/`) so the whole set is one shareable folder.
   - `lib/clusters.py`'s `FamilyIndex` (mmseqs cluster membership → per-protein family
     index + per-family species set) is shared by `build_payload()`,
     `build_core_payload()` and `build_losses_payload()`; `lib/report_common.py` holds
     the CSS/JS fragments shared by the single-table templates (`core_report_template.py`,
     `losses_report_template.py`) — `report_template.py`'s canvas-heatmap page is not
     wired to it.

9. **Final outputs** in `results/<project>/`: `presence_matrix.tsv`, `presence_matrix.function.tsv`, `candidates.fa`, `tblastn_summary.tsv`, `novelties.<SHORT>.tsv` for each ingroup species, the loss-direction equivalents (`loss_presence_matrix.tsv`, `loss_presence_matrix.function.tsv`, `loss_candidates.fa`, `loss_tblastn_summary.tsv`), `novelties.html`, `core.html`, and `losses.html`. A final
`COLLATE_REPORTS` step additionally writes `view/<project>/` containing copies of the three
reports and a `report.html` landing page (run summary + links).

## Key Parameters (`nextflow.config`)

| Parameter | Default | Notes |
|---|---|---|
| `--config` | (required) | Path to analysis CSV |
| `--data_dir` | (required) | Directory containing all FASTA files; also searched in `pep/`, `dna/`, `genome/`, `scaffolds/` subdirs |
| `--project` | CSV basename | Results subdirectory under `--outdir` |
| `--outdir` | `results` | Root output directory |
| `--cluster_tool` | `pairwise` | Presence-matrix producer: `pairwise` (O(N²) SEARCH) or `mmseqs` (family-profile pathway, ADR-0002). Applies to both novelty and loss directions |
| `--run_tool` | `phmmer` | `phmmer`, `diamond`, or `blast` (pairwise search + self-search paralog lookup) |
| `--family_min_seq_id` | `0.3` | mmseqs family-clustering identity (`--cluster_tool mmseqs`) |
| `--family_cov` | `0.8` | mmseqs family-clustering coverage (`--cluster_tool mmseqs`) |
| `--hmm_presence_evalue` | `1e-3` | Family-HMM presence E-value ceiling (`--cluster_tool mmseqs` and `novelty_discovery`); swept by `bin/run_param_sweep.sh` |
| `--hmm_presence_cov` | `0.5` | Family-HMM min profile coverage *fraction* for presence (`--cluster_tool mmseqs` and `novelty_discovery`); `bin/run_param_sweep.sh` now includes it as a swept dimension (2026-09-03) but no sweep has actually been run against it yet — the shipped 0.5 is still a guess. See `--hmm_presence_min_residues` below and the 2026-09-03 finding that this floor rejects 83% of otherwise-significant hits on HMMs >600aa (expected biology for multi-domain proteins, not promiscuous-domain false positives) |
| `--hmm_presence_min_residues` | `100` | Alternative to `--hmm_presence_cov`: a target qualifies if its aligned span is at least this many residues, even below the coverage fraction (`0` = disabled, fraction-only). Targets the coverage floor's actual purpose (reject a hit explained by one small, promiscuous shared domain) without penalizing long proteins purely for length; the `100` default is data-informed (real-run sensitivity analysis) but not yet swept/validated the way `--family_min_seq_id` etc. are. See `lib/family_presence.py` |
| `--family_chunk_size` | `200` | Families per parallel `BUILD_CHUNK` task (`--cluster_tool mmseqs`) |
| `--evalue` | `1e-5` | Flat significance cutoff applied to every hit in `BUILD_PRESENCE_MATRIX`/`novelty_discovery`'s singleton path (not per-query paralog-derived, see the SEARCH workflow notes above); also TBLASTN significance cutoff |
| `--parse_evalue` | `0.01` | Loose noise ceiling passed to `parse_hits.py`; final filtering applies the flat `--evalue` significance cutoff plus the paralog-competition filter |
| `--ingroup_min_frac` | `0.75` | Fraction of ingroup proteomes that must contain a hit |
| `--outgroup_min_frac` | `0.75` | Fraction of outgroup proteomes that must contain a hit, for `LOSS_SEARCH` (the loss-search direction's own presence threshold) |
| `--loss_ingroup_max_frac` | `0.0` | Max fraction of the ingroup a loss candidate may still be present in (loss-search direction). `0.0` = strictly absent from the ingroup; raise for "nearly missing" losses. Wired through to `build_presence_matrix.py --other-max-frac` and `make_losses_report.py --loss_ingroup_max_frac` |
| `--core_min_frac` | `0.95` | Presence fraction (across all proteomes, ingroup + outgroup) for the CORE genes report |
| `--use_orthofinder` | `false` | Placeholder — OrthoFinder clustering not yet implemented |
| `--pfam_hmm` | `null` | Path to Pfam-A.hmm; skips Pfam annotation if unset |
| `--swissprot_dmnd` | `null` | Path to SwissProt `.dmnd`; skips SwissProt annotation if unset |
| `--modelorgs_config` | `null` | Absolute path to model organisms YAML; skips gene-name lookup if unset |
| `--report_sequences` | `novelties` | Which proteins carry a sequence in `novelties.html`: `novelties`, `all`, or `none`. Sequences dominate the file size |
| `--pdf_report` | `true` | Write `view/<project>/summary.pdf` (matplotlib figures). Set `false` to skip the step (gated via the process `when:` directive) |
| `--hmm_mpi` | `false` | Run hmmsearch with MPI (`mpirun -np <hmm_mpi_tasks> hmmsearch --mpi`) |
| `--hmm_mpi_tasks` | `null` | MPI task count; defaults to `max_cpus` when `hmm_mpi=true` |
| `--max_cpus` | `32` | Cluster-wide CPU cap |
| `--max_memory` | `64.GB` | Cluster-wide memory cap |
| `--max_time` | `24.h` | Cluster-wide time cap |

## Process Labels and Resources

| Label | Default CPUs | Default Memory | Use for |
|---|---|---|---|
| `high_cpu` | 8 (32 on SLURM) | 16 GB (32 GB on SLURM) | Diamond, Blast, mmseqs2, ANNOTATE_MATRIX |
| `med_cpu` | 4 | 8 GB | Self-searches |
| `hmmsearch` | 4 | 8 GB | hmmsearch (MPI-aware on SLURM) |
| `low_cpu` | 1 | 4 GB | Parse steps, matrix building, make_novelties |

All processes activate the pixi environment via `beforeScript` and symlink `db/` and `configs/` into the work directory.

## Running the Pipeline

### Minimal run (local executor)

```bash
nextflow run main.nf \
    --config configs/pezizo4_asco.csv \
    --data_dir /path/to/fastas \
    --run_tool phmmer
```

### SLURM cluster

```bash
nextflow run main.nf \
    -profile slurm \
    --config configs/pezizo4_asco.csv \
    --data_dir /path/to/fastas \
    --run_tool diamond \
    --ingroup_min_frac 0.8 \
    --pfam_hmm db/pfam/Pfam-A.hmm \
    --swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd \
    --modelorgs_config configs/modelorgs.yaml
```

### Built-in test profile

```bash
nextflow run main.nf -profile test
```

Sets `params.config`, `params.data_dir`, and `params.project` from `tests/data/test.csv` and uses the local executor.

### Resuming after failure

```bash
nextflow run main.nf -resume --config configs/... --data_dir ...
```

`-resume` reuses completed tasks from `work/`. The `search_cache/` directory is additionally protected by `storeDir` — those steps are never re-run even across separate invocations.

### Running multiple pipelines concurrently

Nextflow keys a run's session (`.nextflow/history`, `-resume`'s "most recent run"
lookup) to the **launch directory**, not `--project`/`--outdir` — two runs started
from the same directory at the same time race for that state, and a `-resume` on
one can attach to the other's still-running session. Launch each concurrent run
from its own `.nf_launch/<project>/` subdirectory instead (already `.gitignore`d),
with every path-bearing argument converted to an absolute path first (relative
paths would otherwise resolve against the isolated subdir, not the repo root).
See the README's "Running multiple pipelines concurrently" for the full pattern
and a worked example, and `bin/run_param_sweep.sh` / the `run_*_refresh.sh`
scripts in the repo root for real usages. This only isolates Nextflow's own
bookkeeping — concurrent runs still share the SLURM queue, and since Nextflow
always stages `bin/`/`lib/` live from wherever `main.nf` was launched — the
Docker image ships only the external tool binaries, never the repo's own
Python code, so this applies to containerized runs too — do **not** `git
checkout`/`merge`/`rebase` in the repo root while any run is still using it.
See `.living/learnings.md`'s 2026-09-05 entry on the live-checkout race this
surfaced.

### Switching search tools

Pass `--run_tool phmmer|diamond|blast`. Results are keyed by tool name in the cache filenames (`*.phmmer.tblout.gz`, `*.diamond.tsv.gz`, `*.blast.tsv.gz`) so cached results from one tool are not reused when switching.

### Switching the presence-matrix producer (`--cluster_tool`)

`--cluster_tool pairwise` (default) runs the O(N²) `SEARCH`/`LOSS_SEARCH` workflows;
`--cluster_tool mmseqs` runs the scalable family-profile pathway (`PROFILE_SEARCH` /
`PROFILE_LOSS_SEARCH`, ADR-0002) — cluster the seed group into gene families, build family
HMMs (famsa + hmmbuild, parallelised across `BUILD_CHUNK` tasks), scan every proteome. Both
emit the same `presence_matrix.tsv`/`candidates.txt` (and `loss_` equivalents) contract, so
CLUSTER/VALIDATE/ANNOTATE/REPORT are unchanged. `--run_tool` is still used under `mmseqs`
for the self-search paralog lookup. See `README.md` "Two analysis approaches" for the
user-facing comparison and copy-paste launch examples.

`--cluster_tool novelty_discovery` runs a third, target-focused pathway (`NOVELTY_DISCOVERY`
→ `NOVELTY_SCREEN`) for a small known target lineage against a small reference panel, then a
broader clade-vs-distant-lineage screen — see "Two-Phase Targeted Novelty Pipeline" below.
It requires the config CSV to use `DISCOVERY_TARGET`/`DISCOVERY_OUT`/`NEAR_INGROUP`/`BROAD_OUTGROUP` GROUP values
instead of (or alongside) `IN`/`OUT`.

---

## Development Guide

### Implementation workflow (tickets → branch → PR)

When moving from planning to implementation, do **not** edit code directly on `main`.
Follow this order (the "to-tickets" flow):

1. **Turn the work into tickets first.** Before writing code, create a GitHub issue per
   discrete unit of work using the `gh` CLI (see `docs/agents/issue-tracker.md` for the
   conventions). Source the tickets from the plan artifacts already in the repo — the
   `todo/` registry items and the relevant `docs/adr/` decision — and link back to them in
   the issue body (`See docs/adr/000N-*.md`, `todo/<item>.md`). One issue = one PR-sized
   change. For a multi-phase plan, open one issue per phase (or smaller) rather than a
   single mega-issue.
2. **Branch per ticket.** Create a feature branch off `main` named for the issue, e.g.
   `git switch -c NNN-short-slug` where `NNN` is the issue number. Never commit
   implementation work to `main` directly.
3. **Implement on the branch**, committing in small, reviewable steps. Reference the issue
   in commit messages (`... (#NNN)`) so GitHub links them. Keep the branch scoped to that
   one ticket.
4. **Open a PR** with `gh pr create`, body linking the issue with a closing keyword
   (`Closes #NNN`) and pointing at the ADR/todo it implements. The PR is where review and
   `/code-review` happen before merge.
5. **Update tracking on completion.** When the PR merges, the issue auto-closes via the
   keyword; flip the matching `todo/TODO_REGISTRY.md` row to `complete` and, if the change
   resolved or altered a decision, note it in `.living/decisions.md`.

Only create issues, branches, push, or open PRs when the user has asked to start
implementing — planning and design stay in `todo/` + `docs/adr/` until then. Confirm before
the first outward-facing action (issue creation, push, PR) unless already told to proceed.

### Where to add new code

| What you are adding | Where it goes |
|---|---|
| A new external tool | New file `modules/<toolname>.nf` |
| A new analysis stage | New file `workflows/<stage>.nf` |
| A Python helper called from a Nextflow process | `bin/<script>.py` (must be executable) |
| Shared Python logic | `lib/<module>.py` |
| Unit tests | `tests/test_<module>.py` |
| Test fixtures | `tests/data/` |
| A new analysis config | `configs/<name>.csv` |
| New pixi-managed dependencies | `pixi.toml` under `[workspace]` |

### Adding a module

1. Create `modules/<toolname>.nf`. One `process` block per file.
2. Name the process `SCREAMING_SNAKE_CASE` matching the file basename.
3. Add a `label` (`high_cpu`, `med_cpu`, `hmmsearch`, or `low_cpu`) so `nextflow.config` resource rules apply.
4. Use `storeDir` for expensive pairwise/self-search outputs; use `publishDir` only for final user-facing results.
5. Write a minimal test in `tests/` covering the Python helper if one exists.

### Adding a workflow

1. Create `workflows/<stage>.nf`. Import only the modules it directly calls.
2. Declare `take:`, `main:`, and `emit:` blocks — never omit `emit:`.
3. Import and call the new workflow from `main.nf`.

### Modifying `bin/` scripts

- Run `chmod +x bin/<script>.py` after creation.
- All arguments must be `argparse` flags. Positional arguments break when Nextflow reorders tokens.
- Import shared logic from `lib/` via `sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))`.
- The common parsed-hit TSV columns are: `query_id, target_id, evalue, bitscore, query_proteome, target_proteome`.
- `candidates.txt` format: `source_proteome_short::protein_id` (one per line).

### Modifying `lib/`

- `lib/` is a proper Python package (`lib/__init__.py` exists).
- Current modules:
  - `hits.py` — `PARSERS` dict (keys: `phmmer`, `diamond`, `blast`), `open_input()` (handles `.gz`), `Hit` namedtuple.
  - `config_parser.py` — `parse_config()` → `list[Sample]`; `short_to_group()`.
  - `fasta.py` — FASTA utilities.
  - `model_organisms.py` — `ModelOrgAnnotator`; three ID-transform strategies: `direct`, `strip`, `diamond_fasta`.
  - `report_data.py` — `build_payload()` → the report's embedded JSON; `ROW_FIELDS` names the row array layout.
  - `report_template.py` — `HTML_TEMPLATE`; substitutes `__PROJECT_TITLE__` and `/*__PAYLOAD__*/`.
  - `lineage.py` — `RANK_NAMES`; `lineage_match()` (deepest shared rank NAME between two fixed-width lineages, missing-rank-safe).
  - `master_pool.py` — `MasterSample`/`load_master_pool()`; `make_short()` (the canonical Short-assignment rule — `bin/convert_bfd_samples.py` imports it rather than keeping its own copy) + `assign_shorts()`; `load_representative_picks()` (repr_assignments.tsv join, hard-errors on zero/multiple `is_representative=True` rows for a species the table actually covers).
  - `trait_data.py` — `load_trait_definitions()`/`load_traits()`; hard-errors on an undeclared `(trait, value)` pair or a `none` value coexisting with another value for the same species+trait.
  - `targeted_selection.py` — `select_nearest()`/`select_trait()`; `DEFAULT_SCOPE_RANK`.
- Keep `lib/` free of I/O side effects at import time.

---

## The interactive report (`novelties.html`)

Regenerate standalone from an existing results directory without re-running the pipeline:

```bash
bin/make_report.py \
    --matrix results/pezizo4_asco/presence_matrix.function.tsv \
    --config configs/pezizo4_asco.csv \
    --tblastn_summary results/pezizo4_asco/tblastn_summary.tsv \
    --novelties results/pezizo4_asco/novelties.*.tsv \
    --candidates_fa results/pezizo4_asco/candidates.fa \
    --output results/pezizo4_asco/novelties.html
```

### Constraints to preserve when editing the page

- **Self-contained.** No CDN scripts, external stylesheets, fonts, or fetches — the page
  must open from `file://` offline. Inline anything new.
- **Untrusted strings.** Protein IDs and product descriptions come from FASTA headers and
  SwissProt/Pfam text. Insert them with `textContent` (the `el()` helper), never
  `innerHTML`. `make_report.py` escapes `</` in the payload so an annotation cannot close
  the `<script>` block early — `tests/test_report_data.py` covers this.
- **Colour lives in `lib/skins.py`, never in a template.** All four pages emit
  `skins_css()`; a raw hex in a template is a test failure
  (`tests/test_report_templates.py::test_report_has_no_hardcoded_theme_colours`). Adding
  a palette means adding a `SKINS` entry — it must define every token in
  `REQUIRED_TOKENS`, and `tests/test_skins.py` enforces the WCAG floors
  (4.5:1 body text, 3:1 data marks and `--on-series` chip labels).
- **Colour encodes evidence type, not group.** series-1 = presence from the protein
  search; series-2 = TBLASTN genome hit. Ingroup vs outgroup is carried by column
  position and the band labels, so no hue does double duty. The pair must stay separable
  under red-green CVD — which is *hue* separation on the blue-yellow axis, something no
  contrast test can check, so reason about it explicitly when adding a skin (the neon
  skin pairs cyan with amber rather than the obvious cyan/magenta for exactly this
  reason). Re-validate before editing the shipped `paper`/`dark` values.
- **The canvas repaints on a skin change** via `window.onSkinChange` — it reads the same
  tokens through `getComputedStyle` (`palette()`), so a new skin needs no canvas code.
- **A skin owns its type and effects too**, through `--font-ui`/`--font-mono`/`--glow`/
  `--overlay`/`--shadow`. `--overlay` (the neon skin's scanlines) is a **background layer
  on `body`**, not an overlay element: an element's background always paints below all of
  its descendants, so the texture can only ever appear in the page's negative space. Do
  not reimplement it as a positioned overlay — the first attempt did, and a `z-index: 50`
  overlay painted scanlines across the detail panel (`position: sticky`, `z-index: auto`)
  and the heatmap tooltip (`z-index: 20`), because raising only the two scroll containers
  left every other content region below it. Webfonts are impossible here — the pages open
  from `file://` — so every font stack is system fonts.
- **Rows are virtualised** on a canvas (~20k rows). `drawGrid()` renders only the visible
  window; the scroll height comes from the spacer div.
- **The table tab is the accessibility twin** of the heatmap and must keep every value
  reachable without hovering.
- **External links come from one builder**, `externalLinksNode()` in
  `lib/report_common.py`. All three reports call it, so a protein resolves to the same
  records everywhere — the previous per-template copies had drifted. Per-species gene
  records come from the config CSV's optional `SourceDB` column, not from the annotation
  source. A candidate with neither a Pfam domain nor a SwissProt hit gets the
  remote-homology cluster (HHpred/Foldseek/InterProScan) instead of an ID-based NCBI
  search that would return nothing — that population is the point of the pipeline.
- **A `SourceDB` value is only semi-trusted.** Config CSVs get copied between users and
  projects, and the free-form `{gene}` template goes straight into an `href`, so
  `genomeDbLink()` requires an `http(s)` scheme before building a link — otherwise a
  `javascript:` value becomes a clickable link in the report.
  `tests/test_report_js_behaviour.py` covers this with a hostile fixture.

### Testing template JS changes

Stage 1 below now runs as part of the suite —
`tests/test_report_templates.py::test_page_javascript_parses` does the `node --check`
pass over every page's `<script>` bodies and skips where node is unavailable. That matters
more than it used to: the pages are assembled from shared fragments in
`lib/report_common.py` and `lib/skins.py`, so one bad fragment breaks all four at once.
Stage 2 (jsdom) now has a home too: `tests/test_report_js_behaviour.py` generates real
report pages from a fixture and drives them through `tests/js/drive_reports.mjs`,
covering the skin picker's three states and the whole external-link fallback chain. It
**skips unless jsdom resolves**, since that is an npm dependency in a pixi project — set
`NOVINVENIO_JSDOM` to a jsdom install, or `npm install jsdom` somewhere Node will find
it. Extend that driver rather than writing a new throwaway script.

When editing the JS inside `HTML_TEMPLATE` — or any shared fragment — verify it in two
stages:

1. **Syntax check with `node --check`** (or just run `pytest tests/test_report_templates.py`).
   To do it by hand, extract the `<script>` body (the second one — the first is the
   `application/json` payload) and check it parses:

   ```bash
   python3 -c "
   from lib.report_template import HTML_TEMPLATE
   import re
   m = re.search(r'<script>\n(.*)</script>\n</body>', HTML_TEMPLATE, re.S)
   open('/tmp/report_js_check.js', 'w').write(m.group(1))
   "
   node --check /tmp/report_js_check.js
   ```

2. **Drive a real generated report with jsdom** to exercise DOM wiring end-to-end (filters,
   selects, detail-panel buttons, table rendering). This project has no JS dependencies
   checked in, so install jsdom into a scratch dir rather than the repo:

   ```bash
   npm install --prefix /tmp/novinv_check jsdom --no-audit --no-fund --silent
   ```

   Generate a small `novelties.html` from a hand-written fixture matrix/config (same shape as
   the fixtures in `tests/test_report_data.py`) via `bin/make_report.py`, then load it with
   `JSDOM({ runScripts: 'dangerously', resources: 'usable', pretendToBeVisual: true })`.
   jsdom does not implement `<canvas>` or `matchMedia`, and the report page calls both
   unconditionally on load — stub them in the `beforeParse(window)` hook *before* the
   page's inline script runs, or every event handler that calls `drawGrid()` throws:

   ```js
   beforeParse(window) {
     window.matchMedia = () => ({ matches: false, addEventListener() {}, removeEventListener() {} });
     const ctxStub = {
       setTransform() {}, clearRect() {}, fillRect() {}, strokeRect() {},
       measureText: (t) => ({ width: (t || '').length * 6 }),
       fillText() {}, save() {}, restore() {}, translate() {}, rotate() {},
       beginPath() {}, moveTo() {}, lineTo() {}, stroke() {},
     };
     window.HTMLCanvasElement.prototype.getContext = () => ctxStub;
   }
   ```

   From there, drive the page like a user: set `select.value` / `input.value`, dispatch a
   plain `new window.Event('change'|'input'|'click', { bubbles: true })` (jsdom doesn't need
   real `MouseEvent`/`InputEvent` for this page's listeners), await a short `setTimeout` for
   the debounced search box (140ms), and assert on `textContent` of `#count`, `#tbl-body`,
   and `#detail` rather than the canvas heatmap (canvas is stubbed to a no-op, so pixel
   output isn't observable — assert through the table tab instead).

   Two things to stub beyond the list above when driving the *skin* wiring: the page
   writes to `localStorage`, so give jsdom a real `url` (a `file://` document has an
   opaque origin and storage throws), and remember that `--report_sequences` defaults to
   `novelties` — a non-novelty row carries no sequence, so it legitimately shows none of
   the sequence-driven links.

   This caught a real behavioural requirement during the gene-family feature: with default
   filters, a family member that doesn't independently clear its own species' novelty
   threshold stays hidden until the detail panel's "Show family members" button is clicked,
   which must override (not stack with) the novelty-only filter. A payload/unit test
   wouldn't have exercised that filter-override interaction — only driving the actual DOM did.

   Clean up the scratch install and generated fixtures afterward (`rm -rf /tmp/novinv_check`);
   none of it belongs in the repo.

---

## Nextflow Linking Rules

### DSL2 import paths

```nextflow
// From main.nf — workflows are siblings of main.nf
include { SEARCH  } from './workflows/search'

// From a workflow — modules are one level up
include { PHMMER_SEARCH } from '../modules/phmmer'
```

Relative paths are resolved from the file containing the `include`. Never use absolute paths.

### Channel contracts

Every channel element is a tuple `[meta, ...files]`. The `meta` map always contains at least:

```groovy
[id: String, group: 'IN'|'OUT', species: String, strain: String, taxon: String]
```

When a process produces a pair result, the meta is `[id: "${meta_q.id}_vs_${meta_t.id}"]`. Downstream processes that only need the file should `.map { meta, file -> file }` before collecting.

### Emit and intake conventions

- Every workflow `emit:` block names its outputs. Callers reference them as `WORKFLOW_NAME.out.<name>`, never positionally.
- Processes use named outputs (`emit: tsv`, `emit: matrix`, etc.).
- When collecting across samples: `.map { meta, tsv -> tsv }.collect()`.

### storeDir vs publishDir

| Directive | Use for | Behaviour |
|---|---|---|
| `storeDir` | Pairwise + self-search results in `search_cache/` | Skips the process if the output file already exists, even across pipeline runs |
| `publishDir` | Final results in `results/<project>/` | Copies/links on every run; does not affect task execution |

Never apply both to the same process output.

### Parameter access in processes

Access global parameters via `params.*` in the `script:` block. Only pass per-sample-varying values through channel tuples.

---

## Nextflow Conventions

- **DSL2** everywhere (`nextflow.enable.dsl=2` at the top of every `.nf` file).
- Meta map: `[id: row.Short, group: row.GROUP, species: row.Species, strain: row.Strain, taxon: row.TaxonGroup]`.
- Cache filename convention: `${meta_query.id}_vs_${meta_target.id}.${params.run_tool}.<ext>.gz`.
- Self-search cache: `${meta.id}_vs_${meta.id}.${params.run_tool}.<ext>.gz`.

## Python Script Conventions

- All scripts in `bin/` must be executable (`chmod +x`) and have `#!/usr/bin/env python3`.
- Shared logic lives in `lib/`; import via `sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))`.
- Use `--input`/`--output` argparse flags; never hardcoded paths.

## Analysis Config CSV Format

```csv
GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup
OUT,Neolecta irregularis,,Neolecta_irregularis.proteins.fa,Neolecta_irregularis.scaffolds.fa,Neolecta_irregularis.gff3,Nirr,Taphrinomycotina
OUT,Schizosaccharomyces pombe,,Schizosaccharomyces_pombe.proteins.fa,Schizosaccharomyces_pombe.scaffolds.fa,Schizosaccharomyces_pombe.gff3,Spom,Taphrinomycotina
OUT,Saccharomyces cerevisiae,S288C,Saccharomyces_cerevisiae.aa.fa,Saccharomyces_cerevisiae.dna.fa,,Scer,Saccharomycotina
IN,Neurospora crassa,OR74A,Ncrassa.pep.fa,Ncrassa.dna.fa,,Ncra,Pezizomycotina
IN,Aspergillus fumigatus,Af239,Afum.pep.fa,Afum.dna.fa,,Afum,Pezizomycotina
IN,Pyronema omphalodes,CBS 144459,Pomp.pep.fa,Pomp.dna.fa,,Pomp,Pezizomycotina
IN,Coccidioidies immitis,WA_211,Cocci_WA211.pep.fa,Cocci_WA211.dna.fa,,Cimm,Pezizomycotina
```

- `Short` is a ≤8-char unique ID used in filenames and output tables throughout.
- `SourceDB` and `NCBI_TaxID` are optional, report-only, and never affect a
  presence/novelty call. `SourceDB` builds a per-gene database link in the reports'
  detail panels: `fungidb`, `mycocosm:<portal>`, `ensemblfungi:<species>`,
  `veupathdb:<project>`, or any URL template containing `{gene}`. `NCBI_TaxID` turns the
  species' NCBI Taxonomy link from a by-name search into a direct taxid lookup. Omitting
  either (or the whole column) is never an error — see `lib/report_common.py`'s
  `genomeDbLink()`/`taxonomyLink()`.
- `Strain` may be empty.
- The config CSV filename (without `.csv`) becomes the results output subdirectory name.
- `Protein` and `DNA` are basenames resolved relative to `--data_dir` (also checked under `pep/`, `dna/`, `genome/`, `scaffolds/` subdirs).
- `GFF3` is optional (may be an empty cell, or the column may be omitted entirely from
  older config CSVs — `lib/config_parser.py`'s `parse_config()` defaults it to `''`).
  When present, it's a basename resolved relative to `--data_dir` the same way as
  `Protein`/`DNA` (also checked under `gff3/`, in addition to `pep/`, `dna/`, `genome/`,
  `scaffolds/` — see `lib/gff3_genes.py`'s `GFF3_SEARCH_SUBDIRS`). It supplies the
  Chromosome/Start columns in `novelties.html`/`core.html`/`losses.html`
  (`lib/gff3_genes.py`, wired through `lib/report_data.py`'s payload builders). A missing
  or unresolvable `GFF3` value is never an error — that species' report rows simply carry
  no chrom/start.

### GFF3 chrom/start is per-protein-record, not per-gene

The chrom/start columns are resolved per protein/transcript ID (`lib/gff3_genes.py`'s
`lookup_gene_position()`), and the pipeline does not currently deduplicate or filter
alternative splice isoforms. A single gene with multiple annotated transcripts shows up
as multiple report rows — one per protein/transcript ID — each pointing to the same or
very similar chrom/start. A dedup/filtering pass may be worth adding later; see
`todo/TODO_REGISTRY.md`.

## Two-Phase Targeted Novelty Pipeline (`--cluster_tool novelty_discovery`)

A third `--cluster_tool` option, alongside `pairwise` and `mmseqs`: a target-focused
alternative to the whole-ingroup/whole-outgroup design of those two, for when you already
know the small target lineage you care about and want a tighter, faster screen than an
all-vs-all run across a large sample set. Full design rationale and rejected alternatives
are in `todo/novelty-discovery-screen.md`; this section documents what the implementation
actually does today.

### Config format — four new `GROUP` values

The same single analysis CSV (`--config`) carries all four roles together, alongside (or
instead of) `IN`/`OUT`:

| `GROUP` value | Role | Typical size |
|---|---|---|
| `DISCOVERY_TARGET` | The genome(s) to find novelties in | 2-3 |
| `DISCOVERY_OUT` | Small reference panel for the phase-1 absence call | 5-6 |
| `NEAR_INGROUP` | Close relatives of the target, same clade | 20-50 |
| `BROAD_OUTGROUP` | Distant lineages, outside the target clade | 20-100 |

See `configs/novelty_discovery_example.csv` for a minimal example. `lib/config_parser.py`'s
`GROUPS`/`INGROUP_ROLES`/`OUTGROUP_ROLES` constants are the source of truth for which values
are recognized and how they band into "ingroup-like" (`IN`, `DISCOVERY_TARGET`) vs "outgroup-like"
(`OUT`, `DISCOVERY_OUT`, `NEAR_INGROUP`, `BROAD_OUTGROUP`) for the report payload builders.

### Phase 1 — `workflows/novelty_discovery.nf` (`NOVELTY_DISCOVERY`)

1. Cluster `DISCOVERY_TARGET` proteomes into gene families with mmseqs2 (`MMSEQS_FAMILY_CLUSTER`).
2. Build a family HMM per multi-member family (famsa + hmmbuild scatter-gather,
   `BUILD_FAMILY_PROFILES` — the same module `--cluster_tool mmseqs` uses).
3. `hmmsearch` every family HMM against `DISCOVERY_TARGET` + `DISCOVERY_OUT` proteomes
   (`FAMILY_HMMSEARCH`).
4. `CALIBRATE_FAMILY_HMMS` still runs (Nextflow process interface stability) but its output
   is **not** used to gate presence (2026-09-03 fix): it used to set a family's presence
   threshold to the best E-value that family's own HMM scored against `DISCOVERY_OUT`, which
   is circular — no `DISCOVERY_OUT` proteome could ever satisfy `E < threshold` against
   itself, making the outgroup-absence filter a structural no-op (confirmed against a real
   run: three genes with highly significant, annotated `DISCOVERY_OUT` hits down to 3e-160
   were all called outgroup-absent). See `lib/family_presence.py`'s module docstring.
5. Build the phase-1 presence matrix and candidate list (`NOVELTY_PRESENCE_MATRIX` /
   `bin/novelty_presence_matrix.py`): presence uses the flat `--hmm_presence_evalue` /
   `--hmm_presence_cov` for every proteome (families) and the flat `--evalue` plus the
   paralog-competition filter (singletons) — a family/protein is a candidate when present in
   `>= --ingroup_min_frac` of `DISCOVERY_TARGET` and absent from every `DISCOVERY_OUT` proteome.
6. TBLASTN family representatives against `DISCOVERY_OUT` genomes for genomic validation
   (reporting-only, same rationale as `make_novelties.py --skip_tblastn_filter`).

Singleton target proteins (no multi-member family) are extracted (`EXTRACT_SINGLETONS`) and
searched with the "hybrid" pairwise-vs-DISCOVERY_OUT branch (issue #52 — see
`workflows/novelty_discovery.nf`'s `SINGLETON_*_SEARCH` processes and
`lib/singleton_presence.py`), so single-copy target-specific genes are not invisible to this
pathway.

**Known gap:** the design doc's single-`DISCOVERY_TARGET`-genome pairwise-search branch (skip clustering
entirely when `|DISCOVERY_TARGET| == 1`) is likewise not implemented — the pipeline always runs the
mmseqs family-clustering path regardless of `DISCOVERY_TARGET` count, so a lone target genome mostly
yields single-member "families" that phase 1 currently can't evaluate for the same reason.
Family HMM `storeDir` caching (keyed by a hash of the family's member IDs, so overlapping
target sets across runs reuse HMMs) is also not yet implemented — every run rebuilds every
family HMM from scratch.

### Phase 2 — `workflows/novelty_screen.nf` (`NOVELTY_SCREEN`)

Always runs after `NOVELTY_DISCOVERY` for this `--cluster_tool` (there is no flag to skip
it). Re-searches the *same* family HMMs against `NEAR_INGROUP` and `BROAD_OUTGROUP`
proteomes and reclassifies every phase-1 candidate (`NOVELTY_SCREEN_CLASSIFY` /
`bin/novelty_screen.py`) into one of three categories, written as a `novelty_category`
column on the screened presence matrix:

- **`target_specific`** — no hit in `NEAR_INGROUP` or `BROAD_OUTGROUP`.
- **`clade_specific`** — hit in `NEAR_INGROUP`, no hit in `BROAD_OUTGROUP`.
- **`false_novelty`** — hit in `BROAD_OUTGROUP` (removed from the screened candidate list, but
  kept — labelled — in the matrix for visibility).

TBLASTN also validates family representatives against `BROAD_OUTGROUP` genomes
(reporting-only), published as `screen_tblastn_summary.tsv`; the interactive report
currently only surfaces `NOVELTY_DISCOVERY`'s `DISCOVERY_OUT`-genome TBLASTN evidence
(`tblastn_summary.tsv`), not this one.

A config with no `NEAR_INGROUP`/`BROAD_OUTGROUP` rows degrades gracefully rather than erroring: zero
proteomes to search means zero hits, so every phase-1 candidate defaults to
`target_specific` — there's no broader-screen evidence to demote it with.

Only the screened candidates (`false_novelty` removed) proceed to `ANNOTATE` (Pfam +
SwissProt) — annotation is expensive, so there's no reason to spend it on families the
screen already ruled out. `novelty_category` survives `ANNOTATE_MATRIX` untouched
(`annotate_presence_matrix.py` only appends columns) and reaches the interactive report.

### Report integration

`novelties.html` renders `novelty_category` as a table column, a detail-panel field, and a
`<select>` filter (hidden unless the payload actually has category data — pairwise/mmseqs
runs never populate it). `core.html` and `losses.html` are unaffected: the loss direction
(present in the outgroup, absent from the ingroup) is an explicitly deferred future
extension for this pathway, not mirrored the way it is for `pairwise`/`mmseqs` — `main.nf`
stubs a valid zero-row `losses.html` (`modules/empty_loss_stub.nf`) purely so
`COLLATE_REPORTS` can still assemble `view/<project>/report.html`.

### Running it

```bash
nextflow run main.nf \
    --config configs/novelty_discovery_example.csv \
    --data_dir /path/to/fastas \
    --cluster_tool novelty_discovery \
    --run_tool phmmer
```

Same `--pfam_hmm`/`--swissprot_dmnd`/`--modelorgs_config`/report flags as the other two
`--cluster_tool` paths apply. `--ingroup_min_frac`/`--hmm_presence_evalue`/
`--hmm_presence_cov` are shared with `--cluster_tool mmseqs` and mean the same thing here.

## Targeted Config Builder (`bin/build_master_pool.py` + `bin/build_targeted_configs.py`)

Generates `configs/<name>.csv` run configs for the common "many small, targeted
comparative-genomics runs" usage pattern (a focal species + a few close relatives vs. a
shared, reusable outgroup pool) directly from BFD data, instead of hand-editing a config
CSV. Full design in `docs/superpowers/specs/2026-09-05-config-builder-design.md`.

1. **`bin/build_master_pool.py`** builds `config_support/master_pool.csv`
   (`Species,Strain,ProteinPath,DNAPath,Lineage,NCBI_TaxID`, absolute FASTA paths, one row
   per species) by joining `Fungi_BFD/samples.csv` against
   `_reuse_assignments/repr_assignments.tsv`'s per-species representative-assembly pick.
   A species with no row at all in `repr_assignments.tsv` is normal (that table only covers
   the annotated/ANI-assessed subset) and is skipped with a stderr count, not an error; a
   species the table *does* cover but with zero or more than one `is_representative=True`
   row is a hard error (a real data-integrity problem); a representative whose annotation
   FASTAs aren't found on disk is skipped with a stderr warning (mirrors
   `bin/convert_bfd_samples.py`'s existing missing-annotation-dir skip pattern); a
   representative dirname absent from `samples.csv` entirely is a hard error (the two files
   disagree). `lib/master_pool.py` holds `MasterSample`/`load_master_pool()`,
   `make_short()`/`assign_shorts()` (Short assigned at render time from the FULL pool, not
   stored — see the design doc's "Short is not stored in the master pool"), and
   `load_representative_picks()`.
2. **`bin/build_targeted_configs.py`** renders one `configs/<focal_short>_<batch>.csv` +
   `.map.tsv` (Species→Short trace) per `studies[]` entry in a `configs/batches/*.yaml`
   batch spec (see `configs/batches/mucoromycota_focal_v1.yaml` for a real, runnable
   example). Each study picks a focal species' ingroup companions via `mode: nearest`
   (lineage-proximity-ranked, `lib/targeted_selection.py::select_nearest`), `mode: trait`
   (same ranking, filtered to a `(trait, value)` pair from `config_support/traits/` —
   `select_trait`), or `mode: explicit` (a literal, hand-picked list — the escape hatch when
   neither automated mode finds a sensible candidate, e.g. a focal species that is the only
   member of its whole phylum with `repr_assignments.tsv` coverage). A named
   `outgroup_pools` entry is reused across studies and is always excluded from that study's
   ingroup candidate pool. Pass `--link-dir <dir>` to also symlink every selected species'
   Protein/DNA files into `<dir>/pep/` and `<dir>/dna/` and write **basenames** (not the
   master pool's absolute paths) into the rendered config — this is what makes the config
   directly usable as `--data_dir <dir>` with `main.nf` (`main.nf`'s `resolve_fa()` cannot
   resolve an absolute path). Omitting `--link-dir` writes the master pool's absolute paths
   straight through, which is not `main.nf`-ready.

```bash
bin/build_master_pool.py \
    --bfd-samples /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv \
    --annotation-dir /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/genome_annotation \
    --repr-assignments /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation/_reuse_assignments/repr_assignments.tsv \
    --output config_support/master_pool.csv

bin/build_targeted_configs.py \
    --master-pool config_support/master_pool.csv \
    --trait-definitions config_support/traits/trait_definitions.yaml \
    --traits config_support/traits/traits.csv \
    --batch-spec configs/batches/mucoromycota_focal_v1.yaml \
    --outdir configs/ \
    --link-dir data/mucoromycota_focal_v1
```

## Model Organisms YAML Format (`configs/modelorgs.yaml`)

Three `id_transform` strategies:

- `direct` — protein_id used as lookup key unchanged.
- `strip` — regex stripped from protein_id (e.g. `-T-p\d+$` → bare gene ID).
- `diamond_fasta` — protein_id → diamond hit → FungiDB protein_id → FASTA header gene field.

```yaml
model_organisms:
  - short: Ncra
    gene_names_csv: db/modelorgs/Neurospora_crassa_gene_names_FungiDB.csv
    id_transform: diamond_fasta
    diamond_hits: db/modelorgs/Ncra_vs_FungiDB_Ncra.diamond.tsv
    protein_fasta: db/modelorgs/FungiDB-68_NcrassaOR74A_AnnotatedProteins.fasta
    fasta_gene_field: gene
  - short: Afum
    gene_names_csv: db/modelorgs/Aspergillus_fumigatus_gene_names_FungiDB.csv
    id_transform: strip
    strip_pattern: "-T-p\\d+$"
```

## Agent skills

### Issue tracker

Issues live in the repo's GitHub Issues via the `gh` CLI. See `docs/agents/issue-tracker.md`.

### Triage labels

Five canonical triage roles mapped to same-named GitHub labels. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context — one `CONTEXT.md` and `docs/adr/` at the repo root. See `docs/agents/domain.md`.


## Living Repository (mycelium)

This repo carries a `.living/` knowledge layer. Treat it as durable project memory:

- `.living/decisions.md` — append non-obvious choices (context, alternatives, rationale).
- `.living/learnings.md` — append gotchas/surprises with `**Tags**:` (feeds `.living/INDEX.md`).
- `.living/conventions.md` — repo-specific overrides to convention-pack defaults.
- `.living/INDEX.md` — auto-generated quick-reference + tag summary (regenerated by `generate_index.py`).
- `.living/log/LOG_REGISTRY.md` — per-session work log.
- `todo/TODO_REGISTRY.md` — future-work items (add via `/mycelium:core todo-idea`).

**Post-action protocol:** after any significant analysis/pipeline change, update the
relevant `.living/` files (decisions/learnings), then regenerate the index. Only the
`.living/` knowledge layer, `todo/`, and `skillpacks/` were scaffolded — the generic
`analysis/ algorithms/ reference_material/` data-analysis dirs were intentionally skipped
since this repo has its own Nextflow layout. SessionStart/Stop/PostToolUse hooks live in
`.claude/settings.json`.

## Installed Convention Packs

- **bioinformatics** — See `.living/conventions/bioinformatics/analysis-conventions.md`

- **idea-generator** — See `.living/conventions/idea-generator/analysis-conventions.md`

- **report-generator** — See `.living/conventions/report-generator/analysis-conventions.md`

- **robust-analysis** — See `.living/conventions/robust-analysis/analysis-conventions.md`
