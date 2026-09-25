# Pangenome gain/loss visualization: island synteny view + trans-module genomic clustering

**Status**: design approved 2026-09-19; spec-reviewed and revised 2026-09-19 (see
"Resolved during spec review"); not yet implemented. Build in two phases --
see "Implementation phasing".
**Scope**: two new visualizations for the pangenome subworkflow, sharing one
report page and one annotation lookup. Does not change any existing analysis
output.

## Motivation

The pangenome subworkflow computes two complementary kinds of gene-family
gain/loss relationship, and neither is currently visualizable:

1. **Accessory islands** — physically contiguous, presence/absence-variable
   loci (`BUILD_ISLANDS` → `islands_with_domains.tsv`). Gene *order* is
   meaningful here, so the natural question is "where does the block break?"
   — which strains lost which stretch of the locus.
2. **Trans-correlated families** — families whose presence/absence covaries
   significantly but which are confirmed *not* physically adjacent
   (`PAIR_CLASSIFICATION`'s `trans` class), grouped into modules by
   `LEIDEN_MODULES`. Gene order between them is by definition not the point;
   the question is whether they nonetheless cluster in genomic neighborhoods.

Existing pangenome figures are all genome-wide aggregates (frequency
histogram, core/shell/cloud pie, one global presence/absence matrix). None
shows locus-level structure or per-module genomic context.

## View A — Island synteny / presence-absence view

### What it shows

For one selected island: a presence/absence grid, rows = strains, columns =
member families **in locus order**, filled cell = family present in that
strain. Deletion breakpoints appear as vertical edges.

### Design decisions

- **Row collapsing is the load-bearing decision.** 530 strains × a 20-56 gene
  island is unreadable as 530 rows, but the number of *distinct presence
  patterns* (haplotypes) is typically tens. Collapse identical rows into one
  row carrying a `×N` count badge. Sort rows lexicographically on the 0/1
  presence string (equivalent to Hamming-distance leaf ordering), which
  groups related haplotypes and makes breakpoints align visually.
- **Alternate row sorts** offered: by species (immitis/posadasii band), by
  presence count, by strain name. Phylogeny-ordered is explicitly *not*
  offered — this pipeline has no strain tree (see
  `todo/pangenome-phylogeny-aware-gain-loss.md`); offering it would imply a
  capability that does not exist.
- **Domain annotation goes in a glyph strip above the grid**, not in cells:
  one colored tick per column keyed to a Pfam class (NLR/incompatibility =
  NACHT/Ank/HET; secondary metabolite = PKS/NRPS-associated; transporter;
  other; unannotated). Hover reveals the full domain string and family ID.
  This keeps annotation adjacent to the data without text clutter.
- **Island navigation**: a sidebar lists the selected islands with size,
  distinct-haplotype count, dominant Pfam class chip, and carrying-strain
  count; sortable/filterable. Selecting one redraws a single canvas.

### Selection

Top *N* islands by size (`--pangenome_viz_top_islands`, default 50), **after
excluding islands carried by fewer than `--pangenome_top_islands_min_strains`
strains** (default 2, the same parameter and default `report.md`'s top-islands
table uses).

That filter is not optional here. Measured on `genus_vs_ureesii`:
**17,218 of 27,836 located islands (62%) are present in exactly one strain**,
and all 20 of the largest were single-strain. An unfiltered size ranking
therefore fills the sidebar with strain-private content -- often an
assembly/annotation artefact, e.g. a 51 kb, 69-gene-model "island" on one
contig of one strain. A single-strain island is also *meaningless in this
view*: its presence/absence grid is one filled row and N-1 empty ones, with
no breakpoint to see. This is the same defect PR #110 fixed in `report.md`;
the viz must not reintroduce it.

Rendering all ~27,800 islands from a genus-scale run is not a goal.

### Explicitly out of scope

Distinguishing "present but rearranged" from "present in expected order" —
the current data model records per-strain gene positions but no synteny-block
assignment, so a rearrangement encoding would be inferred, not measured.

## View B — Trans-module genomic clustering

A hybrid, because a purely reference-anchored view was measured to be
infeasible (see Measured constraints).

### B1 — Primary: per-strain neighborhood statistics (reference-free)

For each module, for each strain carrying ≥2 of its member families, measure
how close those members sit **within that strain's own assembly**, and
compare against a permutation null drawn from that same strain's gene
complement. Each strain's internal coordinates are fully valid; positions are
never compared *across* strains, so no shared coordinate system is needed and
**all 14,770 trans families are usable**.

Reported per module: observed co-localization fraction, permutation-null
expectation, effect size, and an empirical p-value.

**Default metric** (configurable, see Open questions for the empirical pass
that should confirm it): a member pair counts as co-localized in a strain when
both genes are on the **same contig**, **more than 10 genes apart** in rank
order (clearing the `trans` definition's own k=10 window — see constraint 1),
and **within 100 kb**. Exposed as `--pangenome_neighborhood_max_kb` (default
100) and `--pangenome_neighborhood_min_gene_gap` (default 11, i.e. strictly
outside k=10).

**Two constraints this statistic must respect:**

1. **It must operate at a larger scale than the `trans` definition itself, or
   it is circular.** `pangenome_pair_classification.py` classifies a pair as
   `trans` only when its linkage fraction over a **k=10 gene window** is
   ≤ 0.05 — i.e. immediate adjacency is already excluded by construction.
   The neighborhood test must therefore measure at a deliberately coarser
   scale (e.g. "within 100 kb on the same contig", or "same contig, >10 genes
   apart"), and the chosen scale must be stated in the output, not left
   implicit.
2. **"Different contig" ≠ "far apart."** With a median of 283 scaffolds per
   assembly (range 6–3,035), two genes on different contigs may well be
   neighbors in the real chromosome. Cross-contig pairs must be scored as
   **uninformative/excluded**, never as "distant" — otherwise the statistic
   measures assembly fragmentation. The fraction of pairs excluded for this
   reason must be reported alongside the result.

### B2 — Secondary: reference-anchored view, explicitly partial

A scaffold-track view placing families on one reference strain's coordinates
(RS, best available at 19.0% coverage), **with the unplaceable fraction shown
as an explicit count in the figure**, not omitted silently. Labeled as
illustrative for the placeable subset, not as a genome-wide picture.

This is a scaffold view, not a chromosome view: RS is 349 scaffolds, Silveira
202. There is no chromosome-level assembly in this dataset.

## Measured constraints (why B is designed this way)

Measured on the real `genus_vs_ureesii` run, 2026-09-19:

| Question | Measurement |
|---|---|
| Trans families placeable on best single reference (RS) | 2,805 / 14,770 = **19.0%** |
| Silveira / UTAH_20380X16 / UTAH_20380X10 / Phoenix_2 | 17.6–17.8% |
| Union of 5 best-assembled strains | 6,229 / 14,770 = **42.2%** |
| Mean strains carrying a trans family | 93.4 / 530 |
| Trans families in <15% of strains (cloud) | 9,746 / 14,770 = **66%** |
| Assembly fragmentation | median 283 scaffolds (6–3,035) |

**A better reference assembly would not fix the coverage ceiling.** The limit
is gene content, not contiguity: a single genome can only place families that
strain actually carries, and 66% of trans families are cloud-frequency. A
telomere-to-telomere reference would still land near 19%. The placeable
subset is also *biased* toward commoner shell families — precisely the least
dynamic portion of the data — so a reference-only view would quietly
under-represent the signal of interest.

The 5-reference union reaches 42.2%, but each family is then placed in a
different strain's coordinate frame, so positions are not mutually
comparable. Rejected for that reason.

## Shared infrastructure

- **One annotation lookup**, `family → (pfam_class, module_id)`, ~15k rows,
  embedded once and referenced by both views.
- **The payload is valid only for the run that produced it.** Gene-family IDs
  are *not* stable across pipeline re-runs: a family's ID is its mmseqs
  cluster representative, and mmseqs does not pick the same representative
  twice. Measured on two runs of the same study with the same config and
  data: 54,421 vs 54,410 families (a 0.02% difference in the clustering) but
  only **73% of family IDs shared**. So the page must carry the run's own
  identity, and nothing may join this payload to another run's output by
  family ID. See `.living/learnings.md`, 2026-09-19.
- **One categorical palette** for Pfam classes and one for module IDs,
  defined once as a JSON legend object so the two views are visually
  consistent.
- **Pfam class assignment** is a curated keyword→class map (e.g. `NACHT`,
  `HET`, `Ank_*` → NLR/incompatibility; `ketoacyl-synt`, `AMP-binding`,
  `PP-binding`, `Thioesterase` → secondary metabolite; `MFS_1`, `Sugar_tr`,
  `ABC_tran` → transporter; else "other"; no domain hit → "unannotated").
  The map lives in one place in `lib/`, is data not code, and a family with
  domains in two classes takes the first match in a documented precedence
  order. A **family's "dominant Pfam class"** (used for the sidebar chip) is
  the most frequent class among its domains, ties broken by that same
  precedence order.

## B1 implementation and empirical pass (2026-09-25, issue #182)

Implemented in `lib/pangenome_neighborhood.py`, `bin/pangenome_neighborhood.py` and
the `MODULE_NEIGHBORHOOD` process; `report.md` gets one table. Two changes from the
text above:

- **Two statistics, each with its own null.** The fraction above is conditional on
  "same contig" (`obs_frac`). Its effect ratio is capped at 1/null, and on
  fragmented assemblies (many contigs about 100 kb long) almost every same-contig
  pair is within 100 kb. So the output also carries colocalized / all pairs
  (`*_total`). It is not a fragmentation measure: the null comes from the same
  strain and so has the same contig breaks.
- **Null pool.** `accessory` (genes of non-core families, the default) or `all`.

Empirical pass: 529-strain Coccidioides `rescue_structural_genus_vs_ureesii` run,
7 Leiden modules, 200 permutations, 25/50/100/200 kb x accessory/all, 4 min per
setting on 1 CPU (3.3 GB RSS).

| Module (families) | Effect, same-contig statistic, 25 / 50 / 100 / 200 kb (accessory null) | Effect, all pairs, 100 kb |
|---|---|---|
| 0 (8,070) | 0.97 / 0.99 / 1.00 / 1.00 | 1.04 |
| 1 (5,101) | 0.94 / 1.00 / 1.00 / 1.01 | 1.02 |
| 2 (2,184) | 1.16 / 1.52 / 1.27 / 1.09 | 1.63 |
| 3 (1,095) | 13.69 / 4.10 / 2.26 / 1.39 | 13.18 |
| 4 (321) | 0.86 / 0.82 / 0.89 / 0.96 | 0.84 |
| 5 (9) | 29.6 / 11.4 / 3.79 / 2.20 | 28.2 |

- The module calls are the same at every scale, so the 100 kb / 11-gene defaults
  stand. `min_gene_gap` was not swept.
- 91-99% of member pairs are on different contigs and are excluded.
- With the all-genes null, the non-clustered modules 0 and 1 show 1.2-1.3x at 25 kb;
  with the accessory null they stay at 0.94-0.97. The accessory null is the default.
- With millions of pairs, a 1.02-1.04x effect still reaches the permutation floor
  p = 1/201. Read the effect ratio first.

## Output format and storage

One self-contained HTML page, data embedded as JSON, drawn client-side
(canvas/SVG) — matching the existing novelty/loss report pattern
(`lib/report_template.py`, `lib/report_data.py`).

- **No per-island image files.** 50–100 islands as separate PNG/SVG files is
  the storage pattern this design exists to avoid.
- Presence arrays bitpacked: 530 strains × 56 genes ≈ 30 kB raw per island;
  ~100 islands stays under 1 MB.
- Target: whole page well under 10 MB.

This also respects `NovInvenio_Investigations`' class-3 rule — a regenerable
report artifact is published as a release asset, never committed; keeping it
to one file makes that mechanism cheap.

## Integration

New `bin/` scripts plus a Nextflow module following the pattern established by
`modules/pangenome/trans_modules.nf`:

- Data prep scripts emit compact JSON, not figures.
- **Inputs are compressed.** As of issue #111 the pipeline publishes
  `family_positions.tsv.zst`, `gene_positions.tsv.zst`,
  `pair_classification.tsv.zst` and `cooccurring_pairs.tsv.zst`. Every new
  reader uses `lib/compressed_io.py`'s `open_maybe_compressed()`, which
  handles plain/`.gz`/`.zst` alike -- and note that its `zstd -dc -f` flag is
  load-bearing under Nextflow, because `zstd` silently returns an *empty
  stream* for the symlinked inputs Nextflow stages. A reader that shells out
  to a CLI on its own needs the same care, and a symlink test.
- The HTML assembly step embeds the JSON payload.
- Gated with the existing `--pangenome_island_pfam_hmm` block, which is where
  domain annotation (`pfam.domtblout`) becomes available. View B1's statistics
  do not strictly need domains and could run earlier, but keeping both views
  on one page argues for one gate.
- Must degrade gracefully when there are zero trans modules (the
  single-species-ingroup case, already handled this way in
  `pangenome_detect_trans_modules.py`): emit a valid page stating so, not a
  crash.

## Testing

- Unit: row-collapsing/haplotype grouping, bitpacking round-trip, the
  permutation null (fixed seed, synthetic genome with known clustering), and
  cross-contig exclusion accounting.
- Graceful-degradation: zero islands, zero trans modules, a module whose
  members never co-occur in any single strain.
- Real-data: regenerate against `genus_vs_ureesii` and confirm the reported
  reference-coverage figure matches the 19.0% measured here.
- JS: `node --check` on the page's script bodies, per the existing convention
  in `tests/test_report_templates.py`.

## Implementation phasing

The two views share only a palette and one lookup table; they answer
different questions from different data and are independently useful. Build
them as two changes, not one:

**Phase 1 — View A + shared infrastructure.** The Pfam class map, the
palette/legend object, the annotation lookup, and the island synteny page.
Self-contained, no new statistics, immediately useful. This phase is what an
implementation plan should cover first.

**Phase 2 — View B.** B1 introduces a permutation null that needs its own
validation (fixed-seed synthetic genome with known clustering, plus the
cross-contig exclusion accounting), and B2 depends on B1's placement work.
Treat it as a separate spec-to-plan cycle rather than bolting it onto
Phase 1's PR.

## Resolved during spec review (2026-09-19)

- **Neighborhood scale for B1.** The defaults in the body (100 kb,
  >10 genes apart) stand as the *starting* values, explicitly provisional:
  Phase 2's first task is the empirical pass that either confirms them or
  moves them, and whichever value ships must be printed in the output
  alongside the statistic. Shipping a tunable with an unvalidated default and
  no record of what was used is the failure mode
  `--hmm_presence_cov` already demonstrates in this repo.
- **Does B1 also belong in `report.md`?** Yes -- one compact table (module,
  observed co-localization fraction, null expectation, effect size,
  empirical p, fraction of pairs excluded as cross-contig). `report.md` is
  what gets read on the cluster without copying a file off it, and the table
  is small. The HTML page remains where the per-module detail lives.

## Open questions

None outstanding. Phase 2's neighborhood-scale default is a task, not an
unresolved design question.
