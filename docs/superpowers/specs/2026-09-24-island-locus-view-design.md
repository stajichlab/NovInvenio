# Island locus view: exemplar-anchored synteny and shared indel breakpoints

Status: **approved by PI** (2026-09-24), except the > 2 species rule (section 6). Nothing here is implemented.
Replaces "View A" of `2026-09-19-pangenome-gainloss-visualization-design.md` as the
target design.

**Decisions (PI, 2026-09-24):** `F = 5` flank genes; `k = 10` gene window, but a
locus larger than `k` must still work (handled by chaining, section 5); keep the
current page's code, renamed as a deprecated page (section 8); break types get
different colours (section 5).

## Goal

For an accessory locus, show **where insertions/deletions start and end, and which
strains share them**. There is no reference genome, so every view is anchored on
one **exemplar genome's locus**: its genes, in its order, with conserved flanking
genes as anchors. Each other strain is then described *relative to that locus*:
does it have the flanks, does it have the island genes **at that place**, and where
does its content change.

## What the current page shows, and why it cannot answer this

Measured on the 529-strain Coccidioides `rescue_freqpol_immitis_in_posadasii_out`
run (NovInvenio `91e3157`), with
`NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/analysis/island_locus_view_feasibility.py`:

1. **Columns are only the island's non-core genes.** An island is a maximal run of
   consecutive non-core genes on one contig in one strain
   (`bin/pangenome_build_islands.py`), so the core genes around it are never
   drawn. Without flanks there is no anchor, and a boundary between present and
   absent cannot be placed.
2. **A filled cell means "present somewhere in that genome".** It is the presence
   matrix, not a positional test. (In practice a present copy is almost always at
   the locus: 0.4-0.9% of present cells fail the neighbour test below. But the
   page cannot show that.)
3. **Islands are per-strain variants of the same locus.** 10,579 of 11,880 located
   islands (89.0%) have >= 50% of their families inside a larger located island.
   The page lists them as separate entries.
4. **The exemplar is arbitrary.** `example_strain` is the first strain found with
   that exact member set. For the top 50 islands by size, 21/50 exemplars have the
   island touching a contig end on at least one side.
5. **Ranking by size picks the least informative loci.** The top 50 by size are
   mostly 2-strain islands. Their median number of strains with intact flanks is
   61 of 529, and none shows a clean shared indel (see Feasibility).

The 2026-09-19 spec put "present but rearranged vs present in expected order" out
of scope because the data model has no synteny-block assignment. This design does
not add one. It uses a stated, per-strain **neighbour test** on each strain's own
gene ranks (below), and it scores cross-contig cases as uninformative, following
that spec's constraint "different contig != far apart".

## Design

### 1. Loci, not islands

Group located islands into **loci**: an island joins the locus of a larger island
when >= 50% of its families are in that island (the M4 criterion). A locus entry
lists its variant islands (count, strains). The view draws loci; the sidebar shows
the number of variants.

### 2. Exemplar choice

For each locus, the exemplar is the strain that:
1. carries the locus's largest variant, and
2. has at least `F` genes on the same contig on **both** sides of it (interior, not
   at a contig end);
3. ties: highest N50 (`assembly_quality_vs_content.tsv`), then strain name.

If no strain meets (2), use the one with the most flank genes and mark the locus
"exemplar at contig end" on the page.

### 3. Columns

`F` flank genes left + the locus genes + `F` flank genes right, in the exemplar's
gene order (`family_positions` rank). Default `F = 5`. Each column carries family
ID, frequency bin (core/shell/...), Pfam class, and the exemplar's location
(contig:start-end, protein), as PR #179 already does. Flank columns are visually
marked as anchors.

### 4. Cell states (per strain, per column)

Positions are compared **only within one strain's own assembly** (ranks on that
strain's contigs). Let `k` = the gene window (default 10, same as `pair_class_k`).

| state | rule |
|---|---|
| **in place** | the strain has a copy of the family on contig X at rank r, and a copy of at least one other column family on the same contig X within `k` ranks |
| **elsewhere** | present (annotated or rescued) but no such neighbour |
| **rescue** | as "in place" or "elsewhere", but the strain's only copy is a `genome_only` TBLASTN hit (drawn hatched; no annotated gene) |
| **absent** | no copy |
| **contig break** | absent or elsewhere, **and** the strain's nearest in-place column copy is within `k` ranks of that contig's end: the locus runs off the assembly, so absence is not evidence |

"In place" is a neighbourhood test, not an order test. Order and orientation are
Phase 2 (below).

### 5. Row classes and breakpoints

**Chaining (why the window does not depend on locus size).** `k` is only used
between neighbouring genes, never across the whole locus. Per strain, take the
"in place" copies of the column families, sort them by that strain's own rank on
each contig, and walk them in order. Consecutive chained copies with a gap of

- `<= k` genes: normal adjacency;
- `k < gap <= G_max` genes (default `G_max = 100`): an **insertion** of `gap - 1`
  genes that the exemplar does not have at that point (drawn as a marker between
  the two columns, with its size);
- `> G_max`: a **break** in the chain.

This way a locus of any length works, and a strain carrying a large extra
insertion is shown as such instead of being dropped as "flanks not intact"
(the first draft's fixed `locus length + 2k` span would have done that).

Per strain, over the column states:

- **flanks intact**: one chain on a single contig contains at least one
  left-flank and one right-flank column.
- **full locus**: flanks intact, all locus columns in place.
- **empty site**: flanks intact, >= 80% of locus columns absent. This is the
  deletion (or pre-insertion) state.
- **partial**: flanks intact, some locus columns present, some absent. Its
  boundaries are **breakpoints**.
- **uninformative**: flanks not intact (contig break or missing flanks).

Rows = strains collapsed into identical state vectors (as today), grouped by row
class, then by species. Above the grid, a **breakpoint track**: for each boundary
between adjacent columns, stacked counts of strains by **break type**, each in its
own colour (colours from `lib/skins.py`, like every other page colour):

| break type | meaning | reading |
|---|---|---|
| **indel edge** | flank-intact strain changes between in place and absent here | a shared bar is a shared insertion/deletion boundary |
| **insertion** | chain gap of `k < gap <= G_max` genes here | strain-specific extra content at this point |
| **contig break** | the strain's chain ends at a contig end here | assembly, not biology |
| **rearrangement** | chain break (`gap > G_max`) or the family is "elsewhere" | moved or long-range change; order not tested (Phase 2) |

Each type is also split by species (e.g. hatched vs solid), so a bar that occurs
in only one species reads as a lineage-specific event.

### 6. Which loci to show (ranking)

Measured on the top 200 loci by strain count (same run; `F = 5`, `k = 10`,
empty site = >= 80% of locus columns absent; e = empty-site strains, f =
full-locus strains, among flank-intact strains):

- 154/200 have >= 10 empty-site and >= 2 full-locus strains.
- Widely carried loci are small: median 3 locus columns, max 9.
- 128/200 differ between species at Fisher p < 0.05/200 (species x empty/full).
  Strains are not independent (clonal structure), so this p-value is used only
  to **rank**, never reported as a test of significance.
- 20/200 are near-fixed differences (one species >= 95% full, the other >= 95%
  empty, >= 20 strains each); 36/200 are polymorphic in **both** species
  (>= 10 empty and >= 10 full in each).

Candidate rankings and what they surface (top examples):

| ranking | score | top example (e / f per species) | surfaces |
|---|---|---|---|
| A balanced polymorphism | min(e, f) | immitis 155e/10f, posadasii 90e/263f | common, polymorphic loci; mixes species-fixed and within-species patterns |
| **B species-differentiated** | Fisher p, species x (e, f) | immitis 18e/145f, posadasii 355e/0f | indels that separate *C. immitis* from *C. posadasii* (lineage-specific) |
| **E within-species polymorphism** | sum over species of min(e, f) | immitis 61e/53f, posadasii 166e/178f | the same indel polymorphic inside both species (shared among strains, not species-sorted) |
| C most anchored | flank-intact strains | 527 anchored, 332 partial | well-assembled loci, often complex (many partial rows) |
| D largest | locus columns | 9 columns, 443 anchored | big loci; fewer anchored strains |

Top-50 lists of A and B share only 20 loci, so the choice matters.

**Decided (PI, 2026-09-24):** two tabs, because the study's main comparison is
reciprocal *C. immitis* vs *C. posadasii*. A, C, D and a text search stay as sort
options, not tabs. The page shows e/f per species for every locus, so the reader
sees why it ranked where it did.

| setting | decision |
|---|---|
| default tab | **B "Differs between species"** |
| locus scored only when | >= 10 flank-intact strains in **each** species |
| tab B ranking | **effect size**: difference between species in the fraction of flank-intact strains with the full locus. Fisher p shown for information only (strains are clonal, so p is not a test) |
| tab B filter | "near-fixed only" toggle: one species >= 95% full, the other >= 95% empty |
| tab E "polymorphic" in a species | >= 10 empty-site **and** >= 10 full-locus strains in that species |
| tab E scope | polymorphic in **either** species (default), toggle for **both** |
| tab E ranking | sum over species of min(empty-site, full-locus) |
| C "most anchored" | sort option only |
| > 2 species (e.g. B's score) | **not decided**; not needed for this study. Options: largest pairwise difference, or each group vs the rest |

### 7. Page layout

- Sidebar: loci with size, variant count, empty-site / full / partial /
  uninformative strain counts, Pfam class chip.
- Main: title = exemplar locus (strain:contig:start-end); note line stating the
  exemplar rule and `F`/`k`; breakpoint track; column header (labels + class
  strip + anchor marks); grid; legend for the five states.
- Popups: column = family, bin, domains, exemplar location; cell = strain(s),
  state, that strain's own contig:rank (and bp when annotated), and why the state
  was assigned ("neighbour <family> at rank r±d").
- Keep the diagnostics banner at the top. Add one line under it stating that
  accessory content is confounded with assembly quality on runs where
  `assembly_quality_confound` triggered, so "uninformative" rows are expected.

### 8. The current page is kept, renamed

The current page's code (`lib/island_synteny.py`, `lib/island_synteny_template.py`,
`bin/pangenome_island_synteny.py`, the `ISLAND_SYNTENY` process) is not deleted.
Its output is renamed `island_presence_grid.html` ("deprecated: presence anywhere
in the genome, no positional test") and it keeps running behind a param
(`--pangenome_legacy_island_grid`, default `false` once the new page ships). The new
page is `island_locus.html`. The run report and site pages link to the new page;
the legacy page is linked only when it was produced.

### 9. "Core" must not come from the ingroup alone

`frequency_table.tsv` bins count only the dereplicated **ingroup**
representatives (in `rescue_freqpol_immitis_in_posadasii_out`: 133 *C. immitis*
strains; a *C. posadasii*-only family is binned "singleton" even when 98% of
*C. posadasii* strains carry it). `BUILD_ISLANDS` uses those bins for its
core/non-core split in **every** strain, so in outgroup strains the islands
include genes that are core in their own species. For this view:

- flank anchors ("core" in section 3) are taken from a frequency computed over
  **all** strains of the run (or per species), not from `frequency_table.tsv`;
- the per-species e/f counts (section 6) are computed from the matrix directly,
  never from the bins;
- whether `BUILD_ISLANDS` itself should use an all-strain core is a separate
  pipeline question, to be raised as an issue, not changed inside this view.

Found by the HrmA analysis (NII `studies/fungi/coccidioides_pangenome/analysis/
2026-09-24-hrmA-presence-and-neighbourhoods.md`), confirmed from the run
(max `strain_count` = 133; `5 / 0.0376 = 133`).

## Data and wiring

All inputs already exist per run: `islands_with_domains.tsv`,
`family_positions.tsv.zst` (contig, rank per strain/family/copy; includes rescue
positions), `gene_positions.tsv.zst`, tier-1 cluster TSV, `rescue_positions.tsv`,
`presence_matrix.rescued.tsv` (present vs genome_only), `frequency_table.tsv`,
`assembly_quality_vs_content.tsv` (N50, for exemplar ties). Contig end in genes =
max rank per (strain, contig), from `family_positions`.

Computation is per locus, streamed and filtered to the drawn loci's families (the
pattern `pangenome_island_synteny.py` uses since #118). The feasibility script
does the full state computation for 50 loci x 529 strains in 59 s and 90 MB
(single thread, pure Python), so no new heavy step is needed.

Payload: state vectors as 3-bit codes per column, stored once per collapsed
pattern; 50 loci x ~40 columns x a few hundred patterns stays well under the
current page's 0.8 MB.

## Feasibility (measured, 2026-09-24)

Run `rescue_freqpol_immitis_in_posadasii_out`, `F = 5`, `k = 10`, top 50 qualifying
islands (>= 2 strains). The design's neighbour and flank rules as in sections 4-5.

| measure | top 50 by size (today) | top 50 by strain count |
|---|---|---|
| exemplar has >= 5 genes both sides | 21/50 | 30/50 |
| exemplar locus touches a contig end | 21/50 | 7/50 |
| median flank genes left / right | 6 / 7 | 15 / 16 |
| flank genes that are core or soft-core | 66% | 76% |
| locus cells: in place / elsewhere / absent | 29.8% / 0.9% / 69.3% | 69.3% / 0.4% / 30.3% |
| strains with intact flanks (median of 529) | 61 | 500.5 |
| flank-intact strains: empty site / full locus (median) | 0 / 2 | 52.5 / 245.5 |
| loci with >= 10 empty-site and >= 2 full-locus strains | 17/50 (80% rule), 0/50 (100% rule) | 39/50 (same with 80% or 100%) |

Reading: the flank-anchored view works on the loci that many strains carry;
there, most strains are anchored and 39 of 50 loci show a shared presence/absence
state with intact flanks on both sides. The size ranking the current page uses
selects mostly 2-strain loci, where most strains have no anchored flanks.

## Validation plan

1. Unit tests on synthetic strains: in place / elsewhere / absent / contig break,
   flank-intact, each row class, breakpoint counts.
2. Regression on the real run: the feasibility numbers above (within the
   feasibility script's definitions) as a check that the pipeline implementation
   matches the prototype.
3. Spot check 5 "empty site" and 5 "partial" calls by sequence: extract the flank
   region from two strains' DNA and align (minimap2 or blastn) to confirm the
   site is contiguous in the empty-site strain. This tests the neighbour rule
   against the sequence, which the rule itself cannot.

## Phase 2 (not in this design)

- **Order and orientation.** gene_positions has no strand column; rank order
  alone can flag reversed runs ("in place, reversed") but not strand. Needs strand
  from GFF3.
- **Polarity of indels** (insertion vs deletion) needs a strain tree; out of scope,
  as in the 2026-09-19 spec.
- **Joining across contigs** of one strain is not attempted: cross-contig stays
  "uninformative".

## Open questions for review

1. Settled: `F = 5`, `k = 10` with chaining. Still defaults, not swept: `G_max = 100`,
   empty-site threshold 80%, locus grouping at 50% containment.
2. Settled: ranking and tabs in section 6.
3. Settled: current page kept as `island_presence_grid.html`, deprecated.
4. Settled: contig breaks are drawn, in their own colour.
