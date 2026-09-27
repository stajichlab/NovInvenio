# Island locus view: exemplar-anchored synteny and shared indel breakpoints

Status: **reviewed 2026-09-24; approved as a starting design with the `F_min = 3`
flank fallback (section 2).** Section 8 (clinker synteny panel) added 2026-09-26
and approved as a starting design ("let's see how it looks"). Section 4b (DNA
presence check) added 2026-09-26 and approved; it changes the row classes
(section 5) and the ranking (section 6). Nothing here is implemented.
Replaces "View A" of `2026-09-19-pangenome-gainloss-visualization-design.md` as the
target design; the current page (`island_synteny.html`) stays until this ships.

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

If no strain meets (2) at `F`, retry (2) at `F_min` genes per side (default
`F_min = 3`; reviewer decision 2026-09-24). This keeps a flank anchor for a locus
that sits near a contig end in every strain, where 5 genes are not available.
The locus then uses `F_min` flank columns, and the page marks it "short flanks".
If no strain meets (2) even at `F_min`, use the one with the most flank genes and
mark the locus "exemplar at contig end" on the page.

### 3. Columns

`F` flank genes left + the locus genes + `F` flank genes right, in the exemplar's
gene order (`family_positions` rank). Default `F = 5`, or `F_min = 3` for a locus
whose exemplar was chosen at the `F_min` tier (section 2). The per-strain "flanks
intact" test in section 4 uses the same flank columns the locus shows. Each column carries family
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

### 4b. DNA presence check (added 2026-09-26)

**Why.** Sections 4-5 test annotated genes only. They cannot tell a deletion from
a gene-model or annotation difference. The planning spot check (2026-09-26,
plan Task 14) found 0 of 6 empty-site calls on the top 3 Coccidioides loci were
deletions: each "empty" strain carried the locus DNA at 97.0-99.9% identity
(blastn). At L001 one strain has two gene models (277 aa and 446 aa, two
families) where the other has one 705 aa model over the same DNA (a third
family).

**Strains checked.** For each drawn locus: every flank-intact strain with at
least one locus column not "in place" (the strains that would be empty site or
partial).

**Sequences.**
- Query: the exemplar's locus DNA, from the start of its first locus gene to the
  end of its last locus gene, with each exemplar locus gene's coordinates.
- Target: the strain's DNA on its flank contig, between the inner edges of its
  innermost in-place left-flank gene and innermost in-place right-flank gene.
  When the target is shorter than 50 bp, every missing column is DNA absent
  without an alignment.

**Alignment.** `blastn -task megablast`, query against target (subject mode), one
call per (locus, strain). A locus column is **DNA present** in the strain when
HSPs with identity >= `--locus_dna_min_id` (default 90%) cover >= 
`--locus_dna_min_cov` (default 80%) of that exemplar gene's span. Both defaults
are chosen, not validated.

**Cell states (adds to section 4).** An "absent" or "elsewhere" locus cell in a
checked strain becomes one of:
- **absent, DNA present**: gene-model or annotation difference. Drawn hatched grey.
- **absent, DNA absent**: the site lacks the gene's DNA. Drawn as absent.
"In place" and "contig break" cells are unchanged. Strains that are not checked
keep the section 4 states.

**Model difference (adds to section 5).** Row classes use the DNA states:
- **empty site**: flanks intact, >= 80% of locus columns *DNA absent*.
- **partial**: flanks intact, some columns in place and some *DNA absent*.
- **model difference**: flanks intact, and every locus column is in place or
  *absent, DNA present*, with at least one of the latter. Not a deletion.
- **full locus** and **uninformative** as in section 5.
The breakpoint track counts only changes between "in place" and "DNA absent".

**Ranking (changes section 6).** "Informative polymorphism" uses DNA-confirmed
empty-site strains: at least 10 empty-site strains (DNA absent) and at least 2
full-locus strains, ranked by the smaller count. The sidebar also shows each
locus's model-difference count.

**Pipeline.** New process `ISLAND_DNA_CHECK` between the locus computation and
the page. Inputs: the per-strain regions from the locus step and the genome
FASTA files from `--pangenome_data_dir` (`dna/`). Loci are batched per task,
following the HPCC rule of about 1-1.5 h of work per task; the batch size is set
from a measured run, not guessed. Switch: `--pangenome_locus_dna_check`
(default true). With false, the section 4-5 states are used and the page says
that "empty site" is not DNA-confirmed.

### 5. Row classes and breakpoints

Per strain, over the column states:

- **flanks intact**: some left-flank and some right-flank family are "in place" on
  the same contig, within `locus length + 2k` ranks of each other.
- **full locus**: flanks intact, all locus columns in place.
- **empty site**: flanks intact, >= 80% of locus columns absent. This is the
  deletion (or pre-insertion) state. With the DNA check on (section 4b),
  "absent" here means *DNA absent*, and a model-difference class is added.
- **partial**: flanks intact, some locus columns present, some absent. Its
  boundaries are **breakpoints**.
- **uninformative**: flanks not intact (contig break or missing flanks).

Rows = strains collapsed into identical state vectors (as today), grouped by row
class, then by species. Above the grid, a **breakpoint track**: for each boundary
between adjacent columns, the number of flank-intact strains whose state changes
between "in place" and "absent" there, split by species. A tall bar is a
**shared breakpoint**; a bar present in only one species is a lineage-specific
event.

### 6. Which loci to show

Default ranking: loci with the most **informative polymorphism**: at least 10
flank-intact strains in the empty-site class **and** at least 2 in the full-locus
class, ranked by the smaller of those two counts. With the DNA check on
(section 4b), only DNA-confirmed empty-site strains count. Alternative sorts: strain
count, size (today's order), and a text search. `--top_loci` (default 50).

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

### 8. Clinker synteny panel (added 2026-09-26)

The grid shows presence and position states. It does not show the genes
themselves. For each drawn locus, a clinker figure shows the real gene
neighbourhoods of a few strains side by side, with links between similar genes.

**Loci.** The loci the page draws (`--top_loci`, default 50).

**Strains** (at most `--clinker_max_strains`, default 12), chosen in this order:
1. The locus exemplar (section 2).
2. For each row class (full locus, partial, empty site, model difference) and each species: the
   flank-intact strain with the highest N50 (`assembly_quality_vs_content.tsv`),
   ties by strain name.
3. More full-locus or partial strains by N50, until the cap is reached.

Uninformative strains (flanks not intact) are never chosen, because their
region cannot be anchored.

**Region per strain.** On the contig where the strain's flanks are intact
(section 5): from the lowest to the highest rank among the strain's in-place
column copies, then `F` genes further on each side (`F_min` for a short-flank
locus), clipped at the contig end. For an empty-site strain the region is the
left-flank block plus the right-flank block, so the gap between them shows the
deletion. Regions differ in length between strains. The page states each
strain's region (contig:start-end, gene count).

**GenBank slices.** New process `ISLAND_GBK_SLICE`, one task for all drawn loci.
Inputs: the study genome FASTA, GFF3 and protein FASTA (`--pangenome_data_dir`),
`gene_positions.tsv.zst`, the tier-1 cluster TSV, and the per-strain regions
computed by the locus-view step. Output: one `.gbk` per (locus, strain). Each CDS
carries `/locus_tag` = the gene ID, `/translation` from the protein FASTA, and
`/note="family=<tier-1 family>"`.

**Clinker.** New process `ISLAND_CLINKER`, one task per locus:
`clinker <locus>/*.gbk -gf <locus>.groups.csv -p <locus>.html` with
gamcil/clinker **0.0.32 from PyPI** (added to the pixi environment as a PyPI
dependency, and to the container). Bioconda's `clinker` 1.33 is a different,
unrelated RNA-seq tool; do not use it. `-gf` takes a gene-to-group file, so the
tier-1 family is the group and clinker colours match the grid columns (verified
in the spike). After clinker runs, the embedded data are slimmed: the
`sequence` and `translation` fields are removed from every gene object, in
clusters and in links. clinker copies both genes' full sequences into every
link, which is what makes the file large.

**Page.** `island_synteny.html` gets a "Synteny (clinker)" panel below the
grid. On locus selection it loads `clinker/<locus_id>.html` in an iframe, and it
lists the strains shown, the reason each was chosen, and each region. A locus
without a clinker file (for example, the step was skipped) shows "no synteny
figure for this locus". The `clinker/` files are release assets (class 3), like
`island_synteny.html`. `bin/sync_pangenome_report.py` (NII) copies `clinker/`
with the page.

**Switch.** `--pangenome_clinker true|false` (default true). With false, both
processes are skipped and the panel says the step was not run.

**Cost and limits.** clinker compares every pair of clusters, so cost grows with
the square of the strain count; the cap keeps it at most 66 pairs per locus.
Measured by the 2026-09-26 spike (below): about 30 s and 104 MB per locus, so
50 loci are about 25 min of single-task CPU.

**Clinker spike (measured 2026-09-26, throwaway, `$SCRATCH/clinker_spike/`).**
Run `rescue_freqpol_immitis_in_posadasii_out`. Island
`407-0_S_OLD_CPA0002:scaffold_30:125276-129229` (231 strains, 5 member families,
Prenyltrans domain). 11 strains (10 occupied, 1 empty site; both species), 12-16
genes and about 32-40 kb per region, 5 flank genes each side. GenBank slices
built from the funannotate GFF3, genome FASTA and protein FASTA.
- clinker 0.0.32, 6 cores: 31.4 s wall and 104 MB peak RSS with default
  grouping; 25.1 s with `-gf` family grouping. `-gf` grouping works.
- HTML: 4.13 MB, fully self-contained (no external script or style URLs), so it
  works on the static site. Of that, 3.45 MB is 668 links that each embed both
  genes' DNA and protein sequences. Removing `sequence` and `translation`
  everywhere gives 0.81 MB. It is not yet checked that the slimmed page renders
  and behaves the same; the plan must check this in a browser before relying on
  it. 50 slimmed loci would be about 40 MB of release assets per run.
- Only 3 of the island's 5 member families sit together in any strain (a 2-3
  gene, 2-4 kb cassette); the other 2 are near-singleton families co-listed from
  one genome. The locus grouping (section 1) and in-place test (section 4)
  should drop such columns from a strain's region; the plan checks this on real
  loci.
- One empty-site candidate (B12471) had only 5 genes on its contig, too few for
  5-gene flanks, and was skipped. This is the `F_min` / uninformative case.

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
   With section 4b, this becomes a regression test: the 6 planning spot-check
   cases (L001-L003) must all be classed "absent, DNA present" / model
   difference. Also report how many empty-site calls in the top 50 loci of the
   Coccidioides run move to model difference, and the DNA check's run time.
4. Clinker panel: for 3 loci, check by eye that the exemplar's clinker track
   shows the same genes, in the same order, as the grid's columns, and that an
   empty-site strain shows its flanks joined with no island genes between them.

## Phase 2 (not in this design)

- **Order and orientation.** gene_positions has no strand column; rank order
  alone can flag reversed runs ("in place, reversed") but not strand. Needs strand
  from GFF3.
- **Polarity of indels** (insertion vs deletion) needs a strain tree; out of scope,
  as in the 2026-09-19 spec.
- **Joining across contigs** of one strain is not attempted: cross-contig stays
  "uninformative".

## Open questions for review

1. `F = 5` flank genes (`F_min = 3` fallback near contig ends, decided 2026-09-24), `k = 10` window, empty-site threshold 80%, locus grouping
   at 50% containment: defaults from the feasibility run, not swept.
2. Default ranking: informative polymorphism (section 6) vs strain count.
3. Keep the current `island_synteny.html` as a second page, or replace it?
4. Should the breakpoint track also count "contig break" boundaries, as a
   separate colour, to make assembly effects visible?
