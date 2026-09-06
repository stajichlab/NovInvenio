# Targeted config-builder for NovInvenio (design)

Date: 2026-09-05
Status: proposed

## Problem

Building a NovInvenio analysis config CSV by hand is a bottleneck. The
actual usage pattern that motivates this tool is **many small, targeted
comparative-genomics runs**, not a handful of huge whole-clade runs:

- ~4 ingroup species (a focal species + a few close relatives, chosen
  for phylogenetic proximity and/or a shared trait/phenotype the study
  hypothesizes shares underlying genes) vs.
- ~5-8 outgroup species (often the *same* outgroup pool reused across
  many different focal-ingroup studies).

Concrete first batch: one run each with Mucor, Phycomyces, Basidiobolus,
and Spizellomyces as the focal ingroup species, each paired with a
shared outgroup pool, to both answer the biological question and assess
whether the method finds meaningful signal across independent test
cases.

Source data:
- `Fungi_BFD/samples.csv` (23,683 assemblies) — species + semicolon
  lineage string (PHYLUM;...;GENUS) per assembly, no direct genus-level
  branch-distance information.
- `Fungi_BFD_runs/genome_annotation/_reuse_assignments/repr_assignments.tsv`
  (18,000 rows: `out, species, is_representative, representative_out,
  ani_to_representative, reuse_eligible`) — the authoritative per-species
  representative-assembly flag. This replaces any BUSCO-based heuristic
  for strain collapse.

An earlier framing of this design (whole-clade novelty/loss studies
across hundreds of genomes, e.g. all of Mucoromycotina vs all of
Dikarya) turned out not to be the primary need; that framing's
clade-alias and bulk-downsampling ideas are kept as deferred/future
work (see "Deferred" below) since they may still matter for a
future large-scale study, but are not part of this v1.

## Goals

1. Given a focal species, suggest 3-4 sensible ingroup companions,
   ranked by phylogenetic proximity (shared lineage tokens), from the
   master pool — instead of scrolling a flat list of every species in
   the focal species' genus/family/order.
2. Let a companion also be pulled in because it shares a *trait* with
   the focal species. In v1 this means: among the same lineage-scoped
   candidate pool `nearest` already considers, prefer those sharing a
   specified trait. A trait shared only with something phylogenetically
   distant (e.g. an early-diverging fungus and an animal both having
   motile zoospores) is out of scope for the automatic mode and stays a
   manual, reasoned override (see "Trait mode" below).
3. Define an outgroup pool once, name it, and reuse it across many
   focal-ingroup studies.
4. Generate N config CSVs in one pass from a list of focal species +
   one shared outgroup pool (a "batch" or "sweep").
5. Collapse each selected species to its single representative
   assembly using `repr_assignments.tsv`, not guesswork.
6. Every choice (which relatives, which traits drove a pick, which
   outgroup pool) is recorded in a reusable, versioned file so a study
   can be regenerated identically later, and so a colleague can see
   *why* a given companion was picked.

## Non-goals (v1)

- Whole-clade / hundred-plus-genome studies and their associated
  subsampling and workload-estimation machinery (deferred).
- A GUI or TUI (deferred; see below).
- Animal proteome integration for the "shared with animals" variant of
  the chytrid study (needs its own non-Fungi_BFD data source; separate
  ticket).
- A general trait ontology or external trait-database integration, and
  an automated FUNGuild importer (start with a small hand-curated seed
  file; see "Deferred" below).
- `mode: trait` searching beyond a lineage-scoped candidate pool (see
  "Trait mode" below) — cross-phylogeny trait matches stay a manual
  `mode: explicit` override in v1.
- Phylogenetic-tree-backed relative ranking (lineage-string proximity
  is the v1 heuristic; a real tree is future work).

## Data model

**Master pool** (extends `bin/convert_bfd_samples.py`'s output): one
row per species (post strain-collapse), columns:
`Species, Strain, ProteinPath, DNAPath, Lineage, NCBI_TaxID`. `Short` is
*not* stored in the master pool — `convert_bfd_samples.py::make_short`
disambiguates collisions by numeric suffix in the order rows are seen,
so the same Short can denote a different genome across regenerations or
different `--order`/`--class-name` filters. The spec and this pool key
everything on **Species name** instead; a `Short` is assigned only at
render time (deterministically, from the master pool's full species
list, not per-batch) and the Species→Short mapping used is echoed in
the renderer's summary output and in a `.map.tsv` written alongside
each rendered config, so a run is traceable back to exactly which
genome a Short meant. `ProteinPath`/`DNAPath` are resolved absolute
paths (not basenames) — the master pool spans far more species than any
one batch touches, so symlinking every file into one shared `--data_dir`
up front doesn't scale; the renderer creates a small per-batch link
directory containing only the species actually selected for that batch.
`NCBI_TaxID` is carried through from `Fungi_BFD/samples.csv`'s
`NCBI_TAXONID` column (currently dropped by `convert_bfd_samples.py`)
into the config's optional report-only `NCBI_TaxID` column.

Built by:
1. Reading `Fungi_BFD/samples.csv` as today, keeping `NCBI_TAXONID` and
   resolving `Protein`/`DNA` to absolute paths instead of basenames.
2. Joining `repr_assignments.tsv` on its `out` column (the
   `<Species>_<Strain>` dirname `find_annotation()` already builds) to
   the matching assembly row; group by the `species` column and keep
   the one row where `is_representative == True`. Zero `True` rows *or*
   more than one `True` row for a species is a hard error (not silently
   picking the first) — as of 2026-09-05 the real table has neither
   case, but the check must not assume that stays true.

**Trait data** (`config_support/traits/`, new, hand-curated; scaffolded
2026-09-05, see that directory's `README.md`): split into
`trait_definitions.yaml` (the controlled vocabulary — each trait
declares its own fixed set of legal values, plus an optional
`ontology_term` provenance field citing PATO/FAO/GO/FYPO/OMP where a
real mapping is known) and `traits.csv` (long-format species-to-trait
rows: `Species,trait,value,source,notes`). A loader must hard-error on
any `(trait, value)` pair in `traits.csv` not declared in
`trait_definitions.yaml` — no near-miss string matching, and multiple
rows for the same species+trait are valid (e.g. a pathotroph that is
also a saprotroph). Nine traits are seeded: `hyphal_septation`,
`spore_motility`, `thermotolerance`, `light_response`, `habitat`, a
host-association family (`animal_association`, `plant_association`,
`fungal_association`, `algal_association` — sharing one convention:
`none` is a real assertion that never coexists with another value for
that species+trait, and pathogen host-breadth is read by checking the
pathogen-like value across all four rather than a fifth consolidated
trait), and `trophic_mode` (deliberately pinned to FUNGuild's own
Pathotroph/Saprotroph/Symbiotroph split rather than a finer
biotroph/necrotroph/hemibiotroph distinction FUNGuild can't support).
Starts covering only the species relevant to the first batch
(Mucor/Phycomyces/Basidiobolus/Spizellomyces neighborhoods); grows
trait-by-trait and study-by-study.

`config_support/traits/funguild_reference.csv` is a raw, deduplicated
reference dump (1,871 rows, pulled 2026-09-05) from the `funguild`
table in the Fungi_5k project's `functionalDB/function.duckdb` — a real
FUNGuild run over that assembly set. It is source material, not
curated input: `traits.csv` rows derived from it are transcribed by
hand today, labeled genus-level (FUNGuild's own reference database
resolves at genus level in most cases) rather than species-verified,
and a `bin/import_funguild_traits.py` to automate that transcription
(mapping FUNGuild's `guild`/`trophicMode` strings onto this
vocabulary) is future work, not part of this v1.

**Lineage-proximity ranking**: `Lineage` is root-first
(`Fungi;PHYLUM;SUBPHYLUM;CLASS;SUBCLASS;ORDER;FAMILY;GENUS`, with empty
ranks — commonly SUBPHYLUM/SUBCLASS — dropped from the string entirely
rather than kept as a placeholder), so proximity is the longest common
**prefix**, not suffix. Comparing raw list position is unsafe once a
rank is missing for one species but not the other (a same-genus match
could be misjudged as merely same-order if one species' string skips
SUBCLASS and the other doesn't shift the same way). The ranking
function must therefore align on the fixed rank *names*
(PHYLUM/SUBPHYLUM/CLASS/SUBCLASS/ORDER/FAMILY/GENUS), not string
position — reusing `make_config.py`'s lineage handling only for the
underlying species/lineage data structures, not its
`parent_segment`/`child_segment` helpers (those do token lookups for a
single named rank, not a general depth-of-agreement comparison, and
would need to be adapted rather than called directly). Two species tied
at the same rank of agreement are broken deterministically by sorting
on Species name.

**Ingroup/outgroup disjointness**: a nearest-relatives pick and a named
outgroup pool can legitimately overlap in candidate space (e.g. a
Spizellomyces study's "early-diverging" outgroup pool and its
nearest-relatives candidates are drawn from neighboring chytrid
lineages). The renderer excludes every named outgroup pool's members
from the candidate list for that study regardless of `mode`
(`nearest`, `trait`, or `explicit`), and errors out if the focal
species itself is listed inside the outgroup pool it's paired with.

## Selection spec format

A single YAML file per **batch** (a batch may contain just one study):

```yaml
batch: mucoromycota_focal_v1   # required; used to name rendered configs

outgroup_pools:
  early_diverging_v1:
    members: ["Rhizopus oryzae", ...]   # Species names, resolved against the master pool

studies:
  - focal: Mucor circinelloides
    ingroup_extra:
      mode: nearest        # nearest | trait | explicit — see "Trait mode" below
      n: 3
    outgroup_pool: early_diverging_v1
  - focal: Phycomyces blakesleeanus
    ingroup_extra: {mode: nearest, n: 3}
    outgroup_pool: early_diverging_v1
  - focal: Basidiobolus meristosporus
    ingroup_extra:
      mode: trait
      trait: animal_association
      value: pathogen
      n: 2
    outgroup_pool: early_diverging_v1
  - focal: Spizellomyces punctatus
    ingroup_extra: {mode: nearest, n: 3}
    outgroup_pool: early_diverging_v1
```

Every entry — `outgroup_pools[*].members`, `studies[*].focal`, and
`ingroup_extra.members` under `mode: explicit` — is keyed by **Species
name**, not Short, per the master-pool note above.

`mode: explicit` accepts a literal `members: [...]` list when neither
the auto-picker nor `mode: trait` gives the right answer and a human
overrides it by hand — this is how a cross-phylogeny trait-driven
choice (see "Trait mode" below) gets made in v1. An optional `reason:`
free-text string can accompany any `nearest`/`trait`/`explicit` entry
to record *why* a companion was chosen (required in practice for
`explicit`, since there's no algorithm to fall back on for "why" —
e.g. "shares motile-zoospore trait with focal, per
config_support/traits/traits.csv"). `TaxonGroup` in the rendered CSV,
for a `nearest` or `trait` pick, is the shared rank name (e.g.
`Mucoraceae`) that produced the match; for an `explicit` pick it is the
deepest lineage rank for that species (no comparison to a focal is
implied).

Every study renders with the plain `IN`/`OUT` GROUP scheme (focal +
companions = `IN`, outgroup pool = `OUT`) — v1 does not need the
`DISCOVERY_TARGET`/`NEAR_INGROUP`/`BROAD_OUTGROUP` multi-tier scheme,
since these are already small, single-tier ingroup-vs-outgroup studies.
That richer scheme stays available in NovInvenio for a future
whole-clade study if one is ever built on top of the deferred
clade-alias work below.

**Trait mode: `mode: trait` is `mode: nearest`, filtered.** It reuses
`nearest`'s candidate universe unchanged (species sharing some lineage
rank with the focal, scoped to a configurable rank — default `Order` —
not the whole master pool; this is what keeps "small candidate lists"
actually true rather than aspirational: "all saprotrophs" pool-wide
would be hundreds of species, "saprotrophs among Mucorales" stays
small by construction), then:

1. **Filter** that pool to candidates carrying the exact `(trait,
   value)` given (per the exact-equality match rule above) — a species
   present in the pool with no row, or a different value, for that
   trait is simply not a candidate. No partial credit.
2. **Rank** the survivors by the same lineage-proximity function
   `nearest` uses, with the same deterministic tiebreak.
3. Take the top `n`.

So a trait match *filters*, phylogenetic proximity *orders* — trait
matching never overrides phylogeny outright (a trait match three ranks
further away never outranks a closer relative that also has the
trait), and no separate tie-break rule is needed. If the filtered pool
is empty, the renderer errors out rather than silently falling back to
plain `nearest` — a silent fallback would hide that the trait pick
didn't actually work.

This means a companion sharing a trait *only* with something
phylogenetically distant (the original motivating example: an
early-diverging fungus and an animal both having motile zoospores) is
explicitly out of scope for `mode: trait` — that case stays a manual
`mode: explicit` pick with a `reason:` string citing the relevant
`traits.csv` row(s), same as any other human override.

## Renderer

A new script, `bin/build_targeted_configs.py` (sibling to
`make_config.py`, sharing its `MasterSample`/lineage-segment helpers
rather than duplicating them):

1. Load master pool + `trait_definitions.yaml` + `traits.csv`.
2. For each `studies[]` entry: resolve `focal`, compute the candidate
   list for the study's `mode` (`nearest`: lineage-ranked; `trait`:
   lineage-scoped pool filtered to the given `(trait, value)`, then
   lineage-ranked; `explicit`: the literal list), take the top `n` (or
   the explicit list), exclude the resolved `outgroup_pool`'s members
   from that candidate list first (per the disjointness rule above),
   and emit one `configs/<focal_short>_<batch_name>.csv`.
3. Print a per-study summary table (focal, chosen companions with the
   lineage rank or `(trait, value)` that justified each, outgroup pool
   size) so the pick is visible before any pipeline run starts — this
   is the stand-in for a live UI in v1.
4. Fail loudly (not silently skip or fall back) if a focal species, a
   named trait/value, a pool member isn't found in the master pool, or
   a `mode: trait` filter leaves an empty candidate pool.

## Workload estimate

At ~10-12 genomes per study, `--cluster_tool pairwise` is trivially
cheap (tens of pairwise jobs) — no per-study estimator is needed. The
renderer's batch summary reports total study count and total distinct
proteome count across the whole batch (for `search_cache/` sizing
awareness), which is sufficient for v1.

## Testing

- Golden-file test: a small fixture master pool + spec → expected
  config CSV(s), byte-stable (mirrors the existing
  `tests/test_config_parser.py` fixture style).
- Unit tests for the lineage-proximity ranking function (including a
  case where one species is missing a rank the other has) and the
  `repr_assignments.tsv` join (both the "no True row" and the "more
  than one True row" failure cases).
- One test per `mode` (`nearest`, `trait`, `explicit`), plus for
  `trait`: a case where the filter narrows a multi-candidate pool to
  the correct proximity-ranked subset, and a case where the filter
  leaves an empty pool and the renderer errors rather than falling back
  to `nearest`.
- A disjointness test: an outgroup pool member that would otherwise
  rank as a top nearest-relative (or trait-matching) candidate must be
  excluded from the ingroup pick, and a focal species listed inside its
  own outgroup pool must be a hard error.
- A `traits.csv`/`trait_definitions.yaml` loader test: an undeclared
  `(trait, value)` pair is a hard error, and a `none` value coexisting
  with another value for the same species+trait is a hard error.

## Deferred / future work

- Whole-clade studies (e.g. all of Mucoromycotina vs all of Dikarya):
  the clade-alias-with-citation and `representative_per_child`
  bulk-downsampling ideas from the earlier draft of this design remain
  a reasonable approach if/when that scale of study is actually
  needed; not built now.
- A `textual`-based TUI for interactively browsing/adjusting a batch
  spec (would edit the same YAML format, so no rework of the layer
  above).
- A browser-based GUI (phase 2 of the TUI, per the same reasoning).
- Real phylogenetic-tree-backed relative ranking (via `nf_phyling`).
- Animal-proteome pool for the chytrid/Blastocladiomycota-vs-animals
  variant.
- `bin/import_funguild_traits.py`: automating the transcription from
  `config_support/traits/funguild_reference.csv` into `traits.csv`
  (mapping FUNGuild's `guild` substrings onto the host-association
  traits, `trophicMode` onto `trophic_mode`) instead of hand-copying.
  Today's `traits.csv` rows derived from FUNGuild were transcribed by
  hand for four species only.
- `mode: trait` that searches the whole pool rather than a
  lineage-scoped one, which would make cross-phylogeny trait matches
  (the animal/chytrid motile-zoospore case) expressible as an automatic
  mode instead of a manual `explicit` override — deliberately not built
  now; see the "Trait mode" section above for why.
- A real trait ontology / external trait-database integration beyond
  the hand-curated seed file and the optional `ontology_term` citation
  field.
- Timing-calibrated (not just count-based) workload estimation, once
  real batch runs produce timing data to calibrate against.
