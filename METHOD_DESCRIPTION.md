# Method Description: Candidate Detection Strategies

This document is the single reference for *how* NovInvenio decides a gene is
present or absent in a proteome, and how that decision differs across the
three `--cluster_tool` pathways. It exists separately from `README.md` (which
covers running the pipeline) and `CLAUDE.md` (which covers repository layout
and conventions) because the presence/absence logic is the scientific core of
the tool and is easy to lose across three parallel implementations.

All three pathways answer the same underlying question — "is this gene
present in proteome X?" — but they source that evidence differently:
pairwise search score, cluster/HMM profile membership, or a discovery+screen
two-phase HMM search. Whichever pathway is used, the output contract is
identical: `presence_matrix.tsv` (protein × proteome 0/1) and
`candidates.txt` (`source_proteome::protein_id`), so downstream
clustering/TBLASTN/annotation/report stages don't need to know which pathway
produced them.

Every pathway is run **twice** per analysis, symmetrically: once with the
ingroup (or, for `novelty_discovery`, `DISCOVERY_TARGET`) as query, looking
for **novelties** (present in ingroup, absent from outgroup); once with the
outgroup as query, looking for **losses** (present in outgroup, absent from
ingroup). The mechanics below are described in the novelty direction; the
loss direction is the identical algorithm with the query/target group
swapped (`--query-group OUT`, `outgroup_min_frac`, `loss_ingroup_max_frac`).

---

## 1. `pairwise` (default) — all-vs-all search with paralog calibration

**Where:** `workflows/search.nf` → `bin/parse_self_hits.py` →
`bin/build_presence_matrix.py`

### Search step

Every ingroup proteome is searched, protein-by-protein, against every *other*
proteome (ingroup ∪ outgroup, excluding self) using one of `phmmer` /
`diamond blastp` / `blastp` (`--run_tool`). This is O(|ingroup| × (N−1))
pairwise jobs — the direct cost that motivates the `mmseqs` pathway below at
scale.

In parallel, each proteome is searched **against itself** (`-E 100
--max-target-seqs 2`), and the best hit that is *not* the query itself (the
"rank-2 hit") is recorded per protein as its within-proteome paralog.

### Why a paralog cutoff, not a flat E-value

A single global E-value cutoff (e.g. `1e-5`) makes the wrong call for two
different reasons at once: it is too loose for a protein with weak paralogs
(a hit just inside `1e-5` might be homologous to the *paralog*, not
specifically to the query), and too strict for a protein with strong,
divergent true orthologs. The self-hit calibration turns "is this cross-
species hit real" into a question answerable from each protein's own paralog
structure:

1. **Paralog-cutoff filter.** A cross-proteome hit only counts if its E-value
   is strictly better than the query protein's own paralog E-value (from the
   self-search). A protein with no detectable within-proteome paralog falls
   back to `--evalue` (default `1e-5`, `DEFAULT_EVALUE` in
   `build_presence_matrix.py`).
2. **Paralog-competition filter.** Even if filter 1 passes, the hit is
   disqualified if the query's own paralog scores *better* against the same
   target than the query does — i.e. the hit is better explained by the
   shared domain with the paralog than by the query specifically.
   `--paralog-competition-scope` controls the granularity:
   - `proteome` (default, strict): disqualify if the paralog's best hit
     *anywhere in the target proteome* beats the query's hit.
   - `target`: disqualify only if the paralog beats the query on the *same
     target protein*. Looser — preserves calls where the paralog also has a
     real, distinct ortholog in the target genome (documented example:
     HEX-1 vs its paralog eIF5A survives under `target` scope but is
     incorrectly dropped under `proteome` scope).

### Candidate rule

A protein (sourced from a `--query-group` proteome, default `IN`) becomes a
row in `presence_matrix.tsv` for every target proteome it passes both filters
against. It becomes a line in `candidates.txt` when:

```
query_frac >= ingroup_min_frac   (default 0.75)
AND
other_frac <= other_max_frac     (default 0.0 — strictly absent from every other-group proteome)
```

`other_max_frac` is the knob that turns a strict novelty/loss rule into a
"nearly missing" one (`--loss_ingroup_max_frac` for the loss direction). Note
this predicate only shapes `candidates.txt` — the matrix itself always keeps
every scored row, both novel and non-novel, so reports can recompute
different thresholds without re-running search.

### Genomic validation and annotation

Candidates are extracted, clustered with `mmseqs easy-cluster`
(`--min-seq-id 0.3 -c 0.8 --cov-mode 0`, purely to deduplicate near-identical
candidates before TBLASTN — not a presence call), and cluster representatives
are TBLASTN'd (`-evalue` = same `--evalue` param) against outgroup genomes.
In the novelty direction this runs with `--skip_tblastn_filter`: a TBLASTN
hit is reported as a column, not used to disqualify a candidate, because the
paralog-calibrated absence call is already trusted. Pfam (`hmmscan`) and
SwissProt (`diamond blastp`) annotate every row in the final matrix.

### Cost / sensitivity trade-off

Diamond tolerates the O(N²) job count; phmmer does not scale past small-to-
medium configs. The paralog calibration is per-protein and precise, but pays
for that precision with the full pairwise search cost — this is the reason
the `mmseqs` pathway exists (docs/adr/0002).

---

## 2. `mmseqs` — family-profile search (ADR-0002)

**Where:** `workflows/search.nf`'s `PROFILE_SEARCH` equivalent;
`modules/mmseqs_family_cluster.nf`, `modules/build_family_profiles.nf`,
`modules/family_hmmsearch.nf`; docs/adr/0002-family-profile-search-pathway.md

### Core idea

Collapse the redundant ingroup into gene families up front, and turn the
absence call into a single profile-HMM search per target proteome instead of
one search per protein-pair. This drops the job count from ~|IN|×(N−1)
pairwise searches to ~N hmmsearch jobs (one per proteome, searching the whole
family-HMM database at once).

### Pipeline

1. `mmseqs cluster` (cascaded, `-s 7`, `--min-seq-id 0.3 -c 0.8 --cov-mode 0
   --cluster-mode 0`) clusters the **whole seed group** (ingroup for novelty,
   outgroup for loss) into gene families. This is clustering as *seeding*,
   not as the candidate-cluster step used in the pairwise pathway (ADR-0002
   is explicit these are two different concepts sharing the word "cluster").
2. Families with ≥2 members are kept (`--family_min_members`,
   `--family_max_members`); singletons are dropped — a true single-copy
   ingroup-specific gene either gets recovered because a shattered-family
   partner shows up in the all-proteome search, or fails the presence-
   fraction filter as an orphan.
3. Each surviving family is aligned (`famsa`) and turned into an HMM
   (`hmmbuild -n <rep_id>`), scattered across `BUILD_CHUNK` tasks
   (`--family_chunk_size`, default 200 families/chunk) and merged into one
   `family_profiles.hmm`.
4. That single HMM database is `hmmsearch`'d against **every** proteome in
   the analysis (both ingroup and outgroup) — presence is evaluated
   uniformly across all species, not read off cluster membership for the
   seed group. This was a deliberate ADR-0002 resolution: reading ingroup
   presence from membership while reading outgroup presence from a search
   would silently deflate `ingroup_frac` relative to `outgroup_frac`.

### Presence cutoff

A family is called **present** in a proteome when, from `hmmsearch
--domtblout` (fixed `--Z`):

```
per-sequence E-value < hmm_presence_evalue   (default 1e-3)
AND
profile coverage       >= hmm_presence_cov    (default 0.5, i.e. >=50% of the HMM's length is aligned)
```

Both thresholds were chosen deliberately lenient-for-absence: the E-value
ceiling is loose enough not to miss a diverged true ortholog (a family
profile is more sensitive than any single representative), and the coverage
floor guards against the profile hitting only a small, promiscuous shared
domain rather than the whole gene. This is intentionally the mirror image of
the pairwise pathway's per-protein paralog calibration — precision comes
from clustering quality, not per-hit E-value competition. See ADR-0002 item
8 for the swept parameter range (E-value + coverage) and the planned
per-family GA/TC/NC gathering-threshold refinement (not yet implemented).

### Candidate rule

Same shape as the pairwise pathway, evaluated at the family level and then
expanded back to every member protein: a family is a novelty candidate when
present (by the rule above) in `>= ingroup_min_frac` of ingroup proteomes and
absent from every outgroup proteome. `presence_matrix.tsv` and
`candidates.txt` come out in the same shape as the pairwise pathway's, so all
downstream stages are unchanged.

### Validation is a filter here, not just an annotation

Because a clustering-based absence call is more artifact-prone (mis-merged
or mis-split families) than the pairwise pathway's per-protein paralog
calibration, TBLASTN validation is applied here as a **disqualifying
filter**, not merely a reported column — the opposite of the pairwise
pathway's `--skip_tblastn_filter` default. The family *is* its own validation
cluster: the family representative + membership feed `VALIDATE` directly, no
separate candidate-clustering step is needed (unlike the pairwise pathway,
which clusters candidates post hoc purely to deduplicate before TBLASTN).

### Known trade-off

Over-merging two true families reads as a false "present" (masks a real
novelty); under-splitting one true family into two reads as a false "absent"
(spurious novelty/loss call). The `mmseqs` identity threshold (`0.3`) is the
key knob controlling this and is meant to be swept with biological controls
(BUSCO single-copy recovery as the primary curation-free quality metric),
not tuned per-run.

---

## 3. `novelty_discovery` — two-phase targeted screen

**Where:** `workflows/novelty_discovery.nf`, `workflows/novelty_screen.nf`,
`lib/family_presence.py`, `bin/calibrate_family_hmms.py`; see CLAUDE.md
"Two-Phase Targeted Novelty Pipeline" for the config-format and category
semantics this section assumes.

This pathway answers a narrower question than the other two: "is this gene
specific to my *small* target lineage, and if not, how far does its
conservation actually extend?" It trades the O(N²)/O(N) search-cost
trade-off of the other pathways for a two-phase design that scales the
*negative control* instead of the whole outgroup, then only pays for a
broader search once a phase-1 candidate already looks interesting.

### Phase 1 — discovery, calibrated against a small negative-control panel

1. `DISCOVERY_TARGET` proteomes are clustered into families (mmseqs2) and
   family HMMs are built exactly as in the `mmseqs` pathway (same
   `BUILD_FAMILY_PROFILES` module, reused).
2. Each family HMM is searched (`hmmsearch`) against `DISCOVERY_TARGET` +
   `DISCOVERY_OUT` proteomes only — a small reference panel (typically 5-6
   proteomes), not the full outgroup.
3. **Per-family threshold calibration** (`bin/calibrate_family_hmms.py`):
   rather than one fixed E-value for every family, each family's threshold
   is set from its own best (lowest-E) hit against the `DISCOVERY_OUT`
   negative-control panel. If a family hits any `DISCOVERY_OUT` proteome,
   its threshold is tightened to just below that hit's E-value
   (`calibration_source = negative_control`); if it never hits the panel, it
   falls back to the global default (`--hmm_presence_evalue`, `1e-3`,
   `calibration_source = global_default`). This is the same self-calibration
   idea as the pairwise pathway's paralog cutoff, applied at family
   granularity against an explicit negative control instead of a paralog.
4. Presence is then called the same way as the `mmseqs` pathway —
   `full-sequence E < per-family threshold AND profile coverage >=
   hmm_presence_cov` (shared logic in `lib/family_presence.py`, used by both
   phase 1 and phase 2 so "present" means the same thing throughout this
   pathway).
5. Candidate rule: family present in `>= ingroup_min_frac` of
   `DISCOVERY_TARGET` and absent from every `DISCOVERY_OUT` proteome —
   structurally identical to the other pathways' novelty rule, just against
   a smaller reference set.

**Known gap:** phase 1 currently only evaluates multi-member families.
Singleton target proteins have no family HMM to search with, so a true
single-copy target-specific gene is invisible to this pathway today (see
CLAUDE.md for the tracked gap and the unimplemented single-target-genome
pairwise fallback).

### Phase 2 — screen against the broader clade

The same calibrated family HMMs from phase 1 (no rebuilding) are re-searched
against `NEAR_INGROUP` (close relatives) and `BROAD_OUTGROUP` (distant
lineages), using the identical presence rule. Every phase-1 candidate is
reclassified:

| `novelty_category` | Condition |
|---|---|
| `target_specific` | No hit in `NEAR_INGROUP` or `BROAD_OUTGROUP` |
| `clade_specific`  | Hit in `NEAR_INGROUP`, no hit in `BROAD_OUTGROUP` |
| `false_novelty`   | Hit in `BROAD_OUTGROUP` — dropped from the candidate list (label retained in the matrix for visibility) |

Only `target_specific`/`clade_specific` survive to annotation (Pfam +
SwissProt) — `false_novelty` is filtered out before the expensive annotation
step, not after. TBLASTN against `BROAD_OUTGROUP` genomes runs for reporting
only (published as `screen_tblastn_summary.tsv`; not yet surfaced in the
interactive report). A config with no `NEAR_INGROUP`/`BROAD_OUTGROUP` rows
degrades gracefully: zero proteomes to search means every phase-1 candidate
defaults to `target_specific`, since there's no broader evidence available
to demote it with.

### What's structurally different here vs. the other two pathways

- The **negative control is explicit and small** (`DISCOVERY_OUT`), calibrated
  per family, rather than "every proteome not in the query group."
- **Calibration happens once, reused twice** — phase 2 does not recompute
  thresholds; it reapplies phase 1's per-family E-value ceiling against a
  larger proteome set.
- There is **no loss direction** for this pathway (deferred future
  extension) — `main.nf` stubs a valid empty `losses.html` purely so the
  landing-page report can still assemble.

---

## Cross-pathway summary

| | `pairwise` | `mmseqs` | `novelty_discovery` |
|---|---|---|---|
| Search cost | O(\|IN\|×(N−1)) pairwise | O(N) hmmsearch (1 per proteome) | O(\|DISCOVERY_OUT\|) then O(\|NEAR_INGROUP\|+\|BROAD_OUTGROUP\|) |
| Presence unit | per-protein pairwise hit | per-family HMM hit | per-family HMM hit |
| Calibration signal | within-proteome paralog (self-search) | fixed E-value + coverage (swept, not per-family) | per-family negative-control panel (`DISCOVERY_OUT`) |
| Absence-call risk | paralog-competition edge cases (mitigated by `--paralog-competition-scope`) | family over/under-merging | missed singleton target genes (known gap) |
| TBLASTN role | reporting only (`--skip_tblastn_filter`) | disqualifying filter | reporting only (`BROAD_OUTGROUP`, not yet in report) |
| Candidate-cluster step | separate mmseqs pass over candidates only | none — family *is* the cluster | none — family *is* the cluster |
| Loss direction | yes (`LOSS_SEARCH` mirror) | yes (`PROFILE_LOSS_SEARCH` mirror) | not implemented (stubbed empty) |

All three converge on the same `presence_matrix.tsv` / `candidates.txt` /
cluster-TSV contract, which is what lets `workflows/validate.nf`,
`workflows/annotate.nf`, `workflows/summarize.nf` and the three HTML report
generators stay pathway-agnostic.
