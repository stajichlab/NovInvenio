# Novelty Discovery & Screen — Two-Phase Targeted Novelty Plan

| Field | Value |
|-------|-------|
| **Date** | 2026-07-26 |
| **Author** | Jason Stajich |
| **Priority** | high |
| **Status** | in-progress |
| **Category** | feature |
| **Related ADRs** | ADR-0002 (family-profile search pathway) |
| **Related todos** | `profile-search-pathway.md` (Phase 1), `cross-method-support-column.md` (Phase 2) |
| **GitHub issues** | #24, #25, #26, #27, #28, #29 |

## Overview

A two-phase, target-focused novelty discovery pipeline that identifies
lineage-specific genes through progressive phylogenetic screening:

1. **`novelty_discovery`** — Start with a small target set (2-3 genomes) and a
   modest outgroup panel (5-6 species). Cluster the target proteomes into gene
   families using mmseqs2, build family HMMs (famsa + hmmbuild) for ALL
   multi-member families, and use hmmsearch against the outgroup panel to
   identify in-group specific families. Calibrate family HMM thresholds using
   the outgroup panel as a negative control.

2. **`novelty_screen`** — Search the calibrated family HMMs against a broader
   set of close relatives (near-ingroup, 20-50 species) and distant lineages
   (broad-outgroup, 20-100 species). Use the results to classify each family
   into one of three novelty categories. TBLASTN validates absence at the
   genomic level.

The plan addresses the family HMM reuse concern by implementing `storeDir`
caching on family HMM building steps, keyed by a hash of the family's sorted
member IDs. If the same family appears in a subsequent run (same target, or
overlapping target), the HMM is a cache hit. Adding a new target species that
doesn't change family X's membership → family X's HMM is a cache hit.

## Config Format

Extend the `GROUP` column in the analysis CSV with four new values:

| GROUP value | Role | Phase |
|-------------|------|-------|
| `DISCOVERY_TARGET` | The genome(s) we're looking for novelties in (2-3) | Discovery |
| `DISCOVERY_OUT` | Small reference panel for the absence call (5-6) | Discovery |
| `NEAR_INGROUP` | Close relatives, part of the same clade (20-50) | Screen |
| `BROAD_OUTGROUP` | Distant lineages, outside the target clade (20-100) | Screen |

Example:

```csv
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
DISCOVERY_TARGET,Neurospora crassa,OR74A,Ncrassa.pep.fa,Ncrassa.dna.fa,Ncra,Pezizomycotina
DISCOVERY_TARGET,Podospora anserina,,Pans.pep.fa,Pans.dna.fa,Pans,Pezizomycotina
DISCOVERY_OUT,Aspergillus fumigatus,Af293,Afum.pep.fa,Afum.dna.fa,Afum,Pezizomycotina
DISCOVERY_OUT,Saccharomyces cerevisiae,S288C,Scer.pep.fa,Scer.dna.fa,Scer,Saccharomycotina
NEAR_INGROUP,Coccidioides immitis,WA211,Cocci_WA211.pep.fa,Cocci_WA211.dna.fa,Cimm,Pezizomycotina
BROAD_OUTGROUP,Schizosaccharomyces pombe,,Spom.pep.fa,Spom.dna.fa,Spom,Taphrinomycotina
BROAD_OUTGROUP,Candida albicans,,Calb.pep.fa,Calb.dna.fa,Calb,Saccharomycotina
```

The existing config parser (`parse_config()`) will be extended to recognize
the new group labels. The four roles map cleanly to the two phases:

- `novelty_discovery` uses `DISCOVERY_TARGET` (ingroup) + `DISCOVERY_OUT` (outgroup)
- `novelty_screen` uses `NEAR_INGROUP` + `BROAD_OUTGROUP`

## Relationship to Existing Pipeline

The new two-phase plan sits **alongside** the existing `SEARCH`/`PROFILE_SEARCH`
workflows as a new `--cluster_tool` option:

- `--cluster_tool pairwise` (default) — existing O(N²) pairwise search pathway
- `--cluster_tool mmseqs` — existing family-profile search pathway (ADR-0002)
- `--cluster_tool novelty_discovery` — **new** two-phase targeted novelty plan

The existing pathways remain for backward compatibility and cross-method
validation (ADR-0002 Phase 2). The new `novelty_discovery` tool value selects
the two-phase plan.

## Loss Direction

Loss analysis (present in outgroup, absent from ingroup) is a **future
extension**. The two-phase architecture is designed to be mirrorable for
losses (just swap target/outgroup roles), but adding loss now would double the
scope. Keep the plan focused on novelty.

## Caching Strategy

Implement `storeDir` caching on family HMM building steps, keyed by a hash of
the family's sorted member IDs:

- If the same family appears in a subsequent run (same target, or overlapping
  target), the HMM is a cache hit.
- Adding a new target species that doesn't change family X's membership →
  family X's HMM is a cache hit.
- The `storeDir` approach is more granular than a versioned HMM database
  (option b, rejected): it enables per-family reuse across runs with
  overlapping ingroups, not just identical ingroups.
- `-resume` alone (option c, rejected) is fragile and doesn't work across
  separate invocations or when the work directory is cleaned.

This is the most granular and reusable approach.

## novelty_discovery Subworkflow

### Inputs

- `DISCOVERY_TARGET` proteomes (2-3 genomes) — the ingroup
- `DISCOVERY_OUT` proteomes (5-6 species) — small reference panel for absence call
- `DISCOVERY_OUT` genomes (5-6 species) — for TBLASTN genomic validation

### Branching Logic

The subworkflow branches based on the number of `DISCOVERY_TARGET` genomes:

- **1 target genome** → pairwise search (phmmer/diamond/blast) of the single
  ingroup proteome against each `DISCOVERY_OUT` proteome. No clustering, no HMMs.
  This is the existing `SEARCH` pathway with `|IN| = 1`.
- **≥2 target genomes** → mmseqs cluster + family HMM approach (see below).

### Steps (≥2 target genomes)

1. **mmseqs cluster** — Cluster the `DISCOVERY_TARGET` proteomes into gene families
   using mmseqs2 (`-s 7`, `--min-seq-id 0.3 -c 0.8 --cov-mode 0`). This is
   the existing `MMSEQS_FAMILY_CLUSTER` process (ADR-0002 Q2/Q3).
   - Option A: Cluster ingroup only, scan outgroups with family HMMs
     (**selected**).
   - Option B: Cluster everything together (target + disc_out), filter
     post-hoc for target-only families (**divergent path, test later**).

2. **Extract family sequences** — Split the mmseqs cluster into per-family
   member FASTAs (≥2 members per family, ≤100 max members). This is the
   existing `EXTRACT_FAMILY_SEQS` process.

3. **Build family HMMs** — For ALL multi-member families, run famsa (MSA) +
   hmmbuild (HMM) per family. This is the existing `BUILD_FAMILY_PROFILES`
   scatter-gather (ADR-0002 Q4/Q6, issue #14).
   - **Caching:** `storeDir` keyed by hash of family's sorted member IDs.
   - **Singletons:** Use pairwise search (phmmer) of singletons vs `DISCOVERY_OUT`
     proteomes (hybrid approach, option c). Merge results into one presence
     matrix.

4. **hmmsearch against outgroups** — Search all family HMMs against the
   `DISCOVERY_OUT` proteomes using hmmsearch (`--domtblout`, `E < 1e-3` AND
   `coverage ≥ 50%`). This is the existing `FAMILY_HMMSEARCH` scatter-gather
   (ADR-0002 Q1/Q5, issue #22).

5. **Filter to in-group specific families** — Keep families with no hmmsearch
   hit in any `DISCOVERY_OUT` proteome. These are the in-group specific families.

6. **TBLASTN genomic validation** — Search in-group specific family
   representatives against `DISCOVERY_OUT` genomes. Validates absence at the
   genomic level, catching genes missing from proteome annotations. This is
   the existing `TBLASTN` process.

7. **Calibrate family HMMs** — Set per-family thresholds using the `DISCOVERY_OUT`
   panel as a negative control (option c):
   - If the HMM hits any `DISCOVERY_OUT` proteome, set the threshold above the
     highest-scoring outgroup hit.
   - If the HMM doesn't hit any `DISCOVERY_OUT` proteome, fall back to the global
     E-value threshold (default 1e-3).
   - This directly calibrates against the false-positive scenario (calling a
     family present when its distant homolog exists in the outgroup).

### Outputs

- **Calibrated family HMMs** — The in-group specific family HMMs with
  per-family thresholds set. These carry forward to `novelty_screen`.
- **Discovery presence matrix** — protein × proteome 0/1 matrix for the
  discovery phase.
- **In-group specific candidate list** — `source_proteome::protein_id` lines.
- **TBLASTN summary** — protein × outgroup-genome 0/1 matrix.

## novelty_screen Subworkflow

### Inputs

- Calibrated family HMMs from `novelty_discovery`
- `NEAR_INGROUP` proteomes (20-50 species) — close relatives, part of same clade
- `BROAD_OUTGROUP` proteomes (20-100 species) — distant lineages, outside target clade
- `BROAD_OUTGROUP` genomes (20-100 species) — for TBLASTN genomic validation

### Steps

1. **hmmsearch against near-ingroup** — Search calibrated family HMMs against
   `NEAR_INGROUP` proteomes using hmmsearch (`--domtblout`, per-family thresholds
   from calibration, `E < 1e-3` AND `coverage ≥ 50%`). One-way HMM searches
   should be fast.

2. **hmmsearch against broad-outgroup** — Search calibrated family HMMs against
   `BROAD_OUTGROUP` proteomes using hmmsearch (same parameters). Catches false
   novelties where a divergent outgroup homolog exists.

3. **TBLASTN genomic validation** — Search family representatives against
   `BROAD_OUTGROUP` genomes. Validates absence at the genomic level. This is the
   existing `TBLASTN` process, applied to the broad-outgroup genomes.

4. **Three-category refinement** — Classify each family into one of three
   novelty categories based on hmmsearch results:
   - **Target-specific novelty:** Not found in `NEAR_INGROUP` or `BROAD_OUTGROUP`. Unique
     to the target genome(s). Strongest novelty candidate.
   - **Clade-specific novelty:** Found in `NEAR_INGROUP`, not found in `BROAD_OUTGROUP`.
     Shared with close relatives, but absent from distant lineages. Still a
     valid novelty claim for the clade.
   - **False novelty:** Found in `BROAD_OUTGROUP`. The family is widespread, not
     clade-specific. Remove from novelty candidate list.

5. **TBLASTN filter (optional)** — If a family has no TBLASTN hit in any
   `BROAD_OUTGROUP` genome, it's a confirmed absence at the genomic level. If
   `--skip_tblastn_filter` is passed, TBLASTN hits are reported but don't
   filter candidates.

### Outputs

- **Screened presence matrix** — protein × proteome 0/1 matrix for the screen
  phase, with `novelty_category` column added (target-specific / clade-specific
  / false novelty).
- **Three-category candidate list** — Refined novelty candidates, classified by
  category.
- **TBLASTN summary** — protein × broad-outgroup-genome 0/1 matrix.

## Post-Screen: Annotation & Reports

### Annotation

After `novelty_screen`, annotate only the final surviving novelty candidates
(option a). Same as the current pipeline's `ANNOTATE` step, but applied to the
three-category refined set. Annotation is expensive (hmmsearch vs Pfam-A,
diamond vs SwissProt), and there's no reason to annotate families that will be
filtered out as false novelties in the screen phase.

Annotation sources (priority order):
1. Model organism gene names (via `--modelorgs_config` YAML)
2. Pfam domain names (all unique domains per protein)
3. SwissProt description (best hit, `sp|ACCN|ID` prefix stripped)

### Reports

Both (option c):

1. **Reuse existing reports** — produce `novelties.html`, `core.html`,
   `losses.html` as the current pipeline does, but with a new
   `novelty_category` column (target-specific / clade-specific / false novelty).
   Minimal new report code.

2. **New two-phase report** (future enhancement) — a single report showing the
   discovery → screen → refinement pipeline, with the three categories
   visualized. More work but more informative for the new workflow.

## Future Enhancements

1. **Quantitative novelty scoring** — Instead of the three-category refinement,
   use the fraction of near-ingroup and broad-outgroup species with a hit to
   assign a continuous novelty confidence score. This would provide finer
   resolution than the binary present/absent call. (Noted as option c in the
   grilling.)

2. **Option B clustering** — Cluster everything together (target + disc_out),
   filter post-hoc for target-only families. This is a divergent path from the
   selected Option A (cluster ingroup only). Could be tested later to compare
   results. The key difference: Option B asks "does this family contain zero
   outgroup proteins?" (clustering-based membership), while Option A asks "is
   this family's HMM absent from the outgroup proteomes?" (sensitive profile
   search).

3. **Loss direction** — Mirror the two-phase architecture for losses (present
   in outgroup, absent from ingroup). The architecture is designed to be
   mirrorable: just swap target/outgroup roles. (Noted as future extension b
   in the grilling.)

4. **Per-family GA/TC/NC gathering thresholds** — Compute per-family score
   cutoffs from seed alignment scores and/or decoy distributions (ADR-0002 Q5
   future enhancement). More precise than the negative-control calibration
   (option c, selected): a fast-evolving family gets a looser threshold, a
   conserved family gets a tighter one. (Noted as option b in the grilling.)

5. **New two-phase report** — A single report showing the discovery → screen →
   refinement pipeline, with the three categories visualized. (Noted as future
   enhancement in the reports section.)

## Edge Cases

- **Empty target** — If no `DISCOVERY_TARGET` genomes are specified, error out.
- **No multi-member families** — If all target proteins are singletons (e.g.,
  highly diverged genomes), fall back to pairwise search for all proteins.
- **No in-group specific families** — If all families have hits in the
  `DISCOVERY_OUT` panel, the discovery phase produces an empty candidate list. The
  screen phase is skipped.
- **Single target genome** — Use pairwise search (phmmer/diamond/blast) of the
  single ingroup proteome against each `DISCOVERY_OUT` proteome. No clustering, no
  HMMs. This is the existing `SEARCH` pathway with `|IN| = 1`.

## Pipeline Architecture Diagram

```
novelty_discovery:
  DISCOVERY_TARGET proteomes (2-3) ─→ mmseqs cluster ─→ families
  │                                              │
  │                                              ├─≥2 members─→ famsa+hmmbuild ─→ family HMMs
  │                                              │                    │
  │                                              │                    └─→ hmmsearch vs DISCOVERY_OUT ─→ filter
  │                                              │                                                    │
  │                                              └─singletons─→ pairwise vs DISCOVERY_OUT ────────────────→ filter
  │                                                                                                    │
  DISCOVERY_OUT proteomes (5-6) ─────────────────────────────────────────────────────────────────→ hmmsearch targets
  DISCOVERY_OUT genomes (5-6) ───────────────────────────────────────────────────→ TBLASTN targets  │
  │                                                                                              │
  └─→ in-group specific families + calibrated HMMs + presence matrix + TBLASTN summary ←────────┘
       │
       └─→ novelty_screen:
             calibrated family HMMs
                   │
                   ├─→ hmmsearch vs NEAR_INGROUP (20-50) ──→ near-ingroup presence
                   │                                        │
                   ├─→ hmmsearch vs BROAD_OUTGROUP (20-100) ─→ broad-outgroup presence
                   │                                        │
                   └─→ TBLASTN vs BROAD_OUTGROUP genomes ──→ genomic validation
                                                            │
                   three-category refinement ←──────────────┘
                   │
                   ├─ target-specific novelty (not in NEAR_INGROUP or BROAD_OUTGROUP)
                   ├─ clade-specific novelty (in NEAR_INGROUP, not in BROAD_OUTGROUP)
                   └─ false novelty (in BROAD_OUTGROUP) → remove
                                                            │
                   annotation (Pfam, SwissProt, model orgs) ←┘
                                                            │
                   reports (novelties.html, core.html, losses.html) ←┘
                   with novelty_category column
```

## Implementation Notes

- The `novelty_discovery` subworkflow reuses the existing `MMSEQS_FAMILY_CLUSTER`,
  `EXTRACT_FAMILY_SEQS`, `BUILD_FAMILY_PROFILES`, `FAMILY_HMMSEARCH`, and
  `TBLASTN` processes.
- The `novelty_screen` subworkflow reuses the existing `FAMILY_HMMSEARCH` and
  `TBLASTN` processes.
- The config parser (`parse_config()`) needs to be extended to recognize the
  new `GROUP` values (`DISCOVERY_TARGET`, `DISCOVERY_OUT`, `NEAR_INGROUP`, `BROAD_OUTGROUP`).
- The `--cluster_tool novelty_discovery` option needs to be added to
  `nextflow.config` and `main.nf`.
- The `storeDir` caching on family HMM building steps needs to be implemented
  with a hash of the family's sorted member IDs as the cache key.
- The three-category refinement logic needs to be implemented in a new
  `bin/novelty_screen.py` script.
- The `novelty_category` column needs to be added to `lib/report_data.py`'s
  `ROW_FIELDS` and rendered as a filter + column in the existing reports.
