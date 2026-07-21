# Family-profile search pathway as an alternative presence-matrix producer

## Status

Proposed — 2026-07-20; **design finalized by grilling 2026-07-20** (8 decisions
resolved, see "Resolved design" below). No implementation yet. Supersedes nothing;
runs alongside the existing pairwise SEARCH pathway.

## Context

The pairwise pathway (`workflows/search.nf`) queries every ingroup proteome
against every other proteome. On an order-scale config like
`configs/Chaetothyriales.csv` (81 ingroup + 128 outgroup = 209 proteomes) this is
`|IN| × (N−1) ≈ 16,848` pairwise jobs for the novelty direction and another
~26,600 for the loss direction, each searching ~11k query proteins against a target
proteome. Diamond tolerates this; phmmer does not. The combinatorial cost is the
blocker for scaling to large clades.

The redundancy the user wanted to collapse ("reduce the ingroup to a highly
identical representative set, search that, expand back") is the *same operation* as
gene-family clustering: once the ingroup is clustered, cluster membership already
encodes ingroup **presence** per species (exact, no search), the representative set
is the reduced probe, and the family is the natural unit for novelty, loss, and gene
contraction alike.

`modules/hmmbuild.nf` and `modules/hmmsearch.nf` are already imported but unwired —
their intended use is exactly this pathway.

## Decision

Add a second pathway, `PROFILE_SEARCH`, selected by a new `--cluster_tool` param
(`pairwise` = current default, `mmseqs` = this pathway, `orthofinder` = future).
It produces the **same three pivot artifacts** the rest of the pipeline already
consumes — `presence_matrix.tsv`, `candidates.txt`, and the cluster TSV — so
ANNOTATE / VALIDATE / SUMMARIZE / REPORT and the HTML/JSON reports are reused
unchanged.

Pathway:

```
mmseqs cluster (all ingroup proteomes, -s7, 0.3/0.8) → gene families + membership (= seed)
  → keep families with ≥2 members (singletons dropped)
  → famsa per family                        → family MSA
  → hmmbuild                                → family HMM database
  → hmmsearch family-HMM-db vs ALL 209 proteomes (both groups)  (~209 jobs)
  → uniform profile-based presence per family per species
  → expand family verdict to per-member-protein rows    → presence_matrix.tsv (same shape)
```

Presence (both groups, uniform) = per-seq `E < 1e-3` at fixed `--Z` AND profile-coverage
≥ 50% (from `--domtblout`). Novelty candidate = family present in ≥ `ingroup_min_frac` of
ingroup species AND zero outgroup presence; then the existing TBLASTN genomic validation
applies as a **filter** (not merely annotate — see Trade-offs). The family *is* the
validation cluster (its rep + membership feed `VALIDATE` directly). Loss direction mirrors
this with outgroup-seeded families.

## Rationale

- **Scale**: ~16,848 pairwise jobs → ~128 hmmsearch jobs (one per outgroup proteome,
  each scanning the whole family-HMM db in a single run). Ingroup presence needs no
  search at all — it is read from cluster membership.
- **Sensitivity where it matters**: novelty is a claim of *absence* from the
  outgroup. A single representative sequence can miss a diverged outgroup homolog
  that another family member would have detected → false novelty. A family **profile
  HMM** is the sensitive probe for the absence call, and is strictly better than the
  single-rep diamond alternative.
- **Reuse**: wires the dormant hmmbuild/hmmsearch modules; the family MSA it produces
  is also the substrate the gene-contraction analysis needs (ADR forthcoming).
- **Same pivot, free UI**: because the pathway emits the existing matrix/candidates/
  cluster-TSV contract, `lib/report_data.py` and the three report templates render
  without change. A `support` column (which pathway called each candidate) is the
  only report-side addition planned.

## Trade-offs

- **Over/under-merging.** Family boundaries are set by the mmseqs identity/coverage
  threshold. Over-merging two real families → false "present"; under-splitting one
  family → false "absent"/novelty. This is the same risk OrthoFinder carries and is
  the fundamental cost of a clustering-based (vs pairwise) absence call. The identity
  threshold is the key knob and must be swept (robust-analysis
  `sensitivity-analysis.md`).
- **TBLASTN must filter, not just report, in this pathway.** Clustering-based absence
  is more artifact-prone than the paralog-calibrated pairwise absence, so genomic
  validation should disqualify candidates here (unlike the current novelty direction,
  which runs with `--skip_tblastn_filter` and only annotates hits).
- **Two notions of "cluster".** `CONTEXT.md` defines *Cluster* narrowly as a group of
  **candidate** proteins for TBLASTN validation. This pathway clusters the **whole
  ingroup** up front to define families before any candidate exists. These are
  different units and need distinct vocabulary — introduce a **gene family** /
  **family profile** glossary entry (via `/domain-modeling`) rather than overloading
  *Cluster*.
- **Ingroup presence from membership vs search.** Reading ingroup presence from
  cluster membership (not a paralog-calibrated search) means the two pathways define
  "presence" slightly differently. That is acceptable — and the concordance between
  them is itself a validation signal — but it must be documented so results are not
  compared naively.

## Alternatives considered

- **Single representative + diamond probe** (instead of family HMM) — simpler, but
  less sensitive for the absence call; rejected in favour of profiles.
- **OrthoFinder first** — tree-aware, yields novelty + loss + a species tree
  (phylostrata) from one symmetric run, collapsing `LOSS_SEARCH` into a matrix
  re-filter. Deferred to a later phase as a cross-validation pathway; heavier, and its
  absence calls depend on orthogroup boundaries. Both pathways will coexist behind
  `--cluster_tool`.
- **Reduce the target (outgroup) side too** — novelty only needs "absent from the
  outgroup *union*", so outgroups could be concatenated into one DB. Rejected for now:
  per-outgroup-species resolution is needed by the presence matrix and reports;
  hmmsearch per outgroup proteome is already cheap enough.

## Phasing

1. `PROFILE_SEARCH` (mmseqs cluster → famsa → hmmbuild → hmmsearch-all-209), same matrix
   contract; validate matrix shape on a small config (`pezizo5`) before scaling.
2. Validation battery + parameter sweep (see Q8) and cross-method `support` column +
   concordance surfaced in the reports.
3. Gene-contraction analysis on core families (own ADR) — length + Pfam-architecture
   stats over the family MSAs, with the confounder controls in that ADR.
4. `--cluster_tool orthofinder` pathway (both directions from one run) + phylostrata
   via a species tree (PHYling / `nf_phyling`).

## Resolved design (grilling 2026-07-20)

Eight decisions, full rationale in `.living/decisions.md` ("ADR-0002 grilling resolutions"):

1. **Presence semantics** — the family HMM probes **both** groups; presence everywhere is
   one uniform, alignment-based criterion, cluster membership is only the seed (~209
   hmmsearch jobs). Avoids clustering silently deflating `ingroup_frac`.
2. **Clustering scope** — ingroup-seeded families + profile-search-all; mirror
   outgroup-seeded for loss. (Joint all-proteome clustering stays the deferred Phase-4
   OrthoFinder-lite method, preserving method independence for the Phase-2 concordance check.)
3. **Clustering** — `mmseqs cluster` (cascaded, `-s 7`) default `--min-seq-id 0.3 -c 0.8
   --cov-mode 0 --cluster-mode 0`; identity threshold **swept** (0.2–0.5) with biological
   controls (see 8).
4. **Profile construction** — add **famsa** to pixi + `hmmbuild` (chosen over mafft for
   scale; over mmseqs-MSA/native-profiles for HMMER interop + E-value model).
5. **Presence rule** — `hmmsearch --domtblout`, fixed `--Z`, presence = per-seq `E < 1e-3`
   AND profile-coverage ≥ 50%; E-value + coverage swept; lenient-for-absence + coverage
   guard against promiscuous shared domains; TBLASTN backstop. **Future:** per-family GA/TC/NC
   gathering thresholds from seed-score + decoy distributions (`--cut_ga`); HMMER3 needs no
   `hmmcalibrate`.
6. **Singletons** — profile only ≥2-member families (min a param); singletons dropped
   (shattered-family members recovered by family profiles' all-proteome search; true orphans
   fail the clade-frac filter); optional separate orphan-gene side table.
7. **Validation clustering** — family-as-cluster: family rep + membership feed `VALIDATE`
   directly, skipping candidate re-clustering (ADR-0001 stays for the pairwise pathway).
   TBLASTN-rep now; family-HMM-vs-translated-genome a future enhancement.
8. **Validation battery** — BUSCO `fungi_odb12` single-copy recovery as the primary
   curation-free family-quality metric; negatives = BUSCO/core must-not-be-novel; positives =
   silver standard (robust pairwise novelties) **plus optional hand-curated set**
   (`configs/controls/<clade>.controls.csv`, template scaffolded); sweep grid scored on all,
   pick the knee = shipped default.

## Remaining open questions

- Exact fixed `--Z` value (summed target size vs a constant like 1e7) — settle at
  implementation; both give cross-target comparability.
- Contraction null model and phylogenetic non-independence (deferred to the Phase-3 ADR).
