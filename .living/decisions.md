# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

## [2026-07-20] Mycelium init scoped to the knowledge layer only

**Context**: Initialized this repo as a mycelium living repository. NovInvenio is a
Nextflow DSL2 pipeline with its own established layout (`modules/`, `workflows/`,
`bin/`, `lib/`, and a real `data/` holding `cds/dna/pep` FASTA inputs), not a generic
data-analysis project.

**Decision**: Scaffold only the `.living/` knowledge layer, the three core convention
packs (robust-analysis, report-generator, idea-generator), `todo/`, `skillpacks/`, and
the Claude hooks. Skip the generic `analysis/`, `algorithms/`, `reference_material/`, and
`data/raw|processed|metadata` directories the default `init_repo.py` creates. Preserve the
existing 601-line `CLAUDE.md` (added a short "Living Repository" section rather than
generating from the template).

**Alternatives considered**:
- Full mycelium scaffold — rejected: would create empty `analysis/`/`algorithms/`/
  `reference_material/` dirs that clash semantically with the pipeline's own structure and
  with the existing `data/` (genuine FASTA inputs, not mycelium data buckets).
- `--restructure` mode — rejected: unimplemented in `init_repo.py` (audit-only, prints a
  TODO and exits without moving anything).

**Consequences**: `validate_structure.py` reports 9 "errors" — all of them the
deliberately-skipped generic dirs/files; the `.living/` layer itself validates clean.
Hooks were placed in `.claude/settings.json` (not `settings.local.json`) because the
harness owns `settings.local.json` and clobbered an out-of-band hooks write there;
`settings.json` hooks merge with the local file and survive permission updates.

**Tags**: mycelium, init, project-structure, nextflow, hooks, tooling

## [2026-07-20] Add a family-profile search pathway to scale beyond pairwise search

**Context**: Order-scale configs (`Chaetothyriales.csv`: 81 ingroup + 128 outgroup)
make the one-way pairwise SEARCH combinatorially infeasible (~16,848 novelty jobs,
~26,600 loss jobs), especially for phmmer. Need to (a) scale, (b) optionally add
OrthoFinder, and (c) add a gene-contraction analysis — all while reusing the existing
presence-matrix pivot and HTML/JSON reports.

**Decision**: Add a `PROFILE_SEARCH` pathway (`--cluster_tool mmseqs`): cluster the
ingroup once (membership = ingroup presence), build per-family HMMs (mafft + hmmbuild),
and probe each outgroup proteome with hmmsearch (~128 jobs). Emit the same
`presence_matrix`/`candidates`/`cluster_tsv` contract so ANNOTATE/VALIDATE/REPORT are
reused unchanged. Family HMM profile (not a single representative) is the sensitive
probe for the absence call. Full rationale, trade-offs, and 4-phase plan (profile
pathway → cross-method support column → gene contraction → OrthoFinder+phylostrata) in
`docs/adr/0002-family-profile-search-pathway.md`. Chose design-doc-first: no code yet.

**Alternatives considered**:
- Single-representative + diamond probe — rejected: less sensitive for the absence call.
- OrthoFinder as the first new pathway — deferred to phase 4 (heavier; boundary-dependent
  absence); will coexist behind `--cluster_tool` for cross-validation.
- Reduce the outgroup/target side too — rejected: per-species resolution needed by the
  reports; hmmsearch per outgroup is already cheap.

**Consequences**: Wires the dormant `hmmbuild`/`hmmsearch` modules. Introduces a
whole-ingroup "gene family" unit distinct from `CONTEXT.md`'s candidate-level *Cluster* —
needs a glossary entry (`/domain-modeling`). Clustering-based absence is more
artifact-prone than paralog-calibrated pairwise absence, so TBLASTN must *filter* (not
just annotate) in this pathway, and the family-definition mmseqs threshold must be swept.

### [2026-07-20] ADR-0002 grilling resolutions (running — folds into the ADR when complete)

Working through the gene-family decision tree one question at a time to finalize Phase 1
before coding. Resolved answers accumulate here; consolidated into
`docs/adr/0002-family-profile-search-pathway.md` "Open questions" once the grilling ends.

- **Q1 — presence semantics = (a)**: the family HMM probes **both** groups; presence
  (ingroup and outgroup alike) is one uniform, alignment-based criterion, and cluster
  membership is only the *seed*. Rejected (b) membership-for-ingroup + search-outgroup —
  the asymmetry would let clustering silently drop a diverged ingroup ortholog and deflate
  `ingroup_frac`. Cost: ~209 hmmsearch jobs instead of ~128 (negligible). Env fact
  surfaced: pixi has hmmer/diamond/mmseqs2 but **no aligner** — the "mafft per family" step
  needs a new dep or an alt route (mmseqs `result2msa`, or native mmseqs profiles); the
  existing `HMMSEARCH` module (one-HMM-vs-SwissProt) needs reworking for the family-scan
  orientation.
- **Q2 — clustering scope = (a)**: cluster the ingroup only to define/seed families, then
  hmmsearch those profiles against all 209 proteomes; mirror with outgroup-seeded families
  for the loss direction. Rejected (b) all-proteome joint clustering with membership-based
  presence — that is the deferred Phase-4 OrthoFinder-lite approach; folding it in now
  would discard the sensitive HMM absence call (Q1) and destroy method independence for the
  Phase-2 concordance check.
- **Q3 — clustering algorithm + params = accepted default**: `mmseqs cluster` (cascaded,
  `-s 7`; NOT linclust) with shipped default `--min-seq-id 0.3 -c 0.8 --cov-mode 0
  --cluster-mode 0`. Bidirectional 80% coverage is the guard against promiscuous-domain
  over-merging; greedy set-cover avoids single-linkage chaining. **User added**: the
  identity threshold must be *swept* (e.g. 0.2/0.3/0.4/0.5) and backed by **experimental
  confirmation tests** — biological positive controls (known lineage-specific genes that
  SHOULD be called novel) and negative controls (known core/conserved genes that must NOT),
  plus single-copy-ortholog (BUSCO) recovery per threshold. This folds forward into Q8
  (validation) and becomes a Phase-2 deliverable, not just a config knob.
- **Q4 — profile / MSA construction = (a) famsa + hmmbuild**: add `famsa` to pixi.toml
  `[workspace]` (chosen over mafft for scale-safety on large families; conda-available),
  align each family, `hmmbuild` → family HMM. Keeps the dormant HMMER modules and HMMER's
  well-understood E-value model for the Q5 absence call. Rejected mmseqs `result2msa`
  (lower MSA quality weakens the profile) and native mmseqs profiles (drops HMMER / weaker
  interop with the Pfam hmmscan already in ANNOTATE). Dep is added at implementation time
  via the tickets→branch→PR flow, not during design.
- **Q5 — presence rule from hmmsearch = accepted default**: `hmmsearch --domtblout`,
  **fixed `--Z`** (comparable E-values across all 209 targets of differing size), presence =
  per-sequence `E < 1e-3` **AND** profile-coverage ≥ 50% (summed envelope / HMM length from
  domtblout). Lenient E-value because false-*absence* manufactures false novelty (the
  expensive error); the coverage guard stops a single promiscuous shared domain from faking
  family-wide presence (profile-pathway analog of the pairwise paralog-competition filter).
  Both E-value and coverage are swept; TBLASTN-vs-genomes is the backstop.
  **User addition — HMM calibration / per-family gathering thresholds (GA):** HMMER3 needs
  no `hmmcalibrate` (E-values are analytic from hmmbuild's MU/LAMBDA), so no calibration
  step. BUT develop **per-family GA bit-score thresholds** as a Phase-2 refinement replacing
  the single global E-value: derive each family's trusted cutoff (TC) from the score
  distribution of its own seed members and a noise floor (NC) from searching a
  shuffled/reversed decoy DB, set GA between them, write `GA`/`TC`/`NC` lines into each
  `.hmm`, and call presence with `hmmsearch --cut_ga`. This adapts stringency per family
  (fast-evolving families get a lower bar, conserved families higher) — more principled than
  one global threshold. Ships after the global-default version; validated by the same
  sweep + biological controls.
- **Q6 — singletons / orphan families = accepted**: build profiles only for multi-member
  clusters (≥2 members; minimum is a param). Drop singletons from the profile pathway —
  genuine members of a clustering-shattered family are recovered by the multi-member family
  profiles' all-proteome search (Q1), and a true species-unique orphan sits in 1 species
  (`ingroup_frac ≈ 1/81`) so it fails the 0.75 clade-novelty filter anyway. Optionally emit
  singletons as a separate species-restricted / orphan-gene table (future category), not
  forced through the clade-novelty machinery. Only unrecoverable case (a family so diverged
  no two members cluster at 0.3/0.8) is rare at order level and minimized by `-s 7`.
- **Q7 — relation to ADR-0001 candidate clustering + TBLASTN = family-as-cluster**: in the
  profile pathway the family *is* the validation cluster — feed family representative +
  family membership (`rep⇥member`) straight into `VALIDATE`, skipping `EXTRACT_CANDIDATES` +
  `MMSEQS_CLUSTER`. Preserves and improves ADR-0001's rationale (cluster before TBLASTN);
  ADR-0001 stays in force for the pairwise pathway (no upstream families). `--cluster_tool`
  branch: pairwise → cluster candidates; profile → pass the family `cluster_tsv` through
  (matches ADR-0002's "CLUSTER accepts a pre-made cluster tsv" note). Genomic probe: reuse
  **TBLASTN of the family rep vs outgroup genomes** for Phase 1; **future enhancement** =
  family HMM vs 6-frame-translated outgroup genomes (translated-ORF hmmsearch / miniprot) as
  the most sensitive DNA-level absence backstop (deferred — scope).
- **Q8 — family-definition validation & experimental confirmation = accepted battery**:
  (1) **BUSCO `fungi_odb12` single-copy-ortholog recovery** as the primary curation-free
  family-quality metric (recovered as single family, ~1/species → detects over/under-split);
  reuses PHYling/nf_phyling. (2) **Negative controls** = BUSCO/core genes must NOT be called
  novel (specificity / false-novelty rate). (3) **Positive controls** = silver standard from
  robust pairwise-pathway novelties (recall via cross-method concordance on pezizo5 → scale)
  **PLUS an optional hand-curated set** (user wants this option). (4) **Sweep grid**
  `min-seq-id{0.2,0.3,0.4,0.5} × cov{0.5,0.7} × E{1e-3,1e-5}` scored on #families/#novelties/
  BUSCO-recovery/neg-FP-rate/pos-recall → pick the knee = shipped default. (5) TBLASTN removal
  count = annotation-artifact rate.
  **Hand-curated control-set format** (`configs/controls/<clade>.controls.csv`): columns
  `control_id, class{positive|negative}, expected_call{novel|core}, anchor_type{protein_id|
  fasta|busco}, anchor, proteome_short, gene_name, expected_origin, source, notes`. The
  validator resolves each anchor to a family (protein_id → containing family; fasta →
  hmmsearch/assign to nearest family; busco → BUSCO assignment) and checks the family's call
  == expected_call → per-set recall (positives) and FP rate (negatives). Positive *novelty*
  controls are typically anchored by protein_id/sequence (clade-restricted genes are often
  hypothetical/unnamed), negatives by BUSCO id or housekeeping protein_id.

**PHASE 1 BUILT (2026-07-20)** — issue #3, draft PR #4, branch `3-profile-search-pathway`.
Implemented the family-profile pathway per the resolved design: pixi +famsa; modules
mmseqs_family_cluster / build_family_profiles (famsa+hmmbuild) / family_hmmsearch
(`--domtblout`, fixed `-Z`, loose report -E) / profile_presence_matrix; workflows/
profile_search.nf; `main.nf` `--cluster_tool pairwise|mmseqs` branch (pairwise unchanged);
nextflow.config sweepable params; bin/extract_family_seqs.py + bin/profile_to_matrix.py.
Verified: pytest 70 passed (5 new = matrix/candidates contract parity, singleton drop,
coverage filter), `nextflow config` OK, `nextflow lint` 0 errors. **Scope note:** first cut
is novelty-direction only and REUSES the existing CLUSTER→VALIDATE chain (candidate
re-clustering) rather than Q7 family-as-cluster; the loss direction still runs pairwise
LOSS_SEARCH. Both, plus the Phase-2 sweep/BUSCO/controls battery, are documented PR
follow-ups. End-to-end pezizo5 parity is the PR review gate (needs pixi+famsa + data on SLURM).
Followed the tickets→branch→PR workflow (CLAUDE.md).

**PHASE 2 KICKED OFF 2026-07-21 — HANDOFF, PAUSED MID-RECON.** Chose "Phase 2 as planned"
(sweep + BUSCO/controls validation battery + cross-method support column,
`todo/cross-method-support-column.md`). `pixi install` picked up famsa successfully. Before
building further, attempted a smoke test of famsa on the current dev node (`c09`,
`amd,abu_dhabi` — an old Opteron "Abu Dhabi" generation) to close the "never run
end-to-end" gap from Phase 1.

**Finding (important, unresolved): `famsa` SIGILLs (illegal instruction, core dump) on
`c09`.** `/proc/cpuinfo` confirms `c09` has AVX but NOT AVX2 — the bioconda famsa build
requires AVX2+. Confirmed via `sinfo -o "%N %P %f"`: the `short` SLURM partition (the
pipeline's default queue, `nextflow.config` line ~120) mixes `c[01-30]` (`amd,abu_dhabi`,
NO AVX2) with `i[01-62]` (`intel,broadwell`, HAS AVX2) — so a real SLURM run of
`BUILD_FAMILY_PROFILES` (the famsa+hmmbuild process) could land on either and silently
crash depending on scheduler placement. The **`epyc` partition is `ryzen,amd,milan`
(Zen2)** — AVX2/AVX512-capable — confirmed via `sinfo -h -o "%P"` / feature query,
consistent with `MMSEQS_CLUSTER`'s existing `clusterOptions = '-C ryzen'` precedent
(`nextflow.config` line ~163) for the same class of problem.

**User decision: defer the fix — "save for handoff later, we will run this on an avx2
node."** Not fixed yet. On resume: add a SLURM `-C` constraint (ryzen/broadwell/cascade,
excluding `abu_dhabi`) to `BUILD_FAMILY_PROFILES` in `nextflow.config`'s slurm profile
(mirror `MMSEQS_CLUSTER`'s pattern), OR pin the run to an explicit AVX2 partition/node via
`-C` or `--queue`. Local/interactive smoke-testing of the family-profile pathway is NOT
possible on `c09` — needs an AVX2+ allocation (ryzen/epyc, or i-nodes/broadwell, or x-nodes
/cascade). No code changes made this sub-session; `pixi.toml`'s famsa dependency itself is
correct and unchanged — the issue is SLURM node placement, not the dependency choice.

**PHASE 1 MERGED 2026-07-21** — PR #4 merged to main (issue #3 closed), local branch
`3-profile-search-pathway` deleted post-merge. `todo/profile-search-pathway.md` and the
registry flipped to `complete`. Still true: end-to-end `--cluster_tool mmseqs` run on
pezizo5 (shape-parity + reports render) was never executed — only unit-tested — so it
remains an open verification item, naturally folded into Phase 2's validation battery.

**GRILL COMPLETE — all 8 resolved AND consolidated.** Done: folded into `docs/adr/0002`
("Resolved design" section + updated pathway diagram/phasing); Phase-1
(`profile-search-pathway.md`) and Phase-2 (`cross-method-support-column.md`) todos updated
with settled params + validation battery; scaffolded `configs/controls/`
(Chaetothyriales.controls.csv template + README); `mycelium-stop-check.sh` Stop hook
**RESTORED** in `.claude/settings.json`; `last-session.md` final write done. ADR-0002 design
is FINAL — ready to implement Phase 1 via the tickets→branch→PR workflow.

**Session note:** blocking Stop hook `mycelium-stop-check.sh` temporarily removed from
`.claude/settings.json` for this grill (user-authorized); `mycelium-data-lineage-stop.sh`
(non-blocking) retained. **RESTORE the stop-check hook at grill end**, then do the final
session-end write.

**Tags**: adr-0002, gene-family, grilling, mmseqs, hmmer, in-progress

**Tags**: architecture, scaling, nextflow, mmseqs, hmmer, orthofinder, novelty, loss, contraction, adr-0002

## [2026-08-06] Reconciled the 2026-07-22 → 08-01 investigation window into the living layer

**Context**: A reconciliation pass (issue: "living repo initialize / repopulate") found that
none of the substantial analysis and engineering work 2026-07-22 → 2026-08-01 had been
disseminated into `.living/decisions.md` / `learnings.md` / `findings/`. The session logs
(`.living/log/`) were mostly bare "Session started" stubs (except 07-24), the
`LOG_REGISTRY.md` had a single row, and `INDEX.md` had not been regenerated since 07-21.
Root cause of the silent INDEX staleness: the mycelium health hook invokes **bare `python3`
(Python 3.9 on this cluster)** to run `generate_index.py`, which uses `str | None` union
syntax (3.10+), so it crashes and is swallowed by the hook's `|| true`. Regenerating with
`python3.12` works. The substance below was reconstructed from `git log --all` (per-branch,
per-commit), the run logs in `logs/`, and `results/*` counts.

**Decision / actions**: Reconcile the window now so future sessions see the full picture.
Capture the engineering decisions and the quantitative investigation findings here and in
`learnings.md`, backfill `LOG_REGISTRY.md`, and rebuild `INDEX.md`.

**Findings — per-thread (07-22 → 08-01)**:

1. **Novelty-discovery two-phase pipeline (issues #24–#29, #31–#42; merged)** — third
   `--cluster_tool novelty_discovery` pathway; four config GROUP roles
   (`DISCOVERY_TARGET`/`DISCOVERY_OUT`/`NEAR_INGROUP`/`BROAD_OUTGROUP`, with aliases
   `TARGET`/`DISC_OUT`/`NEAR_IN`/`BROAD_OUT` normalized on both the Python (`GROUP_ALIASES`)
   and Groovy (`normalizeGroup()`) sides). `NOVELTY_DISCOVERY` (phase 1: family cluster →
   HMM → hmmsearch → singleton extract → negative-control calibration → TBLASTN) and
   `NOVELTY_SCREEN` (phase 2: three-category re-classification
   `target_specific`/`clade_specific`/`false_novelty`; `false_novelty` removed before
   expensive ANNOTATE). Also surfacing `novelty_category` in the report.
2. **#33 wiring class-bug** — the discovery workflow "ran to completion but never rendered":
   `Channel.empty().collect()` silently produces nothing and starves every downstream
   process; made value channels correctly. Also TBLASTN was queried with the HMM DB instead
   of family-representative protein FASTA; `make_novelties.py`/`report_data.py` only
   recognized strict IN/OUT.
3. **Singleton screening (PR #52/#53; merged)** — closed two gaps: (a) singleton search
   lacked a paralog-competition guard (`NCU08332`/HEX-1-vs-eIF5A false positive); fixed by
   self-vs-self paralog calibration + extended singleton query; (b) phase-2 screen had no
   singleton evidence, so singletons always defaulted `target_specific`; re-search fixed it.
   Bonus: `extract_singletons.py` mis-counted reps (`rep\trep` self-line) as singletons.
4. **Context search (PR #48–#51; merged)** — NEAR_INGROUP/BROAD_OUTGROUP as **report-only**
   context for `--cluster_tool pairwise` (candidate + paralog searched, never affects strict
   IN/OUT novelty filtering). Reverted #47's empty-placeholder matrix columns that shadowed
   it. `stageAs` fix for stub-file input filename collisions.
5. **Roles fix (PR #47; merged)** — presence-matrix columns + control scoring moved from
   strict `group == 'IN'/'OUT'` to coarse `INGROUP_ROLES`/`OUTGROUP_ROLES` banding.
6. **Hit e-values (PR #44/#46; merged)** — `build_presence_matrix.py` now emits an e-value
   **sidecar** (`--output-evalues`, phase-1 novelty-discovery only) surfaced as `'ev'` in the
   report detail panel. **Gap:** mmseqs/PROFILE_SEARCH still lacks e-values (stubbed
   `EMPTY_EVALUES_STUB`).
7. **Containerization + HPCC (PR #41–#43, #45; merged)** — Dockerfile + `novinvenio.def`
   singularity def, container/directive on all 38 processes, refactored `nextflow.config`
   process.container, OIDC keyless cosign build+sign workflow, `conf/ucr_hpcc_slurm.config`.
   Simpler/faster local runs no longer activate pixi `beforeScript` per process.
8. **Pfam annotation speedup (07-31 → 08-01; committed direct to `main`)** —
   `hmmsearch`→`hmmscan` orientation swap (large DB on the parallelized target side; ~2.8→real
   threading), then moved Pfam out of `ANNOTATE_MATRIX` into a chunked scatter-gather
   `ANNOTATE_PFAM` (`SPLIT_CANDIDATES` → `PFAM_HMMSCAN_CHUNK` → `MERGE_PFAM_TBLOUT`), with
   `scratch = true` (2.2× wall-clock from staging pressed DB on node-local disk), and
   `SLURM_CPU_BIND=none` for MPI in `modules/hmmsearch.nf`.

**Findings — run/investigation evidence (`logs/`, `results/*`)**: see `learnings.md`
"investigation-run observations" for the full list; the headline numbers were:
- **mmseqs family-profile yields FEWER novelty candidates than pairwise on identical data**
  (pezizo5: 1851 vs 2544). Form/identity presence-calling diverges — cross-method concordance
  (Phase 2 todo) is warranted.
- **Singleton search roughly doubles recovered target genes** (PR #52): sordario per-species
  novelties ~422→907 (Afum), 329→693 (Cimm), 419→1173 (FgraPH1), 429→956 (Ncra);
  `target_specific` 1,599→3,729. The singleton pathway has a large, previously-invisible
  payoff.
- **Antarctolithicomycotina is novelty-rich** and shows an unusual **loss_candidates (1,697)
  > novelty candidates (1,530)** from just 2 ingroup genomes.
- **`preempt` did NOT rescue the mmseqs BUILD_CHUNK work** — the 07-24 handoff ("wait for
  SLURM job 26715986") outcome was `completed=5063 failed=45 cached=57` from famsa 600s
  timeouts (exit 140), 100 ABORTED + 28 FAILED BUILD_CHUNK tasks, plus a MERGE_PROFILES cache
  collision. The outcome was never recorded in `.living/` — now captured.

**Consequences**: The living layer now reflects the 07-22→08-01 window. The `INDEX.md`
regeneration still fails under system `python3` (3.9) — see `learnings.md` "mycelium
generate_index.py needs Python 3.10+"; use `python3.12` (present in the pixi env and
`/usr/bin`) until the hook is fixed upstream.

**Tags**: mycelium, reconciliation, decisions, novelty-discovery, singleton, context-search,
roles, ev-score, docker, pfam, perf, findings, in-progress

---

## Report presentation: one skin registry, one linkout builder, one landing-page design (2026-09-05)

**Context**: A code audit of the report/`view/` layer, driven by a request for a
cyberpunk/Neuromancer look that could be *swapped* rather than hardcoded, plus better
external-database linkouts for a broad audience reading these off GitHub Pages.

**What we found**: the presentation layer had four independent copies of the same colour
tokens (inline in `report_template.py`, `report_common.THEME_VARS_CSS`, and twice as raw
hex in `view/generate_index.py`), plus a hardcoded light/dark pair for `losses.html`'s
`.badge.warn`. Adding a fifth palette meant editing four files and remembering the badge.
Two *different* generators also produced `view/<project>/report.html`
(`bin/make_index_report.py` from the config CSV; `view/generate_index.py` from the report
payloads) with different markup and CSS, so whichever ran last decided what a shared
folder looked like — and the top-level gallery matched neither.

**Decision**:

1. `lib/skins.py` is the sole source of colour. Selection is three-state — no `data-skin`
   attribute means follow the OS; an explicit value overrides it; `@media print` always
   pins Paper. Persisted in `localStorage` under one key so a choice carries across a
   project's three reports, and applied from `<head>` before first paint.
2. Skins own type and effects too (`--font-ui`/`--font-mono`/`--glow`/`--overlay`), so a
   neon look needs no `[data-skin=...]` selectors scattered through the templates.
3. `lib/report_common.py`'s `externalLinksNode()` is the sole external-link builder for
   all three reports.
4. `lib/index_page.py` renders both landing pages; the two generators keep their different
   *inputs* but no longer differ in what they draw.

**Alternatives rejected**:

- *Cyan/magenta for the neon skin* — the obvious cyberpunk pair, but magenta desaturates
  toward blue-grey under deuteranopia and converges with cyan, breaking the report's
  colour-role contract (series-1 = search presence, series-2 = TBLASTN hit). Cyan/amber is
  separated on the blue-yellow axis, which every common form of red-green CVD preserves.
- *A full `report_template.py` refactor onto `report_common.py`* — still rejected for the
  heatmap-specific chrome (large, well covered, no counterpart in the single-table pages).
  Only the two fragments that had actually **drifted** were shared: the skin CSS and the
  linkout helpers.
- *Recolouring the shipped `paper`/`dark` palettes* to widen the series-pair luminance gap.
  CLAUDE.md requires re-validating CVD separation first, and this change was meant to be
  pixel-neutral for the existing skins; the greyscale guard is calibrated to "no worse than
  what ships" (1.10) instead.
- *An `amber` monochrome VT220 skin* — cannot carry two evidence hues without a
  fill-pattern fallback. Dropped in favour of `contrast`.

**Consequences**: `tests/test_skins.py` makes the contrast floors a build failure rather
than an aspiration — it immediately caught white-on-neon presence-chip labels, which is
why `--on-series` exists. `tests/test_report_templates.py` runs the `node --check` pass
CLAUDE.md documented but nothing enforced; that matters more now that one bad shared
fragment breaks all four pages at once.

**Tags**: report, skins, accessibility, wcag, colour-blind, linkouts, github-pages,
refactor, duplication, config-csv
