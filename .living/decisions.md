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

9. **Fixed circular per-family presence threshold in `novelty_discovery` (2026-09-03)** —
   `bin/calibrate_family_hmms.py` set each family's presence E-value threshold to the best
   (lowest) E-value that family's HMM scored against DISCOVERY_OUT, then
   `bin/novelty_presence_matrix.py`/`bin/novelty_screen.py` required `evalue < threshold` to
   call presence. Since the threshold *is* the minimum observed DISCOVERY_OUT E-value, no
   DISCOVERY_OUT proteome could ever satisfy `E < min(E)` against itself — the
   outgroup-absence filter was a structural no-op. Confirmed against a real run
   (`results/sordario/`): all matrix rows had zero DISCOVERY_OUT presence, 83.8% of 14,779
   "candidates" had a contradicting TBLASTN hit. Flagged by a collaborator manually
   re-BLASTing three "novel" N. crassa genes (NCU00765, NCU00411, NCU01935) against S. pombe
   and finding highly significant (down to 3e-160) annotated hits. Independently verified by
   an agent code review before the fix landed.
   Fix: presence is now gated purely by the flat `--default-family-evalue` /
   `--hmm_presence_evalue` for every proteome (both discovery and screen phases); per-family
   calibration (`--family-thresholds`) is still computed and piped through for Nextflow
   process interface stability but is no longer consumed for gating. Also fixed a related bug
   in `lib/family_presence.py`'s `parse_domtblout()`: HMM-coverage was pooled across *all*
   target sequences in a proteome instead of per-target, letting several unrelated weak
   partial hits fake a passing coverage score — coverage is now computed per (query, target)
   and paired with whichever target gave the best E-value.
   **All prior `--cluster_tool novelty_discovery` run results (including `results/sordario/`)
   should be considered invalid for outgroup-absence and re-run.** The known-but-separate
   over-strict singleton paralog-cutoff filter (`lib/singleton_presence.py`, filter 1 — cutoff
   derived from a within-genome self-vs-self search, can become unreachable when a close
   in-paralog exists) was reviewed but left unfixed — it did not cause the three reported
   genes (they are family members, not singletons) and is out of scope for this fix.

**Tags**: novelty-discovery, bug-fix, correctness, family-hmm, calibration, presence-matrix

10. **Fixed H2: paralog-cutoff filter (filter 1) replaced with a flat significance floor,
    project-wide (2026-09-03)** — `lib/singleton_presence.py`, `bin/build_presence_matrix.py`,
    and `bin/context_presence.py` each independently implemented the same "paralog-cutoff"
    filter 1: hit e-value must beat the query's own within-genome paralog e-value (from a
    self-vs-self search reported down to E=100), falling back to a flat default otherwise.
    Measured on real N. crassa self-search data (`results/pezizo4/self_hits/Ncra.paralog_cutoffs.tsv`,
    8,293 proteins): 28.2% of cutoffs were looser than the intended default (1e-5) —
    self-search noise standing in for "no real paralog" — and 40.0% were tighter than 1e-50,
    an unreachable bar for any real cross-species ortholog. Confirmed concretely against real
    Ncra-vs-Spom phmmer AND diamond hits: **actin (NCU04173/Act1) and a P-type ATPase
    (NCU07966)** — both near-universally conserved — would have been wrongly called ABSENT
    from S. pombe under the old filter (their own near-identical in-genome paralog e-value
    was tighter than the real ortholog hit), in both tool outputs. All three genes from the
    original H1 report (NCU00765/00411/01935) were unaffected by this specific filter in
    this data (their own paralog cutoffs happened to be loose enough) — consistent with H1
    (family-HMM calibration circularity, decision #9) being the sole explanation for those.
    Fix: filter 1 is now a flat `--default-evalue`/`--singleton-evalue` applied to every hit
    (no per-query paralog derivation); filter 2 (paralog-competition, a direct head-to-head
    comparison) is unchanged and remains the real paralogy test. `load_paralog_info()` now
    returns only `paralog_of` (the e-value column is unused) in all three files.
    `bin/calibrate_family_hmms.py`-style self-referential calibration does not appear reused
    elsewhere (the general `--cluster_tool mmseqs` pathway already uses a flat
    `--hmm_presence_evalue`, not per-family calibration, per CLAUDE.md) -- a whole-project
    Fable consistency review was dispatched after this fix landed to verify that
    independently; see the follow-up entry once it returns.

**Tags**: novelty-discovery, bug-fix, correctness, paralog-cutoff, singleton, context-search,
build-presence-matrix

11. **Follow-up: Fable consistency review found and fixed a real gap in fix #9 (coverage
    gated on wrong target) (2026-09-03)** — a dispatched Fable review verified fixes #9/#10
    and audited the whole project for other pockets of either bug class. Findings:
    - **Confirmed correct and complete:** H2's flat-significance-filter fix (decision #10) in
      all three files; the general `--cluster_tool mmseqs` pathway (`bin/profile_to_matrix.py`)
      already used a flat threshold with no calibration step; TBLASTN thresholding is flat;
      the loss-search direction reuses the same (fixed) `build_presence_matrix.py`, no separate
      logic; the paralog-competition filter (filter 2) itself is sound and unaffected by either
      fix, in all four implementations.
    - **Real gap found in fix #9:** `lib/family_presence.py`'s `parse_domtblout()` (from
      decision #9) fixed the coverage-pooled-across-targets bug but introduced a narrower one --
      it evaluated the coverage gate only on whichever single target gave the BEST E-value,
      so a family whose best-E hit happened to be a coverage fragment (but a different,
      slightly-worse-E target had full coverage) was still wrongly called absent. Confirmed
      against real `results/sordario/` data: NCU00411 (ATG43) and part of NCU01935's
      DISCOVERY_OUT presence were still miscalled post-fix #9, not due to circularity anymore
      but due to this best-target-only coverage gate. Also found: `bin/profile_to_matrix.py`
      (the general `mmseqs` pathway) had its own, INDEPENDENT coverage bug -- summed raw
      domain spans without merging overlaps, so overlapping domains on the same target could
      double-count coverage past 1.0.
    - **Fix:** `parse_domtblout(path, default_evalue, min_coverage)` now returns, per query,
      the best E-value among targets where THAT SAME target independently clears both the
      E-value and coverage gates (present if ANY target qualifies -- matching
      `bin/profile_to_matrix.py`'s original "any target" semantics, just with correct merged-
      interval coverage). `bin/profile_to_matrix.py`'s own duplicate `parse_domtblout` was
      deleted entirely and replaced with an import of `lib/family_presence.py`'s version --
      the two family-presence pathways now share one implementation instead of two
      independently-maintained (and previously divergent) copies. Re-verified against real
      sordario data: NCU00765/00411/01935 are now all correctly present in at least one
      DISCOVERY_OUT proteome (so correctly excluded as novelty candidates); NCU00411 remains
      genuinely absent from S. pombe specifically (its only Spom hit has intrinsically low
      HMM coverage, ~26% -- a real biological/parameter-tuning question about a short protein,
      not a code bug) but is now correctly caught via its Ylip hit.
    - Also fixed for consistency: `--paralog-competition-scope`'s Python argparse default was
      `'proteome'` in all four scripts while `nextflow.config`'s pipeline default is `'target'`
      -- silently divergent for any standalone/manual script invocation (pipeline runs were
      unaffected, since Nextflow always passes the flag explicitly). Argparse defaults changed
      to `'target'` to match.
    - Updated stale documentation describing the removed paralog-cutoff/calibration logic as
      current: CLAUDE.md, README.md, METHOD_DESCRIPTION.md, nextflow.config comments,
      `bin/parse_self_hits.py`, `workflows/novelty_discovery.nf`, `workflows/search.nf`.
    202/202 tests pass (added regression tests for the any-target coverage fix and the
    scope-default change), ruff lint clean.

**Tags**: novelty-discovery, bug-fix, correctness, family-hmm, coverage, profile-to-matrix,
paralog-competition-scope, documentation, fable-review

12. **Added an absolute-residue alternative to the family-HMM coverage floor, and sweep
    infrastructure to actually validate it (2026-09-03)** — following decision #11's
    finding that `--hmm_presence_cov 0.5` (a flat coverage *fraction*) rejects 83% of
    otherwise-significant hits on HMMs >600aa (vs 7% for HMMs <150aa) as a structural
    side-effect of requiring coverage over the *whole* profile length, regardless of how
    much of a multi-domain protein's ortholog could realistically align:
    - Added `--hmm_presence_min_residues` (default 100, data-informed from a real-run
      sensitivity table but not yet swept/validated): a target now qualifies if EITHER
      the coverage fraction clears `--hmm_presence_cov` OR its merged aligned span is at
      least this many absolute residues — targeting the floor's actual stated purpose
      (reject a hit explained by one small, promiscuous shared domain) without penalizing
      long proteins purely for being long. Wired through `lib/family_presence.py`,
      `bin/novelty_presence_matrix.py`, `bin/novelty_screen.py`, `bin/profile_to_matrix.py`
      (whose own independent, buggy coverage-double-counting implementation was deleted
      and replaced with `lib/family_presence.py`'s shared, tested one), the corresponding
      Nextflow process/workflow wiring, and `nextflow.config`. A real bug was caught before
      landing: the naive `covered >= min_residues` with default `min_residues=0` is
      trivially always true, silently disabling the coverage floor entirely — guarded
      explicitly (`min_residues > 0 and ...`) and now covered by a dedicated regression
      test in every consuming file.
    - `--hmm_presence_cov` was previously never swept at all (unlike `--hmm_presence_evalue`)
      — its 0.5 default was a design-doc guess (ADR-0002 Q5), not empirically validated the
      way `--family_min_seq_id`/`--family_cov` were. `bin/run_param_sweep.sh` now includes
      `--hmm_presence_cov`/`--hmm_presence_min_residues` as swept dimensions, but the
      existing sweep quality metrics couldn't actually test them: `busco_recovery`
      (`bin/busco_family_recovery.py`) measures mmseqs *clustering* quality and is
      insensitive to presence-calling parameters entirely (they act downstream, after
      clustering); `score_controls.py`'s recall/fp_rate needs curated per-clade controls
      that may not exist. Added `bin/busco_presence_recovery.py`: a curation-free metric
      that checks, for each BUSCO ortholog recovered as one clean family, whether OUTGROUP
      (non-cluster-member) species BUSCO independently confirms carry a Complete copy are
      actually marked present in `presence_matrix.tsv` — the thing `hmm_presence_cov`/
      `hmm_presence_min_residues` directly control. A real implementation bug was caught
      before landing: the first draft required zero "unmapped" species among a BUSCO's
      copies to count it as scorable, but an unmapped copy (never fed to clustering) is
      *exactly* the outgroup case this script exists to test — smoke-testing against real
      sordario data (0 scorable pairs) surfaced this immediately; fixed by only requiring
      the *mapped* (cluster-member) copies to agree on one family, which unmapped/outgroup
      copies don't disqualify.
    - Addressed the user's methodological concern (are BUSCO genes systematically longer
      than the whole proteome, biasing any BUSCO-based tuning toward long genes only?)
      empirically with real *N. crassa* data (`busco_pezizo5/Ncra.busco/` vs its
      `family_hmmsearch` domtblout-derived protein lengths): BUSCO Complete gene lengths
      (median 460 aa, IQR 304–702, 32.9% >600aa) closely track the whole proteome's
      (median 462 aa, IQR 302–686, 32.7% >600aa) — BUSCO modestly *underrepresents* the
      shortest bucket (<150aa: 3.8% vs 6.7%), the one bucket that was already unaffected by
      the coverage floor, so low-risk for this specific use. Built this check into the tool
      itself (`--reference-domtblout`) rather than treating it as a one-off, so every future
      sweep run reports its own representativeness rather than assuming BUSCO generalizes.
    - `bin/collate_sweep.py`'s scoring extended with a `presence_recovery` admissibility
      gate + composite term, but only when the metric was actually measured (an
      all-missing column, e.g. no `BUSCO_OUTGROUP_TABLES` configured, is distinguished from
      a genuinely-measured-and-worst column and its gate is skipped rather than silently
      failing every grid point).
    - Extracted `bin/busco_family_recovery.py`'s BUSCO-loading helpers into `lib/busco.py`
      (shared with the new script) per this repo's shared-logic-goes-in-lib convention.
    218/218 tests pass, ruff lint clean. **No sweep has actually been executed** (needs
    real cluster time, BUSCO runs for outgroup species, and possibly curated controls) —
    this is infrastructure only; running it and wiring the recommended defaults back into
    `nextflow.config` is separate future work requiring the user's go-ahead to launch.

**Tags**: novelty-discovery, coverage-floor, hmm-presence, sweep, busco, adr-0002,
parameter-tuning

13. **Ran the pezizo5 coverage-floor sweep (scoped, 4 points); found and fixed two
    real bugs surfaced by the real run (2026-09-03)** — launched
    `run_sweep_pezizo5_coverage.sh` on the UCR HPCC `batch` partition (job 28051317),
    sweeping `--hmm_presence_cov` {0.5, 0.3} × `--hmm_presence_min_residues` {0, 100}
    against `configs/pezizo5.csv`, clustering params held at shipped defaults.
    - **Bug found and fixed before trusting results:** `bin/collate_sweep.py` treated an
      entirely-unmeasured optional metric (`fp_rate`/`recall`, no `configs/controls/
      pezizo5.controls.csv` exists) the same as a genuinely-measured worst value
      (`fp_rate=1.0`), which failed every grid point's admissibility gate and printed a
      misleading "no grid point met the gates" warning implying a 100% false-novelty
      rate rather than "not measured." Generalized the same guard already built for
      `presence_recovery` (decision #12) to `recall`/`fp_rate` too, added a `run_ok`
      gate/column so a failed pipeline run's partial metrics can never be admissible,
      and added regression tests for both — one of them (`test_recall_and_fp_rate_gates
      _skipped_when_never_measured`) directly encodes the real numbers from this sweep.
    - **Second bug found and fixed, also before trusting results:** the metrics TSV
      itself was silently corrupted — Python's `csv` module always writes `\r\n` line
      endings unless `lineterminator='\n'` is set explicitly (`newline=''` on `open()`
      does NOT prevent this, a common gotcha), so `bin/busco_family_recovery.py`'s
      `*.summary.tsv` had `\r\n` endings; `awk`'s `$2` extraction in
      `bin/run_param_sweep.sh` left a trailing `\r` baked into the shell variable; and
      Python's `csv.DictReader`, run on the resulting metrics TSV, treats a bare `\r`
      as its own row boundary (regardless of the file's declared dialect) — silently
      splitting every 1 real row into 2 corrupt ones with fields shifted into the wrong
      columns. Fixed at both ends: `lineterminator='\n'` added to every `csv.writer`/
      `csv.DictWriter` call in `busco_family_recovery.py`, `busco_presence_recovery.py`,
      `score_controls.py`, and `collate_sweep.py` itself; `tr -d '\r'` added defensively
      to every `awk`-based field extraction in `run_param_sweep.sh` as a second line of
      defense. The already-collected `sweep_metrics.tsv` was cleaned in place (`tr -d
      '\r'`) rather than re-running the (now-fixed-at-the-source) pipeline jobs.
    - **A separate, unrelated infrastructure failure** (not a bug in this session's code)
      also hit the first launch attempt: the Pfam annotation step's `mpirun -np 32
      hmmsearch --mpi` requested more MPI slots than the SLURM allocation had. Fixed
      pragmatically by adding an `ANNOTATE_SWEEP` toggle (default off) to skip Pfam/
      SwissProt annotation for the sweep entirely — none of the sweep's metrics need it,
      it is the single most expensive step, and it was what broke.
    - **Real (if partial) result, from `tblastn_removed`** — candidates with a
      *contradicting* outgroup TBLASTN genomic hit despite being called protein-level
      absent, the same red flag the original H1/H2 bug reports were built from:
      current shipped default (`hmm_cov=0.5, hmm_residues=0`): 1961 novelty candidates,
      **554 contradicted (28.2%)**. Loosest tested setting (`hmm_cov=0.3,
      hmm_residues=100`): 1107 candidates, **21 contradicted (1.9%)** — a 15-fold drop
      in the contradiction rate for a 44% drop in raw candidate count, i.e. the
      candidates being removed are disproportionately the low-quality (contradicted)
      ones, not a uniform thinning. `busco_recovery` (0.753) was identical across all
      4 points as expected (it measures clustering, not presence-calling, so it can't
      discriminate between these settings — confirming decision #12's design rationale
      for needing `presence_recovery` instead).
    - **What this run does NOT establish**: no outgroup BUSCO tables exist for pezizo5
      (only the 5 IN-group species have been BUSCO-run, under `busco_pezizo5/`) and no
      controls CSV exists either, so `presence_recovery`/`recall`/`fp_rate` are all
      unmeasured this round — ranking rests on `busco_recovery` (insensitive to the
      parameters under test) plus `tblastn_removed` read by eye, not a rigorous
      composite score. Running BUSCO against an outgroup proteome (e.g. Spom, Scer) and
      re-running with `BUSCO_OUTGROUP_TABLES` set is the natural next step before
      actually changing `nextflow.config`'s shipped defaults.
    221/221 tests pass, ruff lint clean. `nextflow.config`'s defaults were NOT changed
    by this decision — `hmm_presence_cov=0.5`/`hmm_presence_min_residues=100` remain as
    set in decision #12; this entry documents the sweep run and what it found, not a
    new shipped default.

**Tags**: novelty-discovery, coverage-floor, hmm-presence, sweep, busco, slurm, mpi,
csv-crlf, bug-fix

14. **Built the broader-grid validation infrastructure for the coverage-floor sweep
    (2026-09-03)** — following decision #12/#13's finding that the pezizo5 sweep was
    only 4 points on one clade, added the machinery for real divergence-depth diversity
    per `todo/validate-hmm-presence-coverage-broader-sweep.md`:
    - `bin/convert_1kfg_samples.py` (new, parallels `bin/convert_bfd_samples.py`): joins
      1KFG's (1000 Fungal Genomes, JGI MycoCosm) 345 flat-directory proteomes under
      `/bigdata/stajichlab/shared/projects/1KFG/genomes/final_combine/` against
      Fungi_BFD's `samples.csv` (23,684 assemblies, ~12,401 annotated) by genus+species
      epithet to recover taxonomy 1KFG doesn't carry on its own (~70% match rate, the
      rest skipped and reported, not guessed) — `config_support/1KFG_samples.csv`, 241
      species with full semicolon-delimited Lineage strings, ready for
      `bin/make_config.py`'s existing `--max-per-outgroup-taxon`/`--random`/`--seed`
      stratified sampling.
    - Two new configs drafted from that pool: `configs/sordariales_shallow_1kfg.csv`
      (shallow divergence — ingroup Sordariales: Chaetomium/Neurospora/Podospora/
      Sordaria, 4 species one order; outgroup sibling Sordariomycetes orders, capped
      3/order, seed 42 — 14 proteomes) and `configs/deep_broad_1kfg.csv` (deep + large —
      ingroup 105 species across 7 filamentous-Ascomycota classes, since BFD doesn't
      populate SUBPHYLUM for most rows so "Pezizomycotina" isn't a matchable token and
      the ingroup was built via explicit `--ingroup-short` pins instead; outgroup a
      stratified random sample, 5/phylum seed 42, across 11 other fungal phyla,
      excluding Microsporidia as a reduced-genome outlier — 130 proteomes total, real
      ADR-0002/Chaetothyriales scale, requires `--cluster_tool mmseqs`).
    - `run_busco_1kfg_array.sbatch`: a SLURM array job (130 tasks, `%25` throttle) running
      BUSCO fungi_odb12 once for the full union of both new configs' species (config 1's
      14 species turned out to be entirely a subset of config 3's 105 ingroup, so one
      130-species BUSCO batch covers both configs' `BUSCO_TABLES`/`BUSCO_OUTGROUP_TABLES`
      needs) — output under `busco_1kfg/<Short>.busco/`.
    - `run_sweep_sordariales_shallow.sh` / `run_sweep_deep_broad_1kfg.sh`: new sweep
      launchers mirroring `run_sweep_pezizo5_coverage.sh` (same scoped 2×2
      `hmm_presence_cov`×`hmm_presence_min_residues` grid, clustering params held at
      shipped defaults, annotation skipped) but building `BUSCO_TABLES`/
      `BUSCO_OUTGROUP_TABLES` dynamically from each config CSV via `awk` rather than
      hand-transcribing 105+25 `SHORT=path` entries. `deep_broad_1kfg`'s launcher sets an
      explicit `--time=14-00:00:00` (`batch` partition default is 7 days, max 30) since
      it is an order of magnitude larger than pezizo5 (130 vs 11 proteomes) and the
      driver's own SLURM job must stay alive for the entire nested `nextflow -profile
      slurm` orchestration, not just its own work.
    All three sweeps (pezizo5/medium, sordariales_shallow/shallow, deep_broad_1kfg/deep+
    large) launched together to compare how the recommended `hmm_cov`/`hmm_residues`
    setting shifts (or doesn't) across divergence depth and panel size — results to be
    logged in a follow-up decision once they complete.

**Tags**: novelty-discovery, coverage-floor, sweep, busco, 1kfg, bfd, taxonomy,
broader-grid, slurm-array

15. **Excluded Penicillium chrysogenum from configs/deep_broad_1kfg.csv: a genuine
    duplicate-ID ambiguity in its own 1KFG proteome FASTA (2026-09-03)** — discovered when
    it broke both its own BUSCO run (`Duplicate of sequence >Pchr|PCH_Pc18g05710.1`) and
    then the deep_broad_1kfg sweep's `EXTRACT_FAMILY_SEQS` step
    (`bin/extract_family_seqs.py` -> `lib/fasta.py`'s `read_fasta()` -> `Bio.SeqIO.to_dict()`,
    which raises `ValueError: Duplicate key` on any repeated FASTA header). Checked before
    "fixing" it by deduping: `Pchr|PCH_Pc18g05710.1` has 2 exact-duplicate records (would
    have been safe to collapse), but `Pchr|PCH_Pc17g01280` has 4 records under the same ID
    that are NOT sequence-identical (same length, different residues) -- a genuine upstream
    annotation-export ambiguity, not a harmless duplicate. Rather than guess which variant
    is "correct" or edit the shared `/bigdata/stajichlab/shared/projects/1KFG/` resource
    (used by other projects), removed the one species (105 -> 104 ingroup species)
    from `configs/deep_broad_1kfg.csv` and let the already-running sweep self-heal via
    `-resume` (Nextflow re-reads the CSV fresh each grid-point invocation within the bash
    loop, so the content change is picked up without restarting the job).

**Tags**: novelty-discovery, data-quality, 1kfg, duplicate-id, bug-fix

14. **Two real infrastructure bugs found and fixed via the broader-grid sweep runs
    (2026-09-04)** — running `sordariales_shallow_1kfg`/`deep_broad_1kfg` at real scale
    (14 and 130 proteomes respectively, from 1KFG genomes) surfaced two genuine,
    previously-latent bugs, distinct from the presence-calling logic bugs (#9-#13):
    - **`lib/fasta.py`'s `read_fasta()` crashed on duplicate sequence IDs** (used
      `Bio.SeqIO.to_dict()`, which requires global uniqueness) when concatenating
      proteomes from independently-sourced genome collections — confirmed real:
      `deep_broad_1kfg`'s loss-direction family building hit `ValueError: Duplicate key
      'Pchr|PCH_Pc18g05710.1'` from two different genomes reusing the same short internal
      locus tag as their protein ID, killing the whole `PROFILE_LOSS_SEARCH` run. Fixed to
      keep the first occurrence and warn to stderr instead of crashing — a real
      possibility whenever proteomes come from heterogeneous sources (JGI/1KFG genomes in
      particular), not just a one-off oddity of this test set.
    - **`conf/ucr_hpcc_slurm.config`'s `hmmsearch` label had a flat `memory = '8.GB'`**,
      not scaled by `task.attempt` like `queue`/`time` already are — so `errorStrategy`'s
      retry (up to `maxRetries: 2`) was a structural no-op for a genuinely memory-hungry
      chunk: all 3 attempts hit the identical 8GB OOM wall (exit 137), confirmed on
      `sordariales_shallow_1kfg`'s `HMMSEARCH_CHUNK`. Fixed to `{ 8.GB * task.attempt }`
      (8/16/24GB across attempts), the standard Nextflow retry-scaling pattern. The
      `high_cpu` label (`memory = '32.GB'` flat) has the same latent risk and wasn't
      touched — no confirmed failure there yet, flagging for awareness rather than
      preemptively changing an untested path.
    Both fixes are infrastructure-only; no test coverage exists (or is practical) for the
    SLURM config scaling behavior itself, but `lib/fasta.py` has new unit tests
    (`tests/test_fasta.py`) for the duplicate-ID tolerance. 223/223 tests pass, lint clean.
    Both already-failed grid points (deep_broad_1kfg gp1, sordariales_shallow_1kfg gp2)
    will be manually resumed with `-resume` once their sweep driver's current pass
    finishes and its launch directory is free (avoids a concurrent Nextflow session-lock
    collision within the same sweep, mirroring the earlier per-CLADE `LAUNCH_DIR`
    isolation fix's rationale).

**Tags**: bug-fix, infrastructure, fasta, slurm, oom, retry, hmmsearch, novelty-discovery,
sweep

16. **Sweep-launcher isolation fixes: per-CLADE Nextflow launch dirs and actually loading
    the site SLURM config (2026-09-04)** — resuming the broader-grid sweep work after a
    session crash/restart surfaced a cluster of related launcher-level bugs, distinct from
    both the presence-calling logic bugs (#9-#13) and the infrastructure bugs already
    logged in #14:
    - **`conf/ucr_hpcc_slurm.config` was never actually loaded by any sweep launcher.**
      `bin/run_param_sweep.sh` invoked `nextflow run main.nf -profile slurm ...` but never
      passed `-c conf/ucr_hpcc_slurm.config`, so none of the site-specific tuning in that
      file (queue routing to `short`/`epyc`, the `hmmsearch` memory-scaling fix from #14,
      AVX2 node constraints for mmseqs/famsa, `BUILD_CHUNK`'s `preempt`-queue placement)
      was ever active during any of the sweeps run this session — `-profile slurm` alone
      only sets `process.executor = 'slurm'` from `nextflow.config` itself. Fixed by adding
      an `NXF_SITE_CONFIG` env var (default `conf/ucr_hpcc_slurm.config`) and passing
      `-c "$NXF_SITE_CONFIG"` whenever it's non-empty. This was the single biggest
      root-cause finding of the session: several already-diagnosed-and-"fixed" problems
      (OOM retries never actually recovering, mmseqs SIGILLs on incompatible nodes) had
      their real fix sitting unused in a file nobody was loading.
    - **Nextflow session-lock collisions between concurrently-running sweeps.** `-resume`
      with no explicit run name resumes "most recent run" keyed to the launch (current
      working) directory, not `--project`/`--outdir` — two sweeps launched from the same
      repo root raced for the same `.nextflow/history` entry, and one's `-resume` tried to
      attach to the other's actively-running session. Fixed by isolating each CLADE's
      launch into its own `.nf_launch/${CLADE}/` subdirectory (all path-bearing variables
      resolved to absolute paths first via new `abspath`/`abspath_tables` bash helpers,
      then `cd` into the isolated dir before invoking `nextflow run "$REPO_ROOT/main.nf"
      ...`). `.nf_launch/` added to `.gitignore`.
    - **Stale/parallel session found mid-recovery.** While resubmitting sweeps after the
      launcher fixes above, discovered evidence of a second, uncoordinated Claude Code
      session having independently touched some of the same files (a divergent
      `run_sweep_deep_broad_1kfg.sh`, a `conf/ucr_hpcc_slurm.config` dated ahead of this
      session's own edits) and a duplicate SLURM submission 9 seconds after this session's
      own resubmission. Confirmed via `ps aux` (3 `claude` processes); user identified and
      stopped the dormant/orphaned one. Take-away for future sessions on this project:
      check `ps aux | grep claude` and `squeue` for pre-existing jobs/processes before
      relaunching a sweep after any session restart, since `.nf_launch/` isolation
      prevents Nextflow-level lock collisions but not duplicate `sbatch` submissions from
      independent sessions.
    Both `sordariales_shallow_1kfg` (job 28100809) and `deep_broad_1kfg` (job 28096310)
    are now running cleanly with the fixed launcher, confirmed via `.nextflow.log` to be
    loading `-c conf/ucr_hpcc_slurm.config`.

**Tags**: novelty-discovery, infrastructure, nextflow, slurm, session-lock, sweep,
site-config, bug-fix

---

### [2026-09-05] Report presentation: one skin registry, one linkout builder, one landing-page design

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
