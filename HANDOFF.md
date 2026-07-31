# Handoff — 2026-07-26

## Right now: a live SLURM run is in flight

`run_Chaeto_subset.sh` (`--config configs/Chaetothyriales_subset.csv --cluster_tool mmseqs
--project Chaeto_subset`) is actively executing on the cluster. As of 07:38 this session,
`PROFILE_SEARCH:BUILD_FAMILY_PROFILES:BUILD_CHUNK` and
`PROFILE_LOSS_SEARCH:BUILD_FAMILY_PROFILES:BUILD_CHUNK` tasks (famsa+hmmbuild per family
chunk) are running/queuing on the `preempt` partition — ~100 tasks in flight, most `PD`
on `AssocGrpCpuLimit` (queue is throttled, not stuck). **Do not resubmit** this job; check
progress with `squeue -u $USER | grep PROFI` and `tail -f logs/chaeto_subset.log` /
`.nextflow.log`, or just resume watching — nextflow's own `-resume` will pick up any
already-completed chunks if it needs restarting.

When it finishes, outputs land in `results/Chaeto_subset/` and get collated to
`view/Chaeto_subset/`. A prior run already produced `view/Chaeto_subset/{report,core,losses,
novelties}.html` (committed in `cdc082a`) — this run is likely a rerun after the family-size
fix in `773ecf6` ("build out family profiles but deal with too large families") and the
`extract_family_seqs.py` fix in `f6f5476`, so expect the outputs to change.

## Uncommitted work — needs a decision before it's lost

`git status` shows modified-but-uncommitted:

- `configs/Chaetothyriales.csv` and `configs/Chaetothyriales_subset.csv` — removed 3 IN
  genomes with no matching data (`Exophiala spinifera JCM 15939`,
  `Herpotrichiellaceae sp. CTeuk-1624`, and one more) and **reclassified 3 species from
  IN → OUT** (`Aspergillus fumigatus Af293`, `Zymoseptoria tritici`,
  `Coccidioides immitis WA_211` — these are Pezizomycotina outgroup markers that were
  miscategorized as Chaetothyriales ingroup). This is exactly the kind of edit that should
  be committed once you're satisfied it's correct — check `git diff configs/Chaetothyriales*.csv`
  before committing, since the live run above was launched against whichever version of
  these files was on disk when nextflow read them (confirm which one that was before
  trusting the run's output).

There are also many **untracked new files** from recent work that have never been committed:
`bin/convert_bfd_samples.py`, `config_support/`, `configs/pezizo5.csv`, `db/modelorgs/*`,
`db/pfam/`, `db/uniprot/`, `docs/make_samples_intermediate.md`, `run_Chaeto.sh`,
`run_Chaeto_subset.sh`, `run_pezizo4.sh`, `steps_Chateo.md`, `tests/test_extract_family_seqs.py`,
plus `.living/` and `.claude/` (the mycelium scaffold + hooks — see CLAUDE.md, these are
meant to be checked in). `example/`, `stories/`, `skillpacks/`, `busco_pezizo5/`, `data/`,
and `logs/` are almost certainly working-tree/data directories, not things to commit — sanity
check each before staging anything broadly (don't `git add -A`).

## Where the ADR-0002 family-profile pathway stands

- **Phase 1 (mmseqs+HMM profile pathway, issue #3) is merged** (PR #4, `.living/decisions.md`
  2026-07-21). `todo/TODO_REGISTRY.md` still lists it `in-progress` — flip it to `complete`
  next time you touch that file.
- **AVX2 famsa crash — already fixed.** `.living/learnings.md` and `.living/decisions.md`
  record `famsa` SIGILLing on the cluster's Abu Dhabi (`c[01-30]`) nodes, deferred as a known
  issue on 2026-07-21. It's since been fixed: `nextflow.config`'s `slurm` profile now has a
  `withName: '.*BUILD_CHUNK'` block pinning to `-C ryzen|broadwell|cascade` on the `preempt`
  queue. No action needed here — but the `.living/` entries are now stale and should get a
  short "resolved" follow-up note pointing at the `nextflow.config` fix, so a future session
  doesn't re-discover the same bug.
- **"Too large families" handling** was being worked on as of the last two commits
  (`773ecf6`, `f6f5476`, both 2026-07-25 20:39) — `modules/build_family_profiles.nf` gained
  oversized-family skip logic, `bin/extract_family_seqs.py` gained matching logic to pull
  family members out of the mmseqs cluster. **This hasn't been validated against a completed
  run yet** — that's what the in-flight Chaeto_subset job above is for. Once it completes,
  check `results/Chaeto_subset/*/oversized_families.tsv` (per `family_max_members` config
  comment in `nextflow.config`) to see how many families got skipped and whether that's
  reasonable for this dataset.
- **Phase 2** (`todo/cross-method-support-column.md` — sweep + BUSCO/controls validation
  battery) is still `open`, not started. The full end-to-end pezizo5/Chaetothyriales parity
  check between `--cluster_tool pairwise` and `--cluster_tool mmseqs` (shape-parity,
  reports render, BUSCO/negative-control recall) called for in Phase 1's acceptance
  criteria has still never been run — this Chaeto_subset job is the closest thing to that
  so far, but it's the mmseqs path only, not a pairwise-vs-mmseqs comparison.

## Suggested next steps, in order

1. **Let the running job finish** (or check on it) before doing anything else that touches
   `results/Chaeto_subset/` or `work/`.
2. **Resolve and commit the `configs/Chaetothyriales*.csv` edits** — confirm the 3
   dropped-for-missing-data rows and 3 IN→OUT reclassifications are intentional and final,
   then commit with a clear message (e.g. "fix: drop genomes with no annotation data, move
   3 misclassified Pezizomycotina markers from ingroup to outgroup").
3. **Triage the untracked files** — commit the real deliverables (`bin/`, `docs/`,
   `run_*.sh`, `tests/`, `configs/pezizo5.csv`, `.living/`, `.claude/`), decide whether
   `db/modelorgs/*`, `db/pfam/`, `db/uniprot/` belong in git or should be `.gitignore`d
   (they look like large reference databases — check size before committing), and leave
   `data/`, `logs/`, `example/`, `stories/`, `skillpacks/`, `busco_pezizo5/` alone unless you
   know they're meant to be tracked.
4. **Update `.living/decisions.md`/`learnings.md`** with a short "AVX2 fix shipped" note and
   flip `todo/TODO_REGISTRY.md`'s Phase 1 row to `complete`, since both are currently stale
   relative to what's actually in `nextflow.config` and merged to `main`.
5. **Once Chaeto_subset is validated**, decide whether to scale the same `--cluster_tool
   mmseqs` run to the full `configs/Chaetothyriales.csv` (`run_Chaeto.sh` currently runs
   `--cluster_tool pairwise` — the default — not mmseqs; you'd need a variant or a flag
   change to actually exercise the profile pathway at full scale).
6. **Start Phase 2** (`todo/cross-method-support-column.md`) once Phase 1 is confirmed
   working end-to-end: the sweep grid, BUSCO/negative controls, and pairwise-vs-mmseqs
   cross-method concordance check are all still outstanding and are what will actually
   validate whether the mmseqs pathway's presence calls are trustworthy.

## Quick reference

- Live run: `squeue -u $USER | grep -i profi`, logs at `logs/chaeto_subset.log` and
  `.nextflow.log`.
- Config diff pending decision: `git diff configs/Chaetothyriales.csv
  configs/Chaetothyriales_subset.csv`.
- ADR: `docs/adr/0002-family-profile-search-pathway.md`.
- Todos: `todo/TODO_REGISTRY.md`, especially `profile-search-pathway.md` (Phase 1, needs
  status flip) and `cross-method-support-column.md` (Phase 2, next up).
- Session log for this session: `.living/log/2026-07-26-001-novinvenio.md`.
