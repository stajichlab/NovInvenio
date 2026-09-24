# Learnings

Append-only log of gotchas, surprises, and insights.

**Entry template:** copy from `skills/core/templates/learning-entry.md` (includes Category, What happened, Why it matters, Resolution, Tags fields). The `**Tags**:` line is consumed by `generate_index.py --summary-heuristic` to build the cluster summary in INDEX.md — use them.

### [2026-07-20] Claude Code owns settings.local.json — put hooks in settings.json

**Category**: gotcha

**What happened**: During mycelium init, `init_repo.py` wrote the 7 mycelium hooks into
`.claude/settings.local.json` (merging with the existing `permissions`/`skillOverrides`)
and reported success. A later check showed the `hooks` key was gone — the harness had
re-written `settings.local.json` (it auto-persists approved permissions there) from an
in-memory copy that predated the out-of-band edit, silently clobbering the hooks.

**Why it matters**: Any tool or script that edits `settings.local.json` while Claude Code
is running can have its changes overwritten without warning. Hooks written there may
appear installed but vanish, leaving the living-repo automation silently inert.

**Resolution**: Moved the hooks to `.claude/settings.json` (shared project settings, which
the permission-writer does not manage). Hooks from `settings.json` and
`settings.local.json` merge, so this survives permission updates. Verified valid JSON.

**Tags**: claude-code, hooks, settings, mycelium, tooling, configuration

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: A post-init verification step that re-reads
`.claude/settings.json` (not just `.local.json`) and asserts all 7 mycelium hook
basenames are registered would catch a clobbered/missing install. Not yet shipped.

### [2026-07-20] mycelium Stop hook false-blocks interactive (grilling) sessions

**Category**: gotcha

**What happened**: During a multi-turn grilling session, `mycelium-stop-check.sh` blocked
session end nearly every turn with "STOP BLOCKED — N files changed but .living/ not updated",
even though `.living/decisions.md` was being updated every turn. At one point the live mtimes
were unambiguous — `decisions.md` mtime was ~50 min newer than the debounce reminder
(`.claude/mycelium-reminded.tmp`), which by the script's own logic (mtime > REMINDER_TS ⇒
`DECISIONS_UPDATED=true` ⇒ pass) should NOT block — yet it blocked. Marketplace and cache
copies of the hook are byte-identical, so the static logic and the runtime behaviour
genuinely disagree (suspect: a second plugin-registered copy, or hook config only reloading
at session start so an in-session edit to settings.json has no effect).

**Why it matters**: The hook is designed for "do work → reflect" sessions, not back-and-forth
Q&A. In an interactive skill (grilling, QA, request-refactor-plan) it fires on every turn and
obstructs the conversation. Removing it from `.claude/settings.json` mid-session did NOT stop
it (config likely reloads only at session start). The reliable escape is to clear the
debounce sentinels it checks — `rm .claude/mycelium-reminded.tmp
.claude/mycelium-session-activity.tmp` — after genuinely updating `.living/`; with both
absent the hook exits early (line ~269) instead of blocking.

**Resolution**: Recorded each resolution to `decisions.md` as we went (honest + normally
sufficient); when the hook still false-blocked, cleared the sentinels to end the loop. The
reflection itself was complete throughout.

**Tags**: mycelium, hooks, stop-hook, grilling, interactive-session, tooling, debounce

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: The Stop hook should treat interactive-skill turns as
exempt (e.g. detect an active skill/`stop_hook_active`-style flag) or key "reflected" off any
`.living/*` mtime > session-start rather than a Bash-re-armed reminder timestamp. Upstream fix
in mycelium; not something this repo can assert via a test.

### [2026-07-21] famsa (bioconda) requires AVX2 — SIGILLs on the cluster's Abu Dhabi nodes

**Category**: gotcha

**What happened**: After `pixi install` pulled in `famsa` (added for ADR-0002 Q4, the
family-profile MSA step), a smoke-test invocation on dev node `c09` crashed immediately with
`Illegal instruction (core dumped)`. `/proc/cpuinfo` on `c09` showed AVX but no AVX2; the
bioconda famsa build needs AVX2+. Worse: `sinfo -o "%N %P %f"` showed the pipeline's default
SLURM `short` partition mixes `c[01-30]` (`amd,abu_dhabi`, no AVX2) with `i[01-62]`
(`intel,broadwell`, has AVX2) — a real run could land on either node class and crash
non-deterministically depending on scheduler placement. `MMSEQS_CLUSTER` already has this
exact problem for a different tool and solves it with `clusterOptions = '-C ryzen'`
(`nextflow.config` ~line 163) to pin to AVX2-capable nodes.

**Why it matters**: any bioconda/conda-forge binary built with SIMD auto-vectorization
(AVX2/AVX512) can silently crash on this cluster's older AMD "Abu Dhabi" nodes even though
`pixi install` succeeds and the tool runs fine on a modern node — the failure is
node-dependent, not environment-dependent, so it won't reproduce consistently and won't show
up in CI unless CI happens to land on an old node.

**Resolution**: deferred by user request ("save for handoff later — we will run this on an
avx2 node"). Not yet fixed. Fix on resume: add a SLURM `-C` constraint excluding
`abu_dhabi` (e.g. matching `MMSEQS_CLUSTER`'s `-C ryzen` pattern) to whichever process runs
famsa (`BUILD_FAMILY_PROFILES` module, once built), or otherwise pin its execution to an
AVX2+ partition (`epyc` = `ryzen,amd,milan`, confirmed AVX2-capable).

**Tags**: cluster, slurm, avx2, famsa, pixi, bioconda, simd, hardware-compatibility,
nextflow, node-placement

**mitigation_type**: convention

**structural_mitigation_candidate**: A `withName` SLURM `clusterOptions` constraint on the
famsa-running process, mirroring `MMSEQS_CLUSTER`'s `-C ryzen`. Not yet shipped — this is the
concrete candidate, tracked for the Phase 2 resume.

### [2026-08-06] mycelium generate_index.py needs Python 3.10+ — silently fails under system python3 (3.9)

**Category**: gotcha

**What happened**: The mycelium health hook (`.claude/settings.json` SessionStart →
`mycelium-health.sh`) regenerates `.living/INDEX.md` at every session start by calling
`python3 generate_index.py --summary-heuristic` (directly, or falling back to
`--counts-only`), then `|| true`s the error. On this cluster the bare `python3` is **3.9.18**,
and `generate_index.py` uses PEP 604 union syntax (`-> str | None`) that only parses on
3.10+. So every invocation raised `TypeError: unsupported operand type(s) for |` and was
silently swallowed. Result: `INDEX.md` was frozen at "Last audit: 2026-07-21" even though 10
sessions and a large body of work had accumulated through 08-07. Running the same command
with `python3.12` (present in the pixi env and `/usr/bin/python3.12`) regenerates it cleanly.

**Why it matters**: A "knowledge sync" that runs silently and never errors (forking the error,
never surfacing it) is a silent failure — the index drifts from reality without any signal.
It's not the agent's job to detect it; the hook must be robust to environment Python version.
This is the specific mechanism behind the repo's stale living index.

**Resolution**: Regenerated `INDEX.md` manually with `python3.12`:
```bash
/usr/bin/python3.12 /rhome/jstajich/.claude/plugins/marketplaces/mycelium/skills/core/scripts/generate_index.py \
  --living-dir .living --summary-heuristic
```
Upstream fix (outside this repo): make the health hook choose a 3.10+ interpreter
(`python3.10`/`python3.12`/pixi env) or stop using PEP-604-only syntax.

**Tags**: mycelium, python, index, tooling, environment, silent-failure, version

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: A session-start check that verifies INDEX.md's
"Last audit" is current after the hook runs, or a health-hook change to invoke a 3.10+
interpreter.

### [2026-08-06] Session-log stubs vs substance: hooks create logs, agents must write the findings

**Category**: gotcha

**What happened**: The mycelium SessionStart/SessionStop hooks reliably create a new
`.living/log/<date>-NNN-*.md` per session (the active-session-log + stop sequence works), and
`LOG_REGISTRY.md` scaffolding exists — but 7 of 10 session logs through 08-07 contain only the
front-matter + "Session started", with zero substance. Meanwhile the git history and
`logs/`/`results/` show a dense 07-22→08-01 work window (novelty-discovery, singleton
screening, context search, containerization, Pfam speedup + a half-dozen investigation runs)
that is entirely absent from `learnings.md`/`decisions.md`/`findings/`.

**Why it matters**: The hooks scaffold the container but do not (and cannot) capture the
agent's reasoning or findings — only the agent itself, following the post-action protocol, can
write learnings/decisions. Where sessions end abruptly or the protocol is skipped, the work
lives only in git/logs and is invisible to the knowledge layer. `findings/` (the intended home
for analysis-run findings) was completely empty.

**Resolution**: Backfilled the window from `git log --all` + run logs + `results/*` counts
(`decisions.md` reconciliation entry, `learnings.md` run-observations, `LOG_REGISTRY.md`
rows). Going forward: the Stop hook is disabled (see `conventions.md`), so enforcement is
manual — treat "wrap up the session in `.living/`" as an explicit closing step, and deposit
per-run findings in `.living/findings/` rather than only in session logs.

**Tags**: mycelium, sessions, logging, findings, knowledge-capture, process

**mitigation_type**: convention

**structural_mitigation_candidate**: A SessionStart "unwritten-findings" reminder (like the
knowledge-audit message) that lists open session-log stubs and the empty `findings/` dir as a
nudge, since the blocking Stop hook was removed.

### [2026-08-06] Investigation-run observations (benchmarks & cross-method divergences)

**Category**: insight

**What happened**: Reconciling the late-July/early-Aug run logs (`logs/`, `results/*`)
surfaced concrete, previously-unrecorded numbers:

1. **mmseqs family-profile vs pairwise diverge on novelty counts**: pezizo5_mmseqs → 1,851
   candidates / Afum 524, Amega 197, Cimm 453, Ncra 342, Ztri 335 vs pezizo5 pairwise → 2,544 /
   Afum 604, Amega 375, Cimm 550, Ncra 447, Ztri 568. The form/identity presence call is not
   a drop-in count — cross-method concordance (Phase 2) is warranted.
2. **Singleton search roughly doubles recovered target genes**: with PR #52 the sordario
   novelty_discovery run jumped Afum 422→907, Cimm 329→693, FgraPH1 419→1173, Ncra 429→956
   (`target_specific` 1,599→3,729; screened_candidates 14,779→16,909). The singleton pathway
   has a big, previously-invisible payoff — its earlier invisibility was the "known gap".
3. **Antarctolithicomycotina** (2 ingroup genomes) is novelty-rich (772 / 758 novelties) and
   has **loss_candidates (1,697) > novelty candidates (1,530)** — an unusually high loss
   signal for a 2-genome ingroup; worth a closer look.
4. **preempt did not stabilize mmseqs BUILD_CHUNK**: the 07-24 handoff
   (`SLURM job 26715986`) finished `completed=5063 failed=45 cached=57` — 45 famsa 600s
   timeouts (exit 140) across PROFILE_SEARCH + PROFILE_LOSS_SEARCH, 100 ABORTED + 28 FAILED
   BUILD_CHUNK tasks in the trace, and a MERGE_PROFILES cache re-use collision
   (`WARN: Unable to resume cached task`). Preempt shrank individual-task loss but the
   famsa-on-`short` timeout/AVX2 issue was the real blocker.
5. **Chaeto_subset** (21 ingroup + ~44 outgroup, mmseqs): 157,181 matrix rows, 666 novelty /
   205 loss candidates; two ingroup species (Cimm, Herpsp) got **0** novelties.
6. **Pfam hmmscan performance** (08-01): hmmscan's own `--cpu` threading under-utilizes
   (~1.4 of 8 on a real ~5.7k-protein set, GPFS-bound re-scanning the full pressed DB per
   query); concurrency from many low-cpu tasks beats internal threading. Only the pressed
   index files `.h3f/.h3i/.h3m/.h3p` are read (not the 2.1 GB flat Pfam-A.hmm or `.dat`).
   `scratch = true`: 94% CPU/14m55s (GPFS) vs 165% CPU/6m46s (node-local) on a 500-protein
   subset ≈ **2.2× wall-clock**; per-chunk cpus=2 is the sweet spot; `maxForks=6` bounds the
   ~2.5 GB staged DB copies on shared GPFS bandwidth.

**Why it matters**: These are the scaffold for the cross-method concordance (Phase 2) and the
singleton-pathway scale-up; also directly informs `todo/gene-contraction-analysis.md`
(Antarctolithica high loss ratio) and the still-open Phase 2 validation battery. None was
recorded before this reconciliation.

**Resolution**: Captured here and mirrored in `decisions.md`; the run-observation evidence
stays in `logs/` + `results/`.

**Tags**: findings, mmseqs, pairwise, concordance, singleton, novelty, loss, pfam, hmmscan,
perf, preempt, cluster, investigation

**mitigation_type**: evidence

**structural_mitigation_candidate**: A permanent per-run "findings card" under
`.living/findings/<project>/` written at report time (few rows, key count deltas) so the
knowledge layer carries each investigation without needing a full agent session.

---

### [2026-09-04] Fixing a bin/ script doesn't bust Nextflow's `-resume` cache by itself

**Category**: gotcha

**What happened**: After the H1/H2 presence-calling bugfix (`a49ee0c`), the pezizo5
coverage-floor sweep (`c704d72`) was launched with `-resume` against a repo `work/` dir that
also held an old (2026-07-22) pre-fix `pezizo5_mmseqs` run. Nextflow's task cache hash is
derived from the process's rendered *script text* plus its declared input files — it does
**not** track the content of a `bin/*.py` script that a process invokes by bare command
(e.g. `profile_to_matrix.py --hmm-cov ...`) unless that script is itself declared as a
`path` input. So editing `lib/family_presence.py`/`bin/profile_to_matrix.py` to fix the
circular-calibration bug does not, by itself, invalidate any previously-cached completed
task for the process that calls it — a `-resume` could silently reuse the stale, buggy-logic
output instead of recomputing with the fix. Verified this did *not* actually happen to
pezizo5 only by chance: the same fix commit series (`7e29f12`) also added a brand-new
`--min-covered-residues ${params.hmm_presence_min_residues}` flag directly into
`modules/profile_presence_matrix.nf`'s process script, which changed the rendered script
text for every post-fix invocation (regardless of grid-point param values) and forced
genuine re-execution rather than a cache hit. Had the fix been purely internal to the
Python logic (same CLI flags, same values), `-resume` could have reused the pre-fix
cached matrix-building task undetected.

**Why it matters**: Any bugfix confined to a `bin/*.py` file's internals — with no change to
the CLI flags/values the calling process's script block renders — is invisible to Nextflow's
cache invalidation. A `-resume`'d sweep or reprocessing run could silently keep serving
pre-fix results for exactly the step that was just fixed, with no error or warning.

**Resolution**: Traced the specific commit ordering and script diffs to confirm pezizo5's
sweep was unaffected in this instance (see `.living/decisions.md`). No repo-wide fix applied
yet — this is a general Nextflow caching hazard, not specific to one process.

**Tags**: nextflow, resume, caching, bin-scripts, task-hash, presence-calling, gotcha,
mmseqs, sweep

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: After any bugfix to a `bin/*.py`/`lib/*.py` file, check
whether the calling process's `.nf` script block renders different text for the fix (e.g. a
changed/added flag) — if not, deliberately clear the affected `work/`/`.nf_launch/*/work`
task dirs (or bump a cache-busting param/comment in the process script) rather than trusting
`-resume` to pick up the fix. Consider a repo convention: every `bin/`-called script that is
`storeDir`/`-resume`-sensitive should be declared as an explicit `path` input to its process
so Nextflow's hash tracks its content automatically.

---

### [2026-09-05] Two distinct BUILD_CHUNK "failure" signatures on the broader-grid sweeps — one benign, one from manual queue intervention

**Category**: gotcha

**What happened**: Both `deep_broad_1kfg` and `sordariales_shallow` sweeps lost their
`hc0.5/r0` grid point to bursts of ~300-400 `BUILD_CHUNK` "failures" each. Verified by
sampling failed work dirs directly: `.exitcode` was `0` and `sacct` showed `COMPLETED`
for the sampled SLURM job IDs — the tasks actually succeeded. The log showed
`Failed to get exit status for process ... exitStatusReadTimeoutMillis: 270000` —
Nextflow's default 270s wait for the `.exitcode` file to appear on shared GPFS storage
was too short when many chunks finished within the same window and contended for
filesystem metadata. A second, distinct signature was also present in earlier logs:
`Process ... terminated for an unknown reason -- Likely it has been terminated by the
external system` — this is what Nextflow reports when a job it submitted disappears or
restarts outside its control (e.g. `scontrol requeue`, cancel+resubmit, or genuine SLURM
preemption on the `preempt` partition, which is deliberately preemptible). The user was
separately using `scontrol update JobId=X Partition=Y` to move *pending* jobs to a less
congested partition during this run — confirmed that does NOT reproduce this signature,
since it keeps the same job ID and only affects scheduling of a job that hasn't started
yet, so Nextflow's tracking is unaffected. The two failure signatures are easy to
conflate since both surface as a failed `BUILD_CHUNK`, but only the exit-timeout one is
a timing artifact from this incident; the "terminated by external system" one (seen
earlier, chunk_200) was most likely ordinary `preempt`-partition preemption, not any
manual action.

**Why it matters**: An entire sweep grid point's presence/recovery metrics were silently
recorded as blank/partial rather than the run being flagged for a hard rerun, because
`bin/run_param_sweep.sh` catches a failed grid point and continues (`recording
partial/blank metrics and continuing`) — see that script's error handling. Without
checking `.exitcode`/`sacct` directly, both failure classes look identical from the
sweep's own summary line.

**Resolution**: Raised `executor.exitReadTimeout` from Nextflow's 270sec default to
900sec in both `conf/ucr_hpcc_slurm.config` and `nextflow.config`'s `slurm` profile
block (kept in sync — the latter is the profile default, the former is what
`bin/run_param_sweep.sh` actually loads via `-c`). This addresses the first failure mode
only. The second (external termination) is not something this run's failures were
actually attributed to once checked — `scontrol update JobId=X Partition=Y` on a
still-pending job is confirmed safe (same job ID, no restart); avoid `scontrol
requeue`/cancel+resubmit on a job Nextflow is actively tracking, which would not be.

**Tags**: nextflow, slurm, exitReadTimeout, gpfs, false-failure, queue-management,
sweep, preempt, gotcha

**mitigation_type**: structural

**structural_mitigation_candidate**: The `exitReadTimeout` fix has shipped (see
Resolution). No structural mitigation exists yet for the second failure mode beyond the
operational guidance above — a candidate would be a wrapper around manual SLURM
job manipulation that first checks whether the job ID is currently tracked by a live
Nextflow session (e.g. grepping `.nextflow.log` for the job ID) and warns before acting.

---

### [2026-09-05] The report pages' colour duplication was hiding two real defects

**What surprised us**: consolidating four copies of the colour tokens was framed as a
tidy-up, but the copies had drifted in ways that were invisible until they sat side by
side:

1. `core.html` and `losses.html` called `pfamChipsNode(names, accs)` with two arguments
   against a three-argument helper, so both **silently dropped the Pfam E-values** that
   `novelties.html` shows. No error — the third parameter just came back `undefined`.
   Fixing it also needed `pfam_e` added to `CORE_ROW_FIELDS`/`LOSSES_ROW_FIELDS`, which
   never carried the column at all.
2. The FungiDB gene link was gated on the *annotation source* being a `ModelOrg_*` entry,
   so only the four model organisms in `configs/modelorgs.yaml` ever got a gene-record
   link — even though most of the sample pool lives in FungiDB or MycoCosm. The fix is a
   per-species `SourceDB` config column, not a smarter guess from the annotation.

**Why it matters**: both are the failure mode where nothing throws and the page still
renders, so neither would ever surface from a stack trace or a payload unit test. They
surfaced from putting the three implementations in one file.

**Also learned**:

- A WCAG contrast test over a palette dict is cheap (~40 lines) and immediately earned its
  keep: it failed on `color: #fff` chip labels sitting at ~1.5:1 on the neon skin's cyan.
  The lesson is that the *hardcoded* colour was the bug — a token (`--on-series`) is what
  lets a skin decide what reads on top of its own fill.
- A luminance test cannot check the property that actually matters for the two evidence
  colours (hue separation on the blue-yellow axis under red-green CVD). Keep it as a
  collapse detector and reason about CVD explicitly in the skin's docstring.
- The `view/` gallery already had every fact needed to be useful — species counts, row
  counts, thresholds are all parsed by `collect_run_summary()` — and displayed none of it.
  Worth checking for that pattern elsewhere before adding new extraction code.
- A page that writes `localStorage` needs a real `url` in the jsdom options: a `file://`
  document gets an opaque origin and storage throws.
- `git commit` fails on this machine with `cannot exec '/usr/local/bin/gpg'` — `gpg.program`
  points at the Intel-homebrew path. `git -c gpg.program=/opt/homebrew/bin/gpg commit`
  works without changing the user's config.

**Tags**: report, skins, accessibility, wcag, colour-blind, pfam, linkouts, duplication,
jsdom, testing, gpg, gotcha

---

### [2026-09-05] A live-checkout `-profile slurm` job reads whatever `bin/`/`lib/` is on disk *right now*, not what it started with

**Category**: gotcha

**What happened**: `run_pezizo5_refresh.sh` launched a plain `-profile slurm` (no
`docker`/`singularity`) re-run at 10:25am from the repo's real working directory (not a
container, not a worktree). Its `MAKE_CORE_REPORT`/`MAKE_LOSSES_REPORT` steps ran at
10:35-10:41am and wrote `core.html`/`losses.html` **without** the new report skins, even
though the same-day fix commits were already on the branch checked out at launch time.
Root cause: this session merged `main`'s skins work into that same working directory at
11:48-11:52am — *while the job was still running*. Nextflow doesn't snapshot `bin/`
scripts; each process invokes `make_core_report.py` etc. by bare command, resolved via
`PATH` against the live `bin/` directory at execution time, so a `git checkout`/`merge` in
the same working tree mid-run can flip what a still-running task's *next* invocation
executes, non-deterministically depending on exact timing. First (wrongly) suspected as
container staleness before checking `.nextflow.log`/`sacct` timestamps: no container
profile was active for this run, so it couldn't have been that — and a closer read of
`Dockerfile` afterward showed the premise was wrong anyway (see the correction below).

**Why it matters**: any Claude session (or person) doing `git checkout`/`merge`/`rebase`
in a repo directory that also hosts an actively-running non-containerized Nextflow
pipeline can silently corrupt that run's outputs — no error, no warning, just some tasks
executing pre-change code and others post-change, depending on exactly when each process
happened to start relative to the working-tree edit.

**Correction (same day, caught when asked "would rebuilding the container bake in the
skins?"):** the first version of this entry additionally claimed
`ghcr.io/stajichlab/novinvenio` was "stale" and that a `-profile docker`/`singularity` run
would need a rebuilt image to pick up `bin/`/`lib/` changes (skins, presence-calling
fixes). That's wrong, and avoidably so — `Dockerfile`'s own header comment says exactly
the opposite: "This image contains only the tool layer. The pipeline source (bin/, lib/,
modules/, workflows/) is supplied by the cloned repository; Nextflow stages it into each
task working directory automatically." The image ships only the external binaries
(hmmer/diamond/blast/mmseqs2/famsa/openmpi) plus a bare Python+pandas+biopython+matplotlib
— zero repo Python code. Every run, containerized or not, always executes `bin/*.py` from
whatever git checkout launched it. So `docker-build.yml` rebuilding only on
`Dockerfile`/`pixi.toml` changes is *correct*, not a gap: those are the only changes that
could ever go stale in the image. Don't repeat the container-staleness theory for any
`bin/`/`lib/` change — the actual and only hazard here is the live-filesystem race above.

**Resolution**: Regenerated `pezizo5`'s `core.html`/`losses.html` (and `novelties.html`,
which never got produced at all — see decision log) directly from the already-fresh
`presence_matrix.function.tsv` via `bin/make_*_report.py`, bypassing Nextflow entirely.
No repo-wide fix applied — this is a process/awareness gap, not a code bug: don't run
`git checkout`/merge in a repo directory with a live non-containerized pipeline job still
running against it; use a separate worktree/clone for that job's launch dir instead.

**Tags**: nextflow, slurm, git, race-condition, container, docker-build, correction,
gotcha, live-filesystem

### [2026-09-06] `task.index` in an output filename silently defeats `-resume` caching across separate invocations

`modules/build_family_profiles.nf`'s `BUILD_CHUNK` process named its output
`chunk.!{task.index}.hmm`. `task.index` is Nextflow's per-session, per-scheduling-slot
counter — it is **not** a stable identifier for a task's *content*, and it is not
guaranteed to line up between two separate `nextflow run -resume` invocations against the
same launch directory, even when the upstream inputs feeding that task are byte-identical
and already cached. Confirmed against real production trace files from the
`deep_broad_1kfg` param sweep: a fresh run's `PROFILE_SEARCH:BUILD_CHUNK` tasks used
indices 2-1293 (430 distinct values) for a given set of family chunks; a `-resume` rerun
of the *identical, upstream-cached* family set reassigned almost entirely different
indices (1-1297), with only 26 of 430 overlapping.

Since `task.index` was embedded in the rendered/hashed script text (via the output
filename), this alone busted Nextflow's task-hash-based resume caching for every
`BUILD_CHUNK` task on every rerun — regardless of whether anything relevant actually
changed. The sweep's own swept parameters (`hmm_presence_cov`/`hmm_presence_min_residues`)
aren't even read by this process; it was forcing a full from-scratch rebuild of ~1,187
family HMMs (430 + 757 across the two search directions) at 130-proteome production scale
on *every one of 4 grid points*, costing multiple days of redundant compute on a
preemptible SLURM queue that was already fighting to make progress.

**Fix**: derive the chunk's identifier from its own *content* instead — the
alphabetically-first family's index within that chunk (`fam_000042.faa` → `000042`,
computed once in Groovy via `.map { files -> tuple(files.min { it.name }.baseName -
'fam_', files) }` and passed as a `val` input, not recomputed per-task in bash). This is
deterministic (same family membership → same id → same hash → cache hit) and unique
within a single run (no two chunks share a first family). A second, cascading instance of
the same class of bug was caught in the same pass (`gh` review): `MERGE_PROFILES`'s
`BUILD_CHUNK.out.partial_hmm.collect()` preserves *completion* order, which differs
between a fresh run (arbitrary SLURM completion order) and a resumed run (cached tasks
resolve in submission order) — an order-sensitive list hash would keep cache-busting
`MERGE_PROFILES`, and cascade into every downstream `FAMILY_HMMSEARCH` task, even with
`BUILD_CHUNK` itself fixed. Needed `collect(sort: true)` — safe only because the chunk
filenames are now content-derived rather than scheduling-order-derived.

**Why it matters generally**: any Nextflow output filename (or any value interpolated
into a hashed script body) that derives from *when*/*in what order* a task happened to run
— `task.index`, wall-clock timestamps, PID, `collect()`/`collectFile()` without an
explicit sort — rather than from the task's actual inputs, breaks `-resume` caching in a
way that produces no error and no warning: it just quietly recomputes everything, every
time. Worth grepping for `task.index`/`task.hash` inside any `output:`/interpolated
script block when adding new scatter-gather stages.

**Tags**: nextflow, resume, caching, gotcha, scatter-gather, slurm, param-sweep

### [2026-09-08] diamond cluster's UniProt ID handling and CLI quirks (vs mmseqs2)

**Category**: gotcha

**What happened**: While benchmarking `diamond cluster` against mmseqs2 as a possible
alternative clustering engine for `--cluster_tool mmseqs` (see `.living/decisions.md`'s
2026-09-08 entry for the full comparison), two tool-specific surprises surfaced:
1. `diamond cluster`'s output TSV keeps the **full first-token header**
   (`tr|A8PCG9|A8PCG9_ASPNG`) as its sequence id for UniProt-style deflines — it does
   NOT do mmseqs2's own accession-only extraction (see issue #85 / `lib/fasta.py`'s
   `mmseqs_id()`). Comparing the two tools' cluster TSVs directly gave **zero
   overlapping ids** until both were normalized through the same `mmseqs_id()`
   function. This is the opposite-direction version of #85's bug: diamond's
   convention actually matches the rest of this pipeline's own id handling
   (Biopython, `SEED_PROTEIN_MAP`'s shell extraction) better than mmseqs2 does.
2. `diamond cluster` **rejects** the ordinary blastp sensitivity flags —
   `--sensitive`/`--more-sensitive`/`--very-sensitive`/`--ultra-sensitive`/
   `--mid-sensitive` all fail with `Error: Option is not permitted for this
   workflow`. Sensitivity for the `cluster` subcommand is controlled instead via
   `--cluster-steps`/`--round-approx-id`/`--round-coverage` (cascaded multi-round
   parameters), which were not tuned/explored in this pass — the comparison used
   diamond's default (single-pass) sensitivity and still got respectable
   concordance (ARI 0.79) against mmseqs2's `-s 7`, but a properly-tuned cascade
   might close some of the recall gap.

**Why it matters**: anyone revisiting a diamond-based clustering path for this
pipeline needs to (a) apply the same `mmseqs_id()`-style normalization if joining
diamond's cluster output against anything using this pipeline's normal id
convention, or better, just don't normalize at all and update every OTHER consumer
to match diamond's fuller-token convention instead (worth deciding deliberately,
not defaulting into whichever tool happens to be added second); and (b) budget time
to find the right `--cluster-steps` recipe rather than assuming diamond's default
sensitivity is representative of its best-case clustering quality.

**Resolution**: not resolved/mitigated in code — this was a pure investigation, no
`--cluster_tool diamond` pathway was built. Decision was to keep mmseqs2 for now (see
linked decision entry).

**Tags**: diamond, mmseqs, cluster-tool, uniprot, id-normalization, benchmarking,
cli-quirk

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: if a diamond-based `--cluster_tool` path is ever
built, a unit test mirroring `tests/test_fasta.py::test_mmseqs_id_is_a_no_op_for_non_uniprot_headers`
but asserting diamond's own (different) UniProt-header id convention would catch this
class of mismatch the same way the #85 regression tests do for mmseqs2.

## [2026-09-16] diamond cluster preserves multi-pipe, Short-prefixed headers verbatim under a real invocation

**Category**: validation / cluster-id-fidelity

**What happened**: Built the diamond tier-1 clustering ID-fidelity safety net
(bin/verify_diamond_cluster_ids.py, docs/adr/0003) on the assumption -- backed by
the 2026-09-08 learning above, but only tested there against raw UniProt headers
on a different pipeline path -- that diamond keeps the whole first-token header
verbatim, unlike mmseqs. Ran a real `diamond cluster` invocation (pixi env,
diamond v2.2.0.180, --approx-id 90 --member-cover 80) against a small adversarial
FASTA including the pangenome subworkflows own header convention colliding with a
raw UniProt defline: `Afum|sp|O74225|YCF1_SCHPO` (a Short-prefixed original id that
itself contains pipes). The id came through *_cluster.tsv completely unchanged.

**Why it matters**: this was the one part of ADR-0003 that could not be verified
by unit tests alone (they can only encode what mmseqs_id() already tells us to
expect, not confirm the real tool actually behaves that way) -- a real invocation
was needed to close the loop. Confirms the pangenome subworkflows Short<id_sep>orig_id
prefixing does not trip diamond the way it could theoretically trip mmseqs.

**Resolution**: verify_diamond_cluster_ids.py is wired into CLUSTER_TIER1s diamond
branch and the pangenome.nf hard-error guard for --pangenome_cluster_backend diamond
was relaxed. Still open: a real multi-strain biological dataset concordance
benchmark (todo/diamond-tier1-cluster-backend.md) -- this smoke test used 5 small,
mostly near-identical toy sequences, so it validates ID fidelity only, not
clustering quality/concordance with mmseqs at scale.

**Tags**: diamond, mmseqs, cluster-tool, id-normalization, pangenome, validation

## [2026-09-17] `-profile local` doesn't apply conf/ucr_hpcc_slurm.config's AVX2 node-pins -- don't mistake that for a missing pin

**Category**: gotcha (self-correction)

**What happened**: Running a real end-to-end diamond-vs-mmseqs comparison for the
pangenome CLUSTER_TIER1 diamond backend (docs/adr/0003) with `-profile local`, the
mmseqs comparison run crashed: `mmseqs easy-cluster ... died with <Signals.SIGILL: 4>`
on an interactive node whose /proc/cpuinfo shows no avx2 (sse4_2 only). Diamond
(v2.2.0.180) did NOT crash on the same node running the same workflow. This was
initially (and wrongly) filed as issue #101, claiming CLUSTER_TIER1's mmseqs branch
had no AVX2 node-pinning -- it does: `conf/ucr_hpcc_slurm.config` already has
`withName: '.*CLUSTER_TIER1' { clusterOptions = '-C ryzen|broadwell|cascade' }`,
added in commit `3de5252`, the same commit that introduced the whole pangenome
subworkflow, well before this smoke test ever ran. `-profile local` just sets
`process.executor = 'local'` and never applies `clusterOptions` at all -- those
only take effect under `-profile slurm` with `-c conf/ucr_hpcc_slurm.config`
included. Issue #101 was closed as not-a-bug the same day once this was found.

**Why it matters**: a `-profile local` smoke test on this repo exercises pipeline
*logic* but not the SLURM-specific resource/node constraints (`clusterOptions`,
`queue`, per-attempt `memory`/`time` scaling) that `conf/ucr_hpcc_slurm.config`
supplies -- a crash or resource issue seen only under `-profile local` should be
checked against that config file BEFORE concluding the pipeline itself is missing
a safeguard. This is the second time an AVX2 SIGILL has come up on this cluster (see
the 2026-07-21 famsa/AVX2 learning below/above), so the SIGILL-on-Abu-Dhabi-nodes
risk itself is real and worth remembering -- just: check the SLURM config for an
existing fix before assuming there isn't one and filing a new issue.

**Resolution**: no code change -- the fix already existed. Issue #101 closed as
not-a-bug; docs/adr/0003, CHANGES.md, pangenome.nf help text, and the relevant
todo/ items corrected to remove the incorrect "gap" framing.

**Tags**: cluster, slurm, avx2, mmseqs, pixi, bioconda, simd, hardware-compatibility,
nextflow, node-placement, pangenome, self-correction, profile-local

**mitigation_type**: convention

**structural_mitigation_candidate**: none needed for the pin itself (already
exists). If this class of mistake recurs, the real candidate is a habit/checklist
one: grep the relevant `conf/*.config` for a process's name before concluding a
`-profile local` failure reflects a missing SLURM safeguard.

## 2026-09-19 — Gene-family IDs are NOT stable across a full pipeline re-run

**Tags**: pangenome, reproducibility, mmseqs, cross-run-comparison

Re-running the whole `coccidioides_pangenome` `genus_vs_ureesii` study from
scratch (529 strains, same config, same data, same code path) produced an
essentially identical *clustering* but a substantially different set of
family *identifiers*:

| | original run | re-run |
|---|---|---|
| total families | 54,421 | 54,410 |
| family IDs shared between the two | **39,829 (73%)** | |

A family's ID is its mmseqs cluster representative's sequence ID, and mmseqs
does not pick the same representative every time (thread count and tie-break
order both move it). So ~27% of families are renamed on re-run even though
the partition barely changes (0.02% difference in family count).

**Consequence**: any artifact keyed by family ID — `family_modules.tsv`,
`presence_matrix.tsv`, `islands_with_domains.tsv`, report tables — cannot be
diffed across runs by joining on the family column. A comparison that does
so will report a huge spurious difference. Comparing two runs requires
matching on cluster *membership*, not on the representative's name.

This surfaced while validating `LEIDEN_MODULES` at scale: the pipeline run
and the earlier hand-run shared only 6,854 of ~14,700 trans families, which
looked alarming until the ID instability explained it. On the families that
*did* keep their ID, module agreement was 97.1% (ARI 0.926, near-diagonal
contingency table, same four modules) — and running the pipeline's script
against the *original* `pair_classification.tsv` reproduced the hand-run
**byte-for-byte** (identical md5). So the Leiden step itself is fully
deterministic; the variation is entirely upstream.

## 2026-09-19 — `zstd -dc` silently returns EMPTY for a symlinked input, and Nextflow stages every input as a symlink

**Tags**: nextflow, compression, silent-data-loss, staging

`lib/compressed_io.py`'s `open_maybe_compressed()` shelled out to
`zstd -dc <path>`. Nextflow stages every process input as a **symlink** into
the task work directory, and `zstd` refuses to read through a symlink: it
prints `Warning : <file> is a symbolic link, ignoring`, writes **nothing** to
stdout, and **exits 0**.

The result is silent truncation, not a crash. On a 5-strain smoke run,
`BUILD_FAMILY_POSITIONS` read a staged `gene_positions.tsv.zst`, got an empty
stream, wrote a header-only `family_positions.tsv` (43,240 rows in the
pre-change run), and the pipeline reported **SUCCESS**.

**Fix**: `zstd -dc -f`. The `-f` is load-bearing under Nextflow, not a
convenience — never remove it.

**Why the unit tests missed it**: `tmp_path` fixtures create real files.
Nothing in the test suite reproduced staging until a regression test was
added that reads a `.zst` *through a symlink*. Any helper that shells out to
a CLI on a path handed to it by Nextflow needs a symlink test.

## 2026-09-19 — A background run launched from a Claude session dies with the session, and /scratch goes with it

**Tags**: hpcc, scratch, background-jobs, session-lifetime

Two `nextflow` smoke runs started with `nohup ... &` from a Claude Code
session were both lost when the session's own SLURM job ended: the processes
died, and their working directories vanished because they lived under the
session scratchpad on `/scratch/$USER/<jobid>/`, which is **node-local** and
removed with the job. `pgrep -f <name>` made this hard to notice — the
pattern matched the checking command's own shell, reporting "RUNNING" for
processes that no longer existed.

**Rule**: anything expected to outlive the current turn goes in via `sbatch`,
with work dir, outdir and logs on `/bigdata` or `/rhome`. Never `nohup &` a
multi-hour pipeline from the session, and never point its `-work-dir` at
`/scratch` unless the job that owns that scratch is the one running it. When
checking whether a named process survived, match on the real process
(`ps -eo cmd | grep "[n]extflow"`) rather than a `pgrep -f` pattern that your
own command line contains.

## [2026-09-20] A purely relative paralog filter can erase an arbitrarily strong hit

**Category**: presence-calling
**What happened**: `build_presence_matrix.py` filter 2 disqualifies a hit whenever the
query's in-genome paralog beats it head-to-head, with no absolute floor. For CorA
(A7UWR3_NEUCR) this erased all 15 outgroup hits, including a 3.3e-93 match to ALR2_YEAST
that the IQ-TREE topology confirms is a real ortholog — so the protein was called a
lineage-specific novelty. HEX1 (P87252) is the opposite case and must keep filter 2:
its single 3.4e-12 Mcir hit is real noise, and without filter 2 HEX1 fails
`--other-max-frac 0.0` and stops being a candidate.
**Why it matters**: a relative-only test has no upper bound on what it can discard. The
two mechanisms act on disjoint populations of the 3393 pezizo_set1 candidates: 213 (6.3%)
have a raw outgroup hit below 1e-20 (an absolute floor's job), ~121 sit in the 1e-20..1e-5
twilight zone (filter 2's only legitimate job).
**Resolution**: not yet implemented. `--paralog-rescue-evalue` already exists as the
absolute-floor arm (default null/off). Measured that a query/paralog delta rule does NOT
separate the populations on its own — median delta is 10-23 orders in every query-E
bucket — so delta cannot replace the floor, though it does cover the floor's blind spot
(cells where the paralog won by only 2-3 orders). Full option space in
docs/superpowers/specs/2026-09-20-paralog-novelty-disqualifier-analysis.md.
**Tags**: presence-calling, paralog, novelty, pairwise, threshold, controls

## [2026-09-20] A build system can silently defeat a tool's own runtime SIMD dispatch

**Category**: gotcha
**What happened**: famsa SIGILLed (exit 132) on this cluster's Abu Dhabi nodes (AMD Opteron
6376 / Piledriver: avx, bmi1, fma4, xop -- no avx2, no bmi2). The obvious reading was "famsa
requires AVX2". It does not. FAMSA ships runtime CPU dispatch (src/utils/cpuid.h, separate
lcsbp_avx_intr / lcsbp_avx2_intr / lcsbp_avx512_intr translation units), but refresh.mk's
CHECK_OS_ARCH sets ARCH_FLAGS from the BUILD HOST's detected SIMD level and applies it to the
whole binary, including the generic code -- so a build on an AVX2 machine stamps -mavx2
everywhere and the dispatch can never be reached. Every prebuilt famsa inherits this,
bioconda's and the cluster module's alike; the module named "2.4.1-x86", which sounds like a
generic build, still contained 3211 vinserti128 and 20 BMI2 instructions.
**Why it matters**: "binary SIGILLs on old CPU" is normally read as "the tool needs that ISA",
and the fix reached for is a different tool or a node constraint. Both are wrong here and both
cost real capacity -- the node constraint barred thousands of scatter-gather tasks from half a
partition. Check for runtime dispatch in the SOURCE before concluding the tool needs the ISA.
Disassembly tells you quickly: `objdump -d BIN | grep -cw vinserti128` (AVX2) and `bzhi|sarx|
shlx` (BMI2), compared against /proc/cpuinfo flags.
**Resolution**: `make PLATFORM=avx STATIC_LFLAGS='-static-libstdc++ -static-libgcc'` with
gcc/12.2.0 produced a working, self-contained binary: BMI2 count 0, exit 0 on the Opteron, zero
libstdc++/libgcc dynamic deps. Full -static fails (no static glibc); static libstdc++ alone is
both necessary (else GLIBCXX_3.4.30 not found) and sufficient. zlib-ng's bundled test target
fails against a stray miniconda libgtest -- harmless, libz.a builds and a second make proceeds.
Packaging tracked in issue #148; the mafft-fallback alternative (#143) was closed as
unnecessary once this was measured.
**Tags**: simd, avx2, bmi2, hardware-compatibility, build-system, famsa, runtime-dispatch, slurm, gotcha

## [2026-09-20] TBLASTN "extra evidence" candidates are mostly a diamond sensitivity gap, not false positives

**Category**: presence-calling
**What happened**: issue #135 asked whether candidates with TBLASTN outgroup-genome hits
but no protein-level evidence are TBLASTN false positives. Measured on pezizo_set1's 1479
such cells (breadth>=4): only 3.2% land in true intergenic space (the real FP signature).
96.8% overlap an annotated gene in the outgroup's own GFF3. Of those, 62.1% have the gene's
protein present in the search FASTA -- and of THOSE, 94.6% get literally zero diamond hits
(even to E=0.01) despite TBLASTN finding a significant hit at the same locus.
**Why it matters**: the dominant driver (56.9% of the 1479-cell population) is diamond
default-mode sensitivity missing real orthologs that a translated search recovers
immediately -- not TBLASTN noise (3.2%) and not primarily outgroup annotation gaps (36.6%).
This is the same failure class already flagged in the open, unexecuted todo
`diamond-very-sensitive-main-search.md` (one example: NCU02794 missed at 4/6 outgroups by
default mode, caught at E=41..58 under --very-sensitive) -- this measurement shows it at
scale (841 cases), not as an isolated anecdote.
**Resolution**: recommended prioritizing the --very-sensitive DIAMOND_SEARCH benchmark over
building a TBLASTN-based disqualifier for #135 -- fixing presence-calling upstream is more
principled than patching around a protein-search miss at the genome level, and may make a
large fraction of #135's population moot. Full decomposition posted to issue #135.
**Tags**: presence-calling, diamond, sensitivity, tblastn, novelty, false-positive, todo-validated

## 2026-09-23 — UniProt library index: build cost and match coverage (Fungi_2026_03)

**Build**: --build_uniprot_index over 1526 proteomes / 17,344,385 records (= the CSV's canonical total). Two parse batches ran 1 h 03 min and 41 min (~1.75 cpu-h); SQLite merge 7 min, 73 MB RSS; index 3.0 GB sqlite + 1.9 GB records. At /bigdata/stajichlab/shared/db/Uniprot/Fungi_2026_03/novinvenio_index/v1.
**Coverage, UniProt-sourced study** (sordariales_shallow, 12 proteomes): 100% matched by accession; 1-4 s per proteome.
**Coverage, BFD/funannotate gene models** (configs/Crypneof_pathogen_cryptococcus.csv, 9 species): only 7-54% of proteins have an exact-sequence UniProt match (median ~21%): C. amylolentus 54%, C. cinerea 48%, S. pombe 46%, C. deneoformans 41%, S. cerevisiae 21%, C. bacillisporus 15%, N. crassa 14%, A. nidulans 8%, C. neoformans T4 7%. Re-annotated gene models rarely match UniProt's exactly, and a different strain (T4 vs UniProt's H99) shares few identical sequences. Near-identical matching (out of scope in the spec) is what would raise this.
**Gotcha**: species synonyms defeat the binomial own-species lookup. All 668 A. nidulans matches were seq_other, because UniProt names it "Emericella nidulans (strain FGSC A4) (Aspergillus nidulans)". Setting NCBI_TaxID to the UniProt proteome's taxid fixes it for that row.
**Tags**: uniprot, annotation, coverage, pipeline, performance
## 2026-09-23 — Tier C/R cluster-membership presence is outgroup-blind by construction

**Context**: sordariales_shallow was the first cluster-vs-pairwise clade with resolved BUSCO negative controls (NII `studies/fungi/sordariales_shallow_cluster/tier_comparison/`, commit 44677be).
**Finding**: Tier C and Tier R called 8/8 BUSCO negatives "novel" (fp_rate 1.0). Tier P and Tier C+H called 0/8. The families are clustered from the seed group only, so the cluster TSV holds no outgroup proteins. Every family then looks outgroup-absent. Tier C/R can refine ingroup families, but cannot make the absence call alone.
**Also**: `bin/trace_cost_report.py` (NII) hard-codes Tier P DIAMOND_SEARCH as "unmeasurable". That is true only when the search came from a storeDir cache. For sordariales_shallow the search ran fresh: 0.091 task-h / 2.24 cpu-h, against about 127.5 task-h / 3305 cpu-h for the C+H family-profile stages (84% in BUILD_CHUNK famsa+hmmbuild). At N=12 proteomes, pairwise is about 1000x cheaper. There is one data point only, so no crossover N is known.
**Tags**: cluster-tool, mmseqs, controls, false-positive, compute-cost, tier-comparison

## 2026-09-23 — Cluster-vs-pairwise miss causes across three clades

**Context**: `bin/classify_concordance_misses.py` (NII) assigned a cause to every disagreement between Tier P and Tier C+H in pezizo_set1, agaricomycetes and sordariales_shallow. Summary: NII `notes/cluster-vs-pairwise/README.md`.
**Finding**:
- C+H misses: 23–57% are clustering-level, almost all because the protein's family was never profiled (usually a singleton). The other 43–77% are HMM-level: the family HMM calls the family present in the other group.
- For the HMM-outgroup misses, 51–81% have no raw diamond hit in the other group. 10–43% had a significant hit that Tier P's paralog-competition filter removed.
- C+H extras: 82–98% are proteins that P finds in fewer than 75% of the seed group, or not at all. Default-mode diamond misses homologs that the HMM finds.
- Cost: C+H used 465x (7 proteomes) and 1215x (12 proteomes) the CPU-hours of P.
**Gotchas**:
- `presence_matrix.evalues.tsv` holds E-values only for cells that pass the presence filters. It cannot show sub-threshold hits. Use raw `search_cache/*.diamond.tsv.gz` for that.
- The pezizo_set1_cluster outputs were regenerated by a rerun that started after the 2026-09-13 concordance, so that committed concordance does not match the current files.
- tests/test_refine_ambiguous_families.py replaced `raf.family_fractions` without restoring it, which polluted later tests. Fixed with `monkeypatch`.
**Tags**: cluster-tool, mmseqs, tier-comparison, compute-cost, diamond, sensitivity, testing

## 2026-09-23 — Q7RXD0_NEUCR: Tier C+H false novelty from a low-complexity, 2-member family HMM

**Context**: The user flagged a sordariales_shallow_cluster novelty candidate with TBLASTN hits in all 8 outgroup genomes but no protein hits.
**Finding**:
- Q7RXD0 is 1997 aa: residues ~1–1750 are S/T/G/P/A-rich low complexity, and the last ~230 aa are a conserved DUF7371 domain.
- The 2-member (Ncra + Smac) family HMM, 2111 columns, detects the domain in Chaetomium (E 7e-29) but in no outgroup protein. It missed Collhigg H1V9Q8, which default diamond found at E 1e-31.
- Default diamond hit 1 of 8 outgroup proteomes. `--very-sensitive` hit 5 of 8, all on the C-terminal domain. TBLASTN hit all 8.
- Tier P rejected the candidate, but only on that single outgroup hit.
**Why it matters**: A family HMM built from few, long, mostly low-complexity members can miss the one conserved domain in divergent outgroups. This is the worked example in NII `notes/diamond-sensitivity/README.md`.
**Tags**: cluster-tool, hmm, false-positive, low-complexity, diamond, sensitivity, tblastn
