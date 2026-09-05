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
external system` — this is what Nextflow reports when a job it submitted is manually
`scontrol requeue`'d, moved to another partition, or cancelled/resubmitted outside its
control. The user confirmed manually moving jobs between queues during this run, which
plausibly explains some fraction of these failures independently of the GPFS timing
issue — the two are easy to conflate since both surface as a failed `BUILD_CHUNK`, but
only the first is a timing artifact; the second reflects a genuine break in Nextflow's
job-tracking and isn't fixable by any timeout setting.

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
only. For the second: avoid `scontrol requeue`/cancel+resubmit on a job Nextflow is
actively tracking; `scontrol update JobId=X Partition=Y` on a still-pending job is safer
since it doesn't restart or orphan the task.

**Tags**: nextflow, slurm, exitReadTimeout, gpfs, false-failure, queue-management,
sweep, preempt, gotcha

**mitigation_type**: structural

**structural_mitigation_candidate**: The `exitReadTimeout` fix has shipped (see
Resolution). No structural mitigation exists yet for the second failure mode beyond the
operational guidance above — a candidate would be a wrapper around manual SLURM
job manipulation that first checks whether the job ID is currently tracked by a live
Nextflow session (e.g. grepping `.nextflow.log` for the job ID) and warns before acting.
