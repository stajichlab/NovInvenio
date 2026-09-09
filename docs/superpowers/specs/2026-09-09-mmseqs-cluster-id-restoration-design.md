# mmseqs cluster.tsv ID restoration (design)

Date: 2026-09-09
Status: approved (rev. 2, post-review)

Revision note: rev. 1 was reviewed against the actual codebase and real
on-disk data (review persona: Fable). The review found rev. 1's "zero
consumers normalize" claim was wrong (two files already compensate for this
bug internally, and would *break* once the root cause is fixed unless also
changed), its `summarize_tblastn.py` claim understated a real
scientific-correctness defect as merely cosmetic, and its algorithm's
prefix-matching rule was too narrow for mmseqs's actual behavior. This
revision incorporates all required fixes from that review.

## Problem

`mmseqs easy-cluster` silently reduces certain FASTA header conventions to
a shorter form in its `*_cluster.tsv` output (rep/member columns), while its
own `*_rep_seq.fasta`/`*_all_seqs.fasta` outputs from the *same run* keep the
full header token. Every other artifact in this pipeline — `candidates.txt`,
`presence_matrix.tsv`, `novelties.<Short>.tsv`, TBLASTN's `qseqid`, the
report payload — treats the full header token as `protein_id`. Any code that
joins a `*_cluster.tsv` row against `protein_id`-keyed data therefore
silently matches nothing, for every UniProt-sourced study (confirmed:
mmseqs recognizes 13 header-prefix conventions this way — `sp|`, `tr|`,
`gb|`, `ref|`, `pdb|`, `bbs|`, `lcl|`, `pir||`, `prf||`, `gnl|`, `pat|`,
`gi|`, `cl|` — not just the two, `sp|`/`tr|`, this repo's existing partial
workaround handles; empirically re-derived against this repo's own
`mmseqs 18.8cc5c` binary during review).

Confirmed empirically (a dedicated audit plus a follow-up review, this
session) across real on-disk data:

| Study | Header convention | `*_cluster.tsv` keys | Affected? |
|---|---|---|---|
| pezizo_set1 | `tr\|F9WWV8\|F9WWV8_ZYMTI ...` | `F9WWV8` (bare) | **Yes** |
| zoosporic_dikarya | `tr\|F4NRB2\|...` | bare accession | **Yes** |
| mushrooms_tremella | `tr\|A8N3V4\|...` | bare accession | **Yes** |
| yeast_filamentous | `sp\|P47149\|...` | bare accession | **Yes** |
| UHM_Akkermansia | MAG/prodigal, no recognized prefix | full header, verbatim | No — mmseqs's header-recognition never triggers on this shape |
| UHM_Koxytoca | MAG/prodigal, no recognized prefix | full header, verbatim | No |

### Confirmed consequences

- **`lib/clusters.py`'s `FamilyIndex.index_of()`** — every UniProt-sourced
  study's report ("novelties"/"core"/"losses") shows 0 rows assigned to any
  gene family, even though mmseqs found real multi-member families (verified
  on `pezizo_set1`'s live report: 656 families computed, 0/35,597 rows
  actually placed in one). Report-cosmetic, not a correctness defect.

- **`bin/make_novelties.py`'s own `family_id`/`family_size` columns** (in the
  *persisted* `novelties.<Short>.tsv`, independent of the report's
  recomputation above) — same root cause, same silent empty result.
  Report-cosmetic.

- **`bin/summarize_tblastn.py` — a real scientific-correctness defect, not a
  cosmetic one** (rev. 1 got this wrong; corrected here). The script's job is
  rep→member expansion: a TBLASTN hit against a cluster *representative*
  should be recorded for every member of that cluster too. Because
  `rep_to_members` is built bare-keyed from `cluster_tsv` while a TBLASTN hit's
  `qseqid` is the full header, `rep_to_members.get(rep, [rep])` always misses
  and falls back to `[rep]` — **the expansion never happens.** Verified
  against real `pezizo_set1` data: 656 multi-member clusters, 1,220 non-rep
  members, 539 of which belong to a cluster whose representative *does* have
  a real TBLASTN hit — all 539 are recorded as "no hit" instead. In the
  persisted `novelties.<Short>.tsv` files, that's **539 novelty rows with an
  empty `tblastn_outgroup_hits` column whose family representative actually
  hit an outgroup genome** (e.g. `tr|Q4WWQ6|Q4WWQ6_ASPFU` in
  `novelties.Afum.tsv`, whose representative `tr|J3KH30|J3KH30_COCIM` hits
  Nirr and Mcir). `pezizo_set1`'s actual pipeline run passed
  `--skip_tblastn_filter` (`workflows/summarize.nf`), so this defect did not
  wrongly *retain* any of those 539 as novelties in this particular run — but
  the `tblastn_outgroup_hits` report column is simply wrong for them, and
  **any run that does NOT pass `--skip_tblastn_filter` would wrongly retain
  all 539 as false novelties.** This is the strongest correctness argument
  for this fix, not a side note.

- **`bin/build_alignment_shards.py` — the "0 shards" observation has a
  different, unrelated root cause** (rev. 1 misattributed this; corrected
  here). `build_alignment_shards.py`'s own rep-level fallback
  (`rep_to_members.get(rep, [rep])`) still resolves to the rep itself, and
  every TBLASTN-hit rep *is* present in `candidates.txt` — so the ID bug
  alone would only lose shards for the 1,220 non-rep members, not all of
  them. The actual reason `pezizo_set1` currently produces 0 shards:
  `results/pezizo_set1/tblastn/*.tblastn.tsv` are stale, 10-column files, but
  `build_alignment_shards.py`'s `_COLS` (matching `modules/tblastn.nf`'s
  current `-outfmt 6` string) expects 13 columns — every row fails to parse.
  **This is a separate, pre-existing data-staleness issue, not something
  this fix addresses** — `pezizo_set1`'s TBLASTN output needs regenerating
  before shards will work regardless of this spec. The ID-space fix still
  matters for shards (it's what makes the 1,220 non-rep members' shards
  possible once TBLASTN is current), just not the reason today's count is
  exactly zero.

### An existing partial workaround already exists, and would break under a naive fix

`lib/fasta.py::mmseqs_id()` (issue #85) already normalizes a full header down
to mmseqs's bare form — but only for the `sp|`/`tr|` prefixes, not the other
11 mmseqs recognizes. It's used by exactly two files to work around this
same bug, **and both would break once the root cause below is fixed, unless
also changed**:

- **`bin/extract_family_seqs.py`** (`by_mmseqs_id = {mmseqs_id(rid): rec for
  rid, rec in records.items()}`, then looks up each cluster member against
  this bare-keyed dict) — once `*_cluster.tsv` holds full tokens (this fix),
  every member lookup here would miss `by_mmseqs_id`'s bare keys, hitting the
  file's own fail-loud guard (`sys.exit("cluster member ... not found")`) —
  **this fix would turn a silent bug into a hard pipeline failure for this
  script specifically, unless it's also updated.**
- **`bin/profile_to_matrix.py`** (`mapping[mmseqs_id(pid)] = short`) — same
  shape; once cluster.tsv holds full tokens, this mapping's bare keys would
  never match, and its own documented "skip defensively" behavior means the
  presence matrix would come out silently empty instead of just wrong.

Both **must** be included in this fix's touched-file list (Decision, below)
— dropping their own `mmseqs_id()` call once the source is corrected.

A repo-wide grep otherwise found **~16-22 files** that read a `*_cluster.tsv`
and build a `protein_id`-keyed map from it with no normalization at all
(distinct from the two files above, which do normalize, just not enough).
None of those ~16-22 need to change under this fix.

Separately audited and confirmed **clean, not in scope for this fix**:
TBLASTN's `sseqid` and genome/scaffold/DNA accession handling (BLAST+ does
not do mmseqs's header reparsing; a real genome's FASTA header, GFF3 seqid
column, and TBLASTN `sseqid` all agreed end to end in the case checked). HMM
builds (`hmmbuild`/`hmmsearch`, `novelty_discovery`'s family-profile path)
introduce no *new* mismatch — `hmmbuild -n "$rep"` names a profile after
whatever ID the cluster.tsv rep column already holds and passes it through
domtblout verbatim (`lib/family_presence.py`'s whitespace-based domtblout
parsing is unaffected either way), so fixing the cluster.tsv fixes the HMM
path too, with no HMM-specific change needed.

## Decision

Fix it once, immediately after mmseqs produces each `*_cluster.tsv`, before
any consumer ever reads it — not by rewriting FASTA headers before mmseqs
(would propagate a new ID discontinuity into TBLASTN's confirmed-clean query
chain), and not by asking every consumer to normalize individually (the
existing two-file partial workaround above is itself proof that approach
drifts: issue #85 already tried exactly that, for exactly two files, with an
incomplete prefix rule, and 16-22 *other* files still don't do it at all).

### 1. Exactly two insertion points, two additional consumer fixes

Only two Nextflow processes ever invoke `mmseqs easy-cluster`:

- `modules/mmseqs_cluster.nf`'s `MMSEQS_CLUSTER` — used only by
  `--cluster_tool pairwise` (`workflows/cluster.nf`'s `CLUSTER`/`LOSS_CLUSTER`,
  confirmed via `main.nf`'s final `else` branch, ~line 217-223).
- `modules/mmseqs_family_cluster.nf`'s `MMSEQS_FAMILY_CLUSTER` — used by
  **both** `--cluster_tool mmseqs` (`workflows/profile_search.nf`, via
  `main.nf` ~line 122) and `--cluster_tool novelty_discovery`
  (`workflows/novelty_discovery.nf:331` calls it directly, via `main.nf`
  ~line 153). Fixing this one module covers both `cluster_tool` values.

`modules/profile_candidate_clusters.nf`'s `PROFILE_CANDIDATE_CLUSTERS`
(confirmed, review) does **not** call mmseqs itself — it only filters an
already-produced `family_cluster_tsv` via `bin/candidate_families.py`. Once
`MMSEQS_FAMILY_CLUSTER`'s own output is corrected, this inherits the fix
automatically.

**Full touched-file list** (revised from rev. 1, which only listed the first
two):

1. `modules/mmseqs_cluster.nf` — append the corrector call.
2. `modules/mmseqs_family_cluster.nf` — append the corrector call.
3. `bin/extract_family_seqs.py` — remove its own `mmseqs_id()` normalization
   (key `by_mmseqs_id` by `rid` directly, not `mmseqs_id(rid)`); rewrite
   `tests/test_extract_family_seqs.py` (currently asserts bare-form output,
   e.g. `>O74225` — must assert the full-token form once this lands).
4. `bin/profile_to_matrix.py` — remove its own `mmseqs_id()` normalization
   (`mapping[pid] = short`, not `mapping[mmseqs_id(pid)] = short`); rewrite
   `tests/test_profile_to_matrix.py` (currently asserts bare `protein_id`
   values — must assert full-token values).
5. `lib/fasta.py::mmseqs_id()` — **repurposed, not deleted**: this is exactly
   the "given a full header, what would mmseqs report as its bare ID"
   direction the new corrector script (below) needs to build its
   bare→full map. Extend its recognized-prefix set from the current 2
   (`sp|`, `tr|`) to all 13 mmseqs recognizes, OR (preferred — future-proof
   against mmseqs adding an unlisted convention) replace the hardcoded
   prefix regex with a general rule: match against **every `|`-delimited
   field of the header**, since mmseqs's own bare-ID output is always one
   such field for every prefix convention observed during review (`sp|`,
   `tr|`, `gb|`, `ref|`, `pdb|`, `bbs|`, `lcl|`, `pir||`, `prf||`, `gnl|`,
   `pat|`, `gi|`, `cl|` all fit this shape; `lcl|NC_000001.1_prot_XP_1_1` →
   mmseqs reports the single field `NC_000001.1_prot_XP_1_1` verbatim, no
   further splitting). `restore_mmseqs_cluster_ids.py` (below) imports this
   function directly — one shared rule, not two.

### 2. The corrector script

New script, name chosen to name the actual problem it fixes (mmseqs's own
header-recognition behavior, not a generic "cluster fixer"):

**`bin/restore_mmseqs_cluster_ids.py`**

```
Usage:
  restore_mmseqs_cluster_ids.py \
      --input-fasta <the FASTA mmseqs clustered> \
      --cluster-tsv <mmseqs easy-cluster's *_cluster.tsv, corrected in place>
```

Algorithm:

1. Parse `--input-fasta` (the exact FASTA fed to `mmseqs easy-cluster` — the
   process's own input). For each record, take the full first-whitespace-
   delimited header token (Python `str.split()`, not `split(' ')` — mmseqs
   splits on any whitespace including tabs) as `full_id`, and derive
   `bare_id = mmseqs_id(full_id)` using `lib/fasta.py`'s (extended, per
   Decision §1.5) function — the *same* rule mmseqs itself applies, so this
   never drifts from mmseqs's actual behavior.
2. Build `bare_to_full: {bare_id: full_id}`. **Fail loudly** (exit 1,
   matching this repo's fail-loud convention) if two distinct `full_id`
   values collapse to the same `bare_id`, naming both headers in the error —
   a genuine accession collision in the source FASTA (plausible for
   `MMSEQS_FAMILY_CLUSTER`'s concatenated multi-proteome seed FASTA;
   `lib/fasta.py`'s own `read_fasta()` docstring already documents this class
   of cross-proteome collision as real, not hypothetical).
3. Read `--cluster-tsv` line by line (`rep\tmember` per mmseqs's format). For
   each column value `x`: if `x in bare_to_full`, replace with
   `bare_to_full[x]` (the collapsed case — the actual fix). If `x` already
   equals some `full_id` exactly (the MAG/prodigal case, or any value mmseqs
   left untouched), leave unchanged (idempotent no-op). If `x` matches
   **neither** — an mmseqs output ID this script cannot reconcile against
   the input FASTA — **fail loudly** rather than silently emit an
   unusable row.
4. Overwrite `--cluster-tsv` in place with the corrected content.

Plain, dependency-free script (stdlib + `lib/fasta.py`, matching this repo's
existing `bin/*.py` style).

**Duplicate exact-ID handling** (not new, not introduced by this fix, noted
for completeness): if the input FASTA itself has two records sharing the
same `full_id` (already a real, tolerated case — `lib/fasta.py::read_fasta()`
keeps the first and warns), the corrector passes both `cluster.tsv`
occurrences through unchanged (same full ID, no "distinct" collision to
detect) — any resulting downstream ambiguity is pre-existing and out of
scope here.

### 3. Nextflow wiring — in-process, no new process

Append the corrector as one more shell command at the end of each of the two
processes' existing `shell` block (both already use `shell:` with `!{...}`
interpolation, confirmed), run against the exact filename `mmseqs
easy-cluster` just wrote, **before** the process's `output:` block takes
effect — so `publishDir` copies out the already-corrected file. Example for
`MMSEQS_CLUSTER` (`modules/mmseqs_cluster.nf`):

```groovy
    mmseqs easy-cluster \
        !{candidates_fa} \
        !{out_prefix} \
        ${tmpdir}/mmseqs_$$ \
        --threads !{task.cpus} \
        --min-seq-id 0.3 \
        -c 0.8 \
        --cov-mode 0

    restore_mmseqs_cluster_ids.py \
        --input-fasta !{candidates_fa} \
        --cluster-tsv !{out_prefix}_cluster.tsv
```

Same pattern for `MMSEQS_FAMILY_CLUSTER` against `families_cluster.tsv`
(input FASTA: `seed_fa`). Confirmed: `bin/` is on `PATH` inside the
container (every other `bin/*.py` script is already invoked bare this way in
these same process blocks); Nextflow's default `bash -ue` semantics mean a
non-zero exit from the corrector fails the task, matching the desired
fail-loud behavior; the empty-input-FASTA early-exit branch in both
processes returns before the corrector call, so an empty-candidates run
still short-circuits cleanly. **No change to either process's `output:`
declaration, no change to any of the ~16-22 unaffected consumer files** —
confirmed via Nextflow's `path()` staging being a same-content symlink/copy,
not a header-altering transform.

### 4. Backfill for already-computed studies

Confirmed rerun-free for the `pairwise` pathway: `candidates.fa`/
`loss_candidates.fa` and the raw `*_cluster.tsv` are exactly what's already
published to `results/<study>/`:

```bash
bin/restore_mmseqs_cluster_ids.py \
    --input-fasta results/pezizo_set1/candidates.fa \
    --cluster-tsv results/pezizo_set1/clusters/clusters_cluster.tsv
```

**Not yet confirmed for `mmseqs`/`novelty_discovery`-pathway studies:**
`MMSEQS_FAMILY_CLUSTER`'s input, `seed_all.faa`, is built via a Nextflow
`collectFile` channel operation (`workflows/novelty_discovery.nf`) — whether
it (or an equivalent) is actually published to a `results/`-visible path for
any already-run study needs checking per-study before promising a
rerun-free backfill there; do this check as the first step of implementing
this section, not assumed.

## Documentation

The strongest motivating "why" is the TBLASTN correctness defect (Problem
section), not the report-cosmetic family-grouping loss — lead with that.

- **`METHOD_DESCRIPTION.md`** (repo root) — add a paragraph after each
  pathway's clustering-step description (`pairwise`'s "clustered with
  `mmseqs easy-cluster`" and `mmseqs`'s "`mmseqs cluster` ... clusters the
  whole seed group") explaining: mmseqs recognizes 13 header conventions and
  reports only a collapsed field as the ID in `*_cluster.tsv` specifically
  (not in `*_rep_seq.fasta`); without correction this silently defeats
  rep→member TBLASTN-hit expansion (`summarize_tblastn.py`) as well as
  family grouping and alignment shards; `restore_mmseqs_cluster_ids.py` runs
  immediately after every `mmseqs easy-cluster` invocation to fix this
  before anything else reads the file.
- **`README.md`**'s "Two analysis approaches: pairwise vs family-profile
  clustering" section — one sentence pointing at `METHOD_DESCRIPTION.md`'s
  fuller explanation.
- **`bin/restore_mmseqs_cluster_ids.py`**'s own module docstring — the full
  technical explanation, following this repo's existing convention
  (matches `bin/build_alignment_shards.py`'s and `lib/fasta.py`'s own
  docstring style).
- **`modules/mmseqs_cluster.nf`** and **`modules/mmseqs_family_cluster.nf`**
  — a one-line comment at the corrector call site pointing back at the
  script's docstring.

## Testing

- `bin/restore_mmseqs_cluster_ids.py`: unit tests against synthetic fixtures
  covering (a) each of mmseqs's 13 recognized prefix conventions (or the
  general `|`-field-matching rule, per Decision §1.5's chosen
  implementation) → corrected output matches full headers; (b) a
  MAG/prodigal-style FASTA + an already-full-header cluster.tsv → output
  unchanged (idempotent no-op); (c) two distinct headers colliding to the
  same bare form → exits non-zero, naming both headers; (d) a cluster.tsv ID
  matching neither a bare nor full form → exits non-zero.
- `bin/extract_family_seqs.py` / `bin/profile_to_matrix.py`: update
  `tests/test_extract_family_seqs.py` and `tests/test_profile_to_matrix.py`
  to assert full-token (not bare) output/behavior, matching the corrected
  cluster.tsv contract.
- Nextflow-level: no existing `.nf.test`/Nextflow-level fixture for
  `MMSEQS_CLUSTER`/`MMSEQS_FAMILY_CLUSTER` was found in `tests/` — this is a
  net-new test to add, not an extension of an existing one; assert the
  *published* `*_cluster.tsv` is full-header-keyed for a UniProt-style input
  fixture.
- Manual/integration check against `pezizo_set1`: run the backfill command
  (Section 4), then re-run `report_data.py`'s payload build (or regenerate
  the report) and confirm rows are now assigned to families (currently
  0/35,597); re-run `summarize_tblastn.py` against the corrected
  `clusters_cluster.tsv` and confirm the previously-539 mis-reported rows
  now carry correct `tblastn_outgroup_hits` data. **Do not** expect
  `build_alignment_shards.py` to produce shards from this alone —
  `pezizo_set1`'s `tblastn/*.tblastn.tsv` files need regenerating first
  (10-vs-13-column staleness, Problem section) — that's a separate,
  pre-existing fix this spec does not cover.

## Out of scope

- TBLASTN `sseqid` / genome-scaffold-DNA accession handling — confirmed
  clean, no change.
- Regenerating `pezizo_set1`'s stale (10-column) `tblastn/*.tblastn.tsv`
  files — a real, separate pre-existing bug this session surfaced as a side
  effect, needed before alignment shards work for this specific study, but
  unrelated to the ID-space fix itself.
- Migrating the ~16-22 files that read `cluster_tsv` with no normalization
  at all — deliberately not touched (Decision explains why: fixing the
  source makes this unnecessary). The two files that already *attempt*
  normalization (`extract_family_seqs.py`, `profile_to_matrix.py`) **are**
  touched, since their existing partial workaround becomes actively wrong
  once the source is corrected.
