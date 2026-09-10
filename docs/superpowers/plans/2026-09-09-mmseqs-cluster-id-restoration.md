# mmseqs Cluster ID Restoration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix mmseqs's silent header-collapsing in `*_cluster.tsv` (bare accession instead of the full token every other artifact uses as `protein_id`) at the two places mmseqs actually produces it, so family grouping, TBLASTN rep→member hit expansion, and alignment shards all work correctly for UniProt-sourced studies — with zero changes to the ~16-22 unaffected consumer files.

**Architecture:** A new script, `bin/restore_mmseqs_cluster_ids.py`, runs as one more shell command appended to each of the two Nextflow processes that invoke `mmseqs easy-cluster`, rewriting the just-produced `*_cluster.tsv` back to full-header IDs before `publishDir` copies it out. Two existing files that already contain a *partial*, narrower workaround for this same bug (`bin/extract_family_seqs.py`, `bin/profile_to_matrix.py`) have that workaround removed, since it becomes actively wrong once the source file is corrected.

**Tech Stack:** Python 3.12 (pytest), Nextflow DSL2 (Groovy `shell:` blocks).

**Spec:** `docs/superpowers/specs/2026-09-09-mmseqs-cluster-id-restoration-design.md` (rev. 2, post-Fable-review)

## Global Constraints

- The corrector script **fails loudly** (exit 1) on: two distinct full headers colliding to the same computed bare form; a `*_cluster.tsv` value matching neither a known full header nor its computed bare form. Never silently drop or guess.
- Split on whitespace with plain `str.split()` (not `split(' ')`) when taking a FASTA header's first token — mmseqs splits on any whitespace, including tabs.
- `lib/fasta.py::mmseqs_id()` is **repurposed, not duplicated** — the corrector imports and reuses it directly; no second, competing normalization rule anywhere.
- No changes to any of the ~16-22 files that read `cluster_tsv` with no normalization at all (`lib/clusters.py`, `bin/build_alignment_shards.py`, `bin/summarize_tblastn.py`, `bin/make_novelties.py`, etc.) — the fix at the source makes changing them unnecessary.
- `bin/extract_family_seqs.py` and `bin/profile_to_matrix.py` **are** touched — their existing `mmseqs_id()` calls must be removed, since keeping them would double-normalize an already-corrected full-token ID and break the lookup.

---

## Task 1: Extend `lib/fasta.py::mmseqs_id()` to mmseqs's full recognized-prefix set

**Files:**
- Modify: `lib/fasta.py:1-27`
- Test: `tests/test_fasta.py`

**Interfaces:**
- Produces: `mmseqs_id(header_id: str) -> str` — same signature as today, extended behavior. Consumed by Task 2's corrector script (new import) and already imported by Tasks 4/5's files (those imports get *removed* in those tasks, not here).

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_fasta.py` (after the existing `test_mmseqs_id_is_a_no_op_for_non_uniprot_headers`):

```python
def test_mmseqs_id_handles_ncbi_style_prefixes():
    # Empirically confirmed against a real mmseqs2 easy-cluster run (this
    # session's spec review): mmseqs recognizes several more single-pipe
    # NCBI-style prefixes beyond sp|/tr|, all taking the second field.
    assert mmseqs_id('gb|AAA12345.1|') == 'AAA12345.1'
    assert mmseqs_id('ref|XP_000002.1|') == 'XP_000002.1'
    assert mmseqs_id('pdb|1ABC|A') == '1ABC'
    assert mmseqs_id('bbs|123456') == '123456'
    assert mmseqs_id('lcl|NC_000001.1_prot_XP_1_1') == 'NC_000001.1_prot_XP_1_1'
    assert mmseqs_id('pat|US|123456') == 'US'
    assert mmseqs_id('cl|some_locus') == 'some_locus'


def test_mmseqs_id_handles_gnl_and_double_pipe_prefixes():
    # "gnl|db|LOCUS" takes the THIRD field (not the second, unlike the
    # single-pipe prefixes above) -- confirmed empirically.
    assert mmseqs_id('gnl|DB|LOCUS_001') == 'LOCUS_001'
    # "pir||ACCESSION" / "prf||ACCESSION" -- the middle field is empty by
    # convention; the accession is the third field.
    assert mmseqs_id('pir||S12345') == 'S12345'
    assert mmseqs_id('prf||1234567A') == '1234567A'


def test_mmseqs_id_handles_ncbi_compound_gi_header():
    # "gi|N|db|ACCESSION|..." -- NCBI's classic compound defline; mmseqs
    # reports the accession (4th field), not the gi number.
    assert mmseqs_id('gi|12345|ref|XP_000001.1|') == 'XP_000001.1'


def test_mmseqs_id_unrecognized_prefix_is_a_no_op():
    # A pipe-containing header whose first field is NOT one of mmseqs's
    # known prefixes is left completely unchanged -- confirmed empirically
    # (mmseqs does not special-case this).
    assert mmseqs_id('Pchr|PCH_Pc18g05710.1') == 'Pchr|PCH_Pc18g05710.1'
    assert mmseqs_id('UniRef90_A0A000') == 'UniRef90_A0A000'
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && pixi run pytest tests/test_fasta.py -v`
Expected: FAIL — the 4 new tests fail (e.g. `mmseqs_id('gb|AAA12345.1|')` currently returns the input unchanged, since only `sp|`/`tr|` are recognized today); the 2 pre-existing tests (`test_mmseqs_id_extracts_uniprot_accession`, `test_mmseqs_id_is_a_no_op_for_non_uniprot_headers`) still PASS (don't break those).

- [ ] **Step 3: Replace `mmseqs_id()`'s implementation**

Edit `lib/fasta.py`. Replace lines 6-27 (the comment block, `_MMSEQS_UNIPROT_ID_RE`, and `mmseqs_id()`) with:

```python
# mmseqs2 (createdb/easy-cluster) recognizes several FASTA defline conventions
# and reports a collapsed field as the sequence ID in its own *_cluster.tsv,
# rather than the whole first-whitespace-delimited header token every other
# tool in this pipeline uses (Biopython's SeqRecord.id, and the
# "grep '^>' | sed 's/[[:space:]].*//'" protein-map extraction in
# workflows/profile_search.nf's SEED_PROTEIN_MAP). BLAST+ independently does
# the same header recognition, so the pairwise (--cluster_tool pairwise) path
# never notices this divergence in its TBLASTN step -- it only surfaces
# wherever a *_cluster.tsv itself is read.
#
# Confirmed empirically against a real mmseqs2 build (issues #85, and this
# session's mmseqs-cluster-ID-restoration design review) -- not
# documentation-derived, since mmseqs's own docs don't spell this out:
#   sp|ACC|NAME, tr|ACC|NAME, gb|ACC|, ref|ACC|, pdb|ACC|CHAIN, bbs|ACC,
#     lcl|ACC, pat|COUNTRY|ACC, cl|ACC          -> second field
#   gnl|db|ACC                                   -> third field
#   pir||ACC, prf||ACC (empty second field)      -> third field
#   gi|N|db|ACC|...  (NCBI's compound defline)    -> fourth field
# Any OTHER pipe-containing header (an unrecognized first field, e.g. this
# repo's own "Pchr|PCH_..." study-specific locus tags) is left completely
# unchanged -- mmseqs does not special-case it, confirmed empirically.
_MMSEQS_SINGLE_PIPE_PREFIXES = {
    'sp', 'tr', 'gb', 'ref', 'pdb', 'bbs', 'lcl', 'pat', 'cl',
}
_MMSEQS_DOUBLE_PIPE_PREFIXES = {'pir', 'prf'}


def mmseqs_id(header_id: str) -> str:
    """Normalize a FASTA header's first token to what mmseqs2 would report as
    that sequence's id in its own *_cluster.tsv output -- a no-op for any
    header whose first field isn't one of mmseqs's recognized prefixes."""
    fields = header_id.split('|')
    if len(fields) < 2:
        return header_id
    prefix = fields[0].lower()
    if prefix == 'gi' and len(fields) >= 4:
        return fields[3]
    if prefix == 'gnl' and len(fields) >= 3:
        return fields[2]
    if prefix in _MMSEQS_DOUBLE_PIPE_PREFIXES and len(fields) >= 3 and fields[1] == '':
        return fields[2]
    if prefix in _MMSEQS_SINGLE_PIPE_PREFIXES:
        return fields[1]
    return header_id
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && pixi run pytest tests/test_fasta.py -v`
Expected: PASS (all 8 tests: the 2 original + 6 new)

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
git add lib/fasta.py tests/test_fasta.py
git commit -m "$(cat <<'EOF'
Extend mmseqs_id() to mmseqs's full recognized-prefix set

Was sp|/tr| only (issue #85's original, UniProt-specific fix); now
covers all prefix conventions mmseqs2 is confirmed to collapse
(gb/ref/pdb/bbs/lcl/pat/cl, gnl, pir/prf, and NCBI's compound gi|
defline), each empirically verified against a real mmseqs2 build.
Sets up the corrector script (next task) to reuse this single rule
rather than reimplementing it.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Create `bin/restore_mmseqs_cluster_ids.py`

**Files:**
- Create: `bin/restore_mmseqs_cluster_ids.py`
- Test: `tests/test_restore_mmseqs_cluster_ids.py`

**Interfaces:**
- Consumes: `lib/fasta.py::mmseqs_id(header_id: str) -> str` (Task 1).
- Produces: a CLI script, `restore_mmseqs_cluster_ids.py --input-fasta <path> --cluster-tsv <path>`, overwriting `--cluster-tsv` in place. Consumed by Task 3 (Nextflow wiring) and Task 7 (manual backfill).

- [ ] **Step 1: Write the failing tests**

Create `tests/test_restore_mmseqs_cluster_ids.py`:

```python
"""Unit tests for bin/restore_mmseqs_cluster_ids.py.

mmseqs easy-cluster silently collapses certain FASTA header conventions
(sp|ACC|NAME, tr|ACC|NAME, and others -- see lib/fasta.py::mmseqs_id()) to a
shorter field in its own *_cluster.tsv output, while every other pipeline
artifact keeps the full header token as protein_id. This script corrects
*_cluster.tsv back to the full-header form immediately after mmseqs
produces it, using the same FASTA it was run on to reconstruct the mapping.
"""
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / 'bin' / 'restore_mmseqs_cluster_ids.py'


def _run(fasta_path, cluster_tsv_path):
    return subprocess.run(
        [sys.executable, str(BIN),
         '--input-fasta', str(fasta_path),
         '--cluster-tsv', str(cluster_tsv_path)],
        capture_output=True, text=True,
    )


def test_restores_uniprot_style_collapsed_ids(tmp_path):
    fasta = tmp_path / 'seed.faa'
    fasta.write_text(
        ">sp|O74225|YCF1_SCHPO Uncharacterized OS=S. pombe\nMKVLA\n"
        ">tr|A8PCG9|A8PCG9_ASPNG Putative OS=A. niger\nMKVLB\n"
    )
    cluster_tsv = tmp_path / 'cluster.tsv'
    # mmseqs's own raw output shape: bare accessions.
    cluster_tsv.write_text("O74225\tO74225\nO74225\tA8PCG9\n")

    r = _run(fasta, cluster_tsv)
    assert r.returncode == 0, r.stderr
    assert cluster_tsv.read_text() == (
        "sp|O74225|YCF1_SCHPO\tsp|O74225|YCF1_SCHPO\n"
        "sp|O74225|YCF1_SCHPO\ttr|A8PCG9|A8PCG9_ASPNG\n"
    )


def test_idempotent_no_op_for_non_uniprot_headers(tmp_path):
    fasta = tmp_path / 'seed.faa'
    fasta.write_text(">UHM102bin47__k141_1_1 desc\nMKVLA\n>UHM102bin47__k141_2_1 desc\nMKVLB\n")
    cluster_tsv = tmp_path / 'cluster.tsv'
    original = "UHM102bin47__k141_1_1\tUHM102bin47__k141_1_1\nUHM102bin47__k141_1_1\tUHM102bin47__k141_2_1\n"
    cluster_tsv.write_text(original)

    r = _run(fasta, cluster_tsv)
    assert r.returncode == 0, r.stderr
    assert cluster_tsv.read_text() == original  # unchanged -- already full-token


def test_running_twice_is_idempotent(tmp_path):
    fasta = tmp_path / 'seed.faa'
    fasta.write_text(">sp|O74225|YCF1_SCHPO desc\nMKVLA\n")
    cluster_tsv = tmp_path / 'cluster.tsv'
    cluster_tsv.write_text("O74225\tO74225\n")

    r1 = _run(fasta, cluster_tsv)
    assert r1.returncode == 0, r1.stderr
    once = cluster_tsv.read_text()
    r2 = _run(fasta, cluster_tsv)
    assert r2.returncode == 0, r2.stderr
    assert cluster_tsv.read_text() == once


def test_fails_loud_on_ambiguous_bare_id_collision(tmp_path):
    fasta = tmp_path / 'seed.faa'
    # Two distinct full headers that mmseqs would both collapse to "O74225".
    fasta.write_text(
        ">sp|O74225|NAME_A desc one\nMKVLA\n"
        ">tr|O74225|NAME_B desc two\nMKVLB\n"
    )
    cluster_tsv = tmp_path / 'cluster.tsv'
    cluster_tsv.write_text("O74225\tO74225\n")

    r = _run(fasta, cluster_tsv)
    assert r.returncode != 0
    assert 'O74225' in r.stderr


def test_fails_loud_on_unreconcilable_id(tmp_path):
    fasta = tmp_path / 'seed.faa'
    fasta.write_text(">sp|O74225|YCF1_SCHPO desc\nMKVLA\n")
    cluster_tsv = tmp_path / 'cluster.tsv'
    # An id that matches neither the full header nor its computed bare form.
    cluster_tsv.write_text("SOMETHING_ELSE\tSOMETHING_ELSE\n")

    r = _run(fasta, cluster_tsv)
    assert r.returncode != 0
    assert 'SOMETHING_ELSE' in r.stderr
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && pixi run pytest tests/test_restore_mmseqs_cluster_ids.py -v`
Expected: FAIL — `bin/restore_mmseqs_cluster_ids.py` doesn't exist yet (`FileNotFoundError`/non-zero exit from the subprocess call, or `ENOENT`-shaped error).

- [ ] **Step 3: Write the script**

Create `bin/restore_mmseqs_cluster_ids.py`:

```python
#!/usr/bin/env python3
"""Restore mmseqs easy-cluster's collapsed FASTA header IDs in *_cluster.tsv.

mmseqs easy-cluster recognizes several FASTA defline conventions (UniProt's
sp|ACC|NAME/tr|ACC|NAME, and several NCBI-style ones -- see lib/fasta.py's
mmseqs_id() docstring) and reports a collapsed field as the sequence ID in
its own *_cluster.tsv output, while its *_rep_seq.fasta/*_all_seqs.fasta
outputs from the SAME run keep the full header token. Every other artifact
in this pipeline (candidates.txt, presence_matrix.tsv, novelties.<Short>.tsv,
TBLASTN's qseqid, the report payload) treats the full header token as
protein_id -- so any code joining a *_cluster.tsv row against protein_id
silently matches nothing for a UniProt-sourced study.

This script corrects *_cluster.tsv in place, immediately after mmseqs
produces it (wired into modules/mmseqs_cluster.nf and
modules/mmseqs_family_cluster.nf), using the exact FASTA fed to mmseqs to
reconstruct {collapsed_id: full_header} and rewrite every cluster.tsv value
back to the full header -- before any of this pipeline's ~16-22 cluster.tsv
consumers ever reads the file. See
docs/superpowers/specs/2026-09-09-mmseqs-cluster-id-restoration-design.md.

Fails loudly (never silently mismatches) on:
  - two distinct FASTA headers collapsing to the same id (a genuine
    accession collision in the input, most likely when concatenating
    multiple independently-sourced proteomes for family clustering)
  - a cluster.tsv value matching neither a known full header nor its
    computed collapsed form (an mmseqs output this script cannot reconcile)

Usage:
  restore_mmseqs_cluster_ids.py \
      --input-fasta <the FASTA mmseqs clustered> \
      --cluster-tsv <mmseqs easy-cluster's *_cluster.tsv, corrected in place>
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from fasta import mmseqs_id  # noqa: E402


def build_collapsed_to_full(fasta_path: str) -> dict[str, str]:
    """{collapsed_id: full_id} for every record in fasta_path, using the same
    header the FASTA itself carries (first whitespace-delimited token) as
    full_id. Fails loud on a collapsed-id collision between two distinct
    full ids."""
    collapsed_to_full: dict[str, str] = {}
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            full_id = line[1:].split()[0]
            collapsed_id = mmseqs_id(full_id)
            if collapsed_id in collapsed_to_full and collapsed_to_full[collapsed_id] != full_id:
                sys.exit(
                    f"ERROR: {fasta_path}: both '{collapsed_to_full[collapsed_id]}' and "
                    f"'{full_id}' collapse to the same mmseqs id '{collapsed_id}' -- "
                    f"cannot unambiguously restore cluster.tsv"
                )
            collapsed_to_full[collapsed_id] = full_id
    return collapsed_to_full


def restore_cluster_tsv(cluster_tsv_path: str, collapsed_to_full: dict[str, str]) -> str:
    """Return the corrected cluster.tsv content (rep\tmember per line)."""
    full_ids = set(collapsed_to_full.values())
    out_lines = []
    with open(cluster_tsv_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            fixed = []
            for x in parts:
                if x in collapsed_to_full:
                    fixed.append(collapsed_to_full[x])
                elif x in full_ids:
                    fixed.append(x)  # already full -- idempotent no-op
                else:
                    sys.exit(
                        f"ERROR: {cluster_tsv_path}: id '{x}' matches neither a known "
                        f"full header nor a computed mmseqs id -- cannot restore"
                    )
            out_lines.append('\t'.join(fixed))
    return '\n'.join(out_lines) + ('\n' if out_lines else '')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input-fasta', required=True, dest='input_fasta',
                     help='The exact FASTA fed to mmseqs easy-cluster')
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv',
                     help='mmseqs easy-cluster *_cluster.tsv, corrected in place')
    args = ap.parse_args()

    collapsed_to_full = build_collapsed_to_full(args.input_fasta)
    corrected = restore_cluster_tsv(args.cluster_tsv, collapsed_to_full)
    with open(args.cluster_tsv, 'w') as fh:
        fh.write(corrected)

    print(f"Restored {len(collapsed_to_full)} header(s) in {args.cluster_tsv}", file=sys.stderr)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
```

Make it executable: `chmod +x bin/restore_mmseqs_cluster_ids.py`

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && pixi run pytest tests/test_restore_mmseqs_cluster_ids.py -v`
Expected: PASS (all 5 tests)

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
git add bin/restore_mmseqs_cluster_ids.py tests/test_restore_mmseqs_cluster_ids.py
chmod +x bin/restore_mmseqs_cluster_ids.py
git add bin/restore_mmseqs_cluster_ids.py
git commit -m "$(cat <<'EOF'
Add bin/restore_mmseqs_cluster_ids.py

Corrects mmseqs easy-cluster's *_cluster.tsv output back to full-header
IDs, using the exact FASTA mmseqs was run on. Fails loud on an
unreconcilable id or a genuine header collision rather than silently
mismatching. Not yet wired into the pipeline (next task).

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Wire the corrector into both mmseqs clustering processes

**Files:**
- Modify: `modules/mmseqs_cluster.nf`
- Modify: `modules/mmseqs_family_cluster.nf`

**Interfaces:**
- Consumes: `bin/restore_mmseqs_cluster_ids.py` (Task 2), invoked bare (already on `PATH` inside the container, matching every other `bin/*.py` script call in these same process blocks).
- Produces: both processes' *published* `*_cluster.tsv` output is now full-header-keyed. No output-declaration change, so no downstream file changes.

- [ ] **Step 1: Edit `modules/mmseqs_cluster.nf`**

Current `shell:` block ends with the `mmseqs easy-cluster ...` invocation (inside the `if [ -s ... ]` guard's else-path). Add the corrector call immediately after it, using the process's own `candidates_fa` input and `out_prefix` val:

```groovy
    mmseqs easy-cluster \
        !{candidates_fa} \
        !{out_prefix} \
        ${tmpdir}/mmseqs_$$ \
        --threads !{task.cpus} \
        --min-seq-id 0.3 \
        -c 0.8 \
        --cov-mode 0

    # mmseqs collapses certain FASTA header conventions to a shorter field in
    # its own *_cluster.tsv (not in *_rep_seq.fasta/*_all_seqs.fasta) -- see
    # bin/restore_mmseqs_cluster_ids.py's docstring. Restore it here, before
    # this process's output is published and any downstream script reads it.
    restore_mmseqs_cluster_ids.py \
        --input-fasta !{candidates_fa} \
        --cluster-tsv !{out_prefix}_cluster.tsv
    '''
```

(The trailing `'''` above is the existing shell-block terminator already in the file — this snippet shows only the new lines inserted just before it, immediately after the `mmseqs easy-cluster` command and its arguments.)

- [ ] **Step 2: Edit `modules/mmseqs_family_cluster.nf`**

Same pattern, using this process's own `seed_fa` input:

```groovy
    mmseqs easy-cluster \
        !{seed_fa} \
        families \
        ${tmpdir}/mmseqs_fam_$$ \
        --threads !{task.cpus} \
        -s !{params.family_sensitivity} \
        --min-seq-id !{params.family_min_seq_id} \
        -c !{params.family_cov} \
        --cov-mode !{params.family_cov_mode} \
        --cluster-mode !{params.family_cluster_mode}

    # See modules/mmseqs_cluster.nf's identical comment / bin/restore_mmseqs_
    # cluster_ids.py's docstring -- same correction, same reason.
    restore_mmseqs_cluster_ids.py \
        --input-fasta !{seed_fa} \
        --cluster-tsv families_cluster.tsv
    '''
```

- [ ] **Step 3: Verify Groovy/shell syntax is well-formed**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && nextflow config -show-matched-files 2>&1 | head -5 || true` — this is a light sanity check that Nextflow can still parse the project after the edit (a full run needs test data/containers this environment may not have; a syntax error in the DSL would still surface here). If a fuller check is available, additionally run: `nextflow run . --help 2>&1 | tail -20` and confirm no Groovy parse error appears (a `MissingMethodException`/parse stack trace naming `mmseqs_cluster.nf` or `mmseqs_family_cluster.nf` would mean the edit broke the `shell:` block).
Expected: no Groovy/Nextflow parse error referencing either modified file.

- [ ] **Step 4: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
git add modules/mmseqs_cluster.nf modules/mmseqs_family_cluster.nf
git commit -m "$(cat <<'EOF'
Wire restore_mmseqs_cluster_ids.py into both mmseqs clustering processes

MMSEQS_CLUSTER (--cluster_tool pairwise) and MMSEQS_FAMILY_CLUSTER
(--cluster_tool mmseqs and novelty_discovery, both routing through this
one shared module) now correct their own *_cluster.tsv output in place
before publishDir copies it out -- covers all three --cluster_tool
values with these two edits.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Remove `extract_family_seqs.py`'s now-obsolete `mmseqs_id()` workaround

**Files:**
- Modify: `bin/extract_family_seqs.py:31,69-95`
- Test: `tests/test_extract_family_seqs.py`

**Interfaces:**
- Consumes: nothing new — after Task 3, this script's `--cluster-tsv` input is already full-header-keyed, matching its own `--fasta` input directly.

- [ ] **Step 1: Update the (already-passing, pre-fix-contract) regression test to the corrected contract**

In `tests/test_extract_family_seqs.py`, replace the existing
`test_uniprot_style_headers_resolve_via_mmseqs_id` test (the one asserting
bare-form `>O74225`/`>A8PCG9` output) with:

```python
def test_uniprot_style_headers_join_directly_no_normalization_needed(tmp_path):
    """Post mmseqs-cluster-ID-restoration fix (2026-09-09): --cluster-tsv now
    holds the SAME full "sp|ACC|NAME"/"tr|ACC|NAME" token every other
    artifact uses as protein_id (bin/restore_mmseqs_cluster_ids.py corrects
    mmseqs's own *_cluster.tsv output immediately after it's produced -- see
    that script's docstring). This replaces the old mmseqs_id()-based
    regression test (issue #85), which encoded the pre-fix, bare-accession
    cluster.tsv contract this script no longer needs to compensate for."""
    cluster_tsv = tmp_path / 'cluster.tsv'
    cluster_tsv.write_text(
        "sp|O74225|YCF1_SCHPO\tsp|O74225|YCF1_SCHPO\n"
        "sp|O74225|YCF1_SCHPO\ttr|A8PCG9|A8PCG9_ASPNG\n"
    )
    fasta = tmp_path / 'seed.faa'
    fasta.write_text(
        ">sp|O74225|YCF1_SCHPO Uncharacterized OS=S. pombe OX=284812 GN=x PE=1 SV=1\nMKVLA\n"
        ">tr|A8PCG9|A8PCG9_ASPNG Putative OS=A. niger OX=5061 PE=4 SV=1\nMKVLB\n"
    )
    outdir = tmp_path / 'out'
    outdir.mkdir()
    cmd = [
        sys.executable, str(BIN),
        '--cluster-tsv', str(cluster_tsv),
        '--fasta', str(fasta),
        '--min-members', '2',
        '--outdir', str(outdir),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    fam_files = sorted(outdir.glob('fam_*.faa'))
    assert len(fam_files) == 1
    assert fam_files[0].read_text() == (
        '>sp|O74225|YCF1_SCHPO\nMKVLA\n>tr|A8PCG9|A8PCG9_ASPNG\nMKVLB\n'
    )
```

- [ ] **Step 2: Run the updated test to verify it fails against the current (unfixed) script**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && pixi run pytest tests/test_extract_family_seqs.py::test_uniprot_style_headers_join_directly_no_normalization_needed -v`
Expected: FAIL — the script still calls `mmseqs_id()` on `rid` when building `by_mmseqs_id`, so the dict is keyed by bare `O74225`/`A8PCG9`, and lookups by the now-full-token cluster.tsv members (`sp|O74225|YCF1_SCHPO`, `tr|A8PCG9|A8PCG9_ASPNG`) miss, hitting the `sys.exit("ERROR: cluster member ... not found")` path (non-zero exit).

- [ ] **Step 3: Remove the `mmseqs_id()` workaround**

Edit `bin/extract_family_seqs.py`:

Change the import (line 31) from:
```python
from fasta import mmseqs_id, read_fasta  # noqa: E402
```
to:
```python
from fasta import read_fasta  # noqa: E402
```

Change the lookup-dict construction (lines 69-75) from:
```python
    # Look up by mmseqs's own id convention (see lib/fasta.py's mmseqs_id
    # docstring, issue #85): a plain header normalizes to itself here, so
    # this is a no-op for non-UniProt-shaped headers. First-id-wins on a
    # normalization collision, matching read_fasta()'s own dedup policy.
    by_mmseqs_id = {}
    for rid, rec in records.items():
        by_mmseqs_id.setdefault(mmseqs_id(rid), rec)
```
to:
```python
    # cluster.tsv member ids now match this FASTA's own headers directly --
    # bin/restore_mmseqs_cluster_ids.py corrects mmseqs's own *_cluster.tsv
    # output to full-header form immediately after it's produced (see that
    # script's docstring). No normalization needed here any more.
    records_by_id = records
```

And update the two remaining references (originally lines 91, 95) from `by_mmseqs_id` to `records_by_id`:
```python
                if m not in records_by_id:
                    # A member id absent from the FASTA means a mismatch between the
                    # cluster tsv and the sequences it was built from — fail loud.
                    sys.exit(f"ERROR: cluster member '{m}' not found in {args.fasta}")
                out.write(f'>{m}\n{records_by_id[m].seq}\n')
```

- [ ] **Step 4: Run the full test file to verify everything passes**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && pixi run pytest tests/test_extract_family_seqs.py -v`
Expected: PASS (all tests, including the rewritten one and the pre-existing `--max-members`/oversized-family tests, which use bare non-pipe IDs and are unaffected by this change since `mmseqs_id()` was always a no-op for those)

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
git add bin/extract_family_seqs.py tests/test_extract_family_seqs.py
git commit -m "$(cat <<'EOF'
Remove extract_family_seqs.py's now-obsolete mmseqs_id() workaround

cluster.tsv is now full-header-keyed at the source (previous two
tasks), matching this script's --fasta input directly -- keeping the
old mmseqs_id() normalization here would double-normalize an
already-correct id and break every lookup.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Remove `profile_to_matrix.py`'s now-obsolete `mmseqs_id()` workaround

**Files:**
- Modify: `bin/profile_to_matrix.py:43,46-67`
- Test: `tests/test_profile_to_matrix.py`

**Interfaces:**
- Consumes: nothing new — after Task 3, this script's `--cluster-tsv` input is already full-header-keyed, matching `--protein-map`'s own full-token entries directly.

- [ ] **Step 1: Update the (already-passing, pre-fix-contract) regression test to the corrected contract**

In `tests/test_profile_to_matrix.py`, replace the existing
`test_uniprot_style_protein_map_resolves_via_mmseqs_id` test with:

```python
def test_uniprot_style_ids_join_directly_no_normalization_needed(tmp_path):
    """Post mmseqs-cluster-ID-restoration fix (2026-09-09): --cluster-tsv now
    holds the full "sp|ACC|NAME"/"tr|ACC|NAME" token (bin/restore_mmseqs_
    cluster_ids.py corrects mmseqs's own *_cluster.tsv output immediately
    after it's produced), matching --protein-map's already-full tokens
    directly. Replaces the old mmseqs_id()-based regression test (issue
    #85), which encoded the pre-fix, bare-accession cluster.tsv contract."""
    (tmp_path / 'config.csv').write_text(CONFIG)
    (tmp_path / 'cluster.tsv').write_text(
        "sp|pA1|NAME_ONE\tsp|pA1|NAME_ONE\n"
        "sp|pA1|NAME_ONE\ttr|pA2|NAME_TWO\n"
    )
    (tmp_path / 'families.tsv').write_text(
        "family_index\trepresentative_id\tn_members\n"
        "fam_000001\tsp|pA1|NAME_ONE\t2\n"
    )
    (tmp_path / 'protein_map.tsv').write_text(
        "sp|pA1|NAME_ONE\tIn1\n"
        "tr|pA2|NAME_TWO\tIn2\n"
    )
    _write_domtblout(tmp_path / 'In1.family.domtblout', [('t1', 'sp|pA1|NAME_ONE')])
    _write_domtblout(tmp_path / 'In2.family.domtblout', [('t2', 'sp|pA1|NAME_ONE')])
    _write_domtblout(tmp_path / 'Out1.family.domtblout', [])
    matrix, cands = _run(tmp_path)
    assert set(matrix['protein_id']) == {'sp|pA1|NAME_ONE', 'tr|pA2|NAME_TWO'}
    assert cands == ['In1::sp|pA1|NAME_ONE', 'In2::tr|pA2|NAME_TWO']
```

- [ ] **Step 2: Run the updated test to verify it fails against the current (unfixed) script**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && pixi run pytest tests/test_profile_to_matrix.py::test_uniprot_style_ids_join_directly_no_normalization_needed -v`
Expected: FAIL — `load_protein_map()` still keys its map by `mmseqs_id(pid)` (bare `pA1`, `pA2`), so `load_cluster_members()`'s now-full-token members (`sp|pA1|NAME_ONE`, `tr|pA2|NAME_TWO`) never match `protein_to_proteome.get(member)`, hitting the "skip defensively" branch — the matrix comes out empty, failing the `set(matrix['protein_id'])` assertion.

- [ ] **Step 3: Remove the `mmseqs_id()` workaround**

Edit `bin/profile_to_matrix.py`:

Change the import (line 43) from:
```python
from fasta import mmseqs_id  # noqa: E402
```
Delete this line entirely (no longer needed).

Change `load_protein_map()` (lines 46-67) from:
```python
def load_protein_map(path):
    """protein_id -> proteome Short, from a two-column TSV (no header).

    Keys are normalized with mmseqs_id() (issue #85): --protein-map comes from
    workflows/profile_search.nf's SEED_PROTEIN_MAP, which extracts a FASTA
    header's first whitespace-delimited token verbatim (no UniProt-header
    awareness) -- but --cluster-tsv's member ids come from mmseqs2 itself,
    which *does* recognize "sp|ACC|NAME"/"tr|ACC|NAME" headers and reports
    the bare accession. Without this normalization, every member lookup
    against this map silently misses for UniProt-sourced proteomes (not a
    crash -- main() below treats a missing map entry as "skip defensively",
    so the presence matrix would just come out empty).
    """
    mapping = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            pid, short = line.split('\t')[:2]
            mapping[mmseqs_id(pid)] = short
    return mapping
```
to:
```python
def load_protein_map(path):
    """protein_id -> proteome Short, from a two-column TSV (no header).

    --protein-map (workflows/profile_search.nf's SEED_PROTEIN_MAP) and
    --cluster-tsv's member ids are now the same full FASTA-header token --
    bin/restore_mmseqs_cluster_ids.py corrects mmseqs's own *_cluster.tsv
    output to match immediately after it's produced (see that script's
    docstring). No normalization needed here any more.
    """
    mapping = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            pid, short = line.split('\t')[:2]
            mapping[pid] = short
    return mapping
```

- [ ] **Step 4: Run the full test file to verify everything passes**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && pixi run pytest tests/test_profile_to_matrix.py -v`
Expected: PASS (all tests — the other tests in this file use `_setup()`'s bare-ID fixtures throughout for both `cluster.tsv` and `protein_map.tsv`, matching each other consistently regardless of this change, so they're unaffected)

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
git add bin/profile_to_matrix.py tests/test_profile_to_matrix.py
git commit -m "$(cat <<'EOF'
Remove profile_to_matrix.py's now-obsolete mmseqs_id() workaround

cluster.tsv is now full-header-keyed at the source (Task 3), matching
--protein-map's own full-token entries directly -- keeping the old
mmseqs_id() normalization here would double-normalize an
already-correct id and silently empty the presence matrix.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Documentation

**Files:**
- Modify: `METHOD_DESCRIPTION.md` (two insertion points, per pathway)
- Modify: `README.md` (one insertion point)

**Interfaces:**
- None (documentation only).

- [ ] **Step 1: Add the explanation to `METHOD_DESCRIPTION.md`'s `pairwise` section**

Find the line (around line 107) reading:
```
Candidates are extracted, clustered with `mmseqs easy-cluster`
```
Add, immediately after the paragraph containing that sentence finishes:

```markdown

mmseqs recognizes several FASTA defline conventions (UniProt's
`sp|ACC|NAME`/`tr|ACC|NAME`, and a handful of NCBI-style ones) and reports a
collapsed field as the sequence ID in its own `*_cluster.tsv` output
specifically (not in `*_rep_seq.fasta`), while every other artifact in this
pipeline keeps the full header token as `protein_id`. Uncorrected, this
silently defeats gene-family grouping, TBLASTN rep→member hit expansion,
and the report's alignment-shard feature for any UniProt-sourced study.
`bin/restore_mmseqs_cluster_ids.py` runs immediately after every `mmseqs
easy-cluster` invocation (both here and in the `mmseqs`/`novelty_discovery`
pathway below) to restore the full header before anything else reads the
file — see that script's docstring, and
`docs/superpowers/specs/2026-09-09-mmseqs-cluster-id-restoration-design.md`
for the full design record.
```

- [ ] **Step 2: Add a shorter pointer to `METHOD_DESCRIPTION.md`'s `mmseqs` section**

Find the line (around line 141) reading:
```
1. `mmseqs cluster` (cascaded, `-s 7`, `--min-seq-id 0.3 -c 0.8 --cov-mode 0
   --cluster-mode 0`) clusters the **whole seed group** (ingroup for novelty,
```
Add, immediately after that numbered step's paragraph:

```markdown

   (Same `*_cluster.tsv` ID-restoration step as the `pairwise` pathway above
   applies here too — see that section.)
```

- [ ] **Step 3: Add the pointer sentence to `README.md`**

Find the "## Two analysis approaches: pairwise vs family-profile clustering" section (around line 46). Add one sentence at the end of that section's introductory paragraph (before the `### --cluster_tool pairwise` subsection begins):

```markdown

(Both pathways cluster with mmseqs2, which has a header-parsing quirk that
needs correcting after every clustering step — see `METHOD_DESCRIPTION.md`
for the full explanation.)
```

- [ ] **Step 4: Verify the docs still render as valid Markdown**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && python3 -c "import pathlib; [print(p, 'OK') for p in [pathlib.Path('METHOD_DESCRIPTION.md'), pathlib.Path('README.md')] if p.exists()]"`
Expected: both files print `OK` (this is a minimal existence/no-corruption sanity check — no Markdown linter is configured in this repo per earlier exploration, so this is what's available)

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
git add METHOD_DESCRIPTION.md README.md
git commit -m "$(cat <<'EOF'
Document the mmseqs cluster-ID restoration step

Explains mmseqs's header-collapsing behavior and why
bin/restore_mmseqs_cluster_ids.py exists, in the two places an
algorithm-level reader would look (METHOD_DESCRIPTION.md's clustering
steps for both pathways) plus a pointer from README.md.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Manual backfill/verification against `pezizo_set1`

This is a manual verification pass against the already-computed `pezizo_set1`
study (in the sibling NovInvenio_Investigations repo,
`/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations`) — no fresh
Nextflow run needed for this part, per the spec's Section 4. No code changes
in this task.

**Files touched (regenerated data, not source):**
- `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/results/pezizo_set1/clusters/clusters_cluster.tsv` and `loss_clusters_cluster.tsv` (corrected in place)
- `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/results/pezizo_set1/tblastn_summary.tsv` (regenerated from the corrected cluster.tsv)

- [ ] **Step 1: Run the backfill against the novelty-direction cluster.tsv**

Run:
```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
/bigdata/stajichlab/jstajich/projects/NovInvenio/bin/restore_mmseqs_cluster_ids.py \
    --input-fasta results/pezizo_set1/candidates.fa \
    --cluster-tsv results/pezizo_set1/clusters/clusters_cluster.tsv
```
Expected: prints `Restored N header(s) in results/pezizo_set1/clusters/clusters_cluster.tsv` to stderr, exit 0.

- [ ] **Step 2: Confirm the file is now full-header-keyed**

Run: `grep -c '|' results/pezizo_set1/clusters/clusters_cluster.tsv` (from the `NovInvenio_Investigations` repo root)
Expected: a count equal to the file's total line count (every row now has pipe-containing, full-header IDs) — compare against `wc -l results/pezizo_set1/clusters/clusters_cluster.tsv` to confirm they match.

- [ ] **Step 3: Run the backfill against the loss-direction cluster.tsv, if present**

Run: `ls results/pezizo_set1/clusters/loss_clusters_cluster.tsv results/pezizo_set1/loss_candidates.fa 2>&1`
If both exist, run:
```bash
/bigdata/stajichlab/jstajich/projects/NovInvenio/bin/restore_mmseqs_cluster_ids.py \
    --input-fasta results/pezizo_set1/loss_candidates.fa \
    --cluster-tsv results/pezizo_set1/clusters/loss_clusters_cluster.tsv
```
Expected: same success pattern as Step 1.

- [ ] **Step 4: Re-run `summarize_tblastn.py` against the corrected cluster.tsv and confirm the rep→member expansion now works**

Run:
```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
pixi run python /bigdata/stajichlab/jstajich/projects/NovInvenio/bin/summarize_tblastn.py \
    --hits results/pezizo_set1/tblastn/*.tblastn.tsv \
    --cluster_tsv results/pezizo_set1/clusters/clusters_cluster.tsv \
    --evalue 1e-5 \
    --output /tmp/tblastn_summary_corrected.tsv
```
**Note:** if you already deleted/regenerated `results/pezizo_set1/tblastn/*.tblastn.tsv` per the separate 10-vs-13-column staleness fix (out of scope for this plan, see the spec's "Out of scope" section), use whatever current TBLASTN tsvs exist. If they're still the stale 10-column files, this step will still show whether rep→member expansion itself now works correctly (compare counts against the old `results/pezizo_set1/tblastn_summary.tsv`), even though the underlying hit data is separately known-stale.

Then check:
```bash
awk -F'\t' 'NR>1{h=0; for(i=2;i<=NF;i++) if($i=="1") h=1; if(h) print $1}' /tmp/tblastn_summary_corrected.tsv | grep -c '|'
awk -F'\t' 'NR>1{h=0; for(i=2;i<=NF;i++) if($i=="1") h=1; if(h) print $1}' /tmp/tblastn_summary_corrected.tsv | wc -l
```
Expected: both counts match (every "has a hit" row is now full-header-keyed, same as before the fix — this part already worked, per the Problem section's traced analysis) AND the total number of "has a hit" rows should now be **higher** than the pre-fix count (1,053 non-member-expanded hits before → should include the previously-lost 539 non-rep-member expansions now, once cluster.tsv correctly maps every member to its rep). Compare: `wc -l < /tmp/tblastn_summary_corrected.tsv` should show more total rows with hits than the old `results/pezizo_set1/tblastn_summary.tsv`'s 1,053.

- [ ] **Step 5: Report back — no commit in this task**

This task only regenerates already-gitignored `results/` output in the
sibling NII repo — nothing here needs a plan-step commit (matches the
pattern already established for this repo's xrefs-linkout plan's own Task
6). If Steps 1-4 pass, tell the user what you confirmed (cluster.tsv is
correctly restored; TBLASTN rep→member expansion now recovers the
previously-lost 539 rows) and that they can regenerate `pezizo_set1`'s
`novelties.<Short>.tsv`/reports whenever convenient to see this reflected
end to end — **do not** regenerate or publish those reports yourself without
their explicit go-ahead, and **do not** attempt to fix the separate
10-vs-13-column TBLASTN staleness issue as part of this task (out of scope,
a distinct pre-existing bug this session surfaced).
