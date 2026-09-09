# UniProt Cross-Reference Linkout Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract FungiDB/NCBI Gene/RefSeq/KEGG/Ensembl cross-references already present in cached UniProt `.dat.gz` records, carry them through the annotation-merge pipeline, and render them generically in nf_NovInvenio's static HTML reports — replacing hardcoded per-DB link blocks with a registry-driven loop — then do a manual rebuild of the already-computed `pezizo_set1` study to confirm it end to end.

**Architecture:** Two repos. NII (`NovInvenio_Investigations`) extracts `DR` cross-reference lines from its cached UniProt flat files into a new `xrefs` column, already present in its per-proteome annotation TSVs, and passes it through its existing `merge_uniprot_annotations.py` step. nf_NovInvenio (`NovInvenio`) carries that one new string field through its three report payload builders (`report_data.py`) and renders it via a small DB→URL-template registry in `report_common.py`'s shared `externalLinksNode`, wired into all three HTML templates.

**Tech Stack:** Python 3.12 (pytest), vanilla JS (no framework — hand-rolled DOM builders in `lib/report_common.py`).

**Spec:** `docs/superpowers/specs/2026-09-09-uniprot-xref-linkout-design.md` (rev. 2, post-review)

## Global Constraints

- `xrefs` packing: `DB:id` pairs joined by `|`; unpacking splits each entry on the **first colon only** (`id` may itself contain colons — `KEGG:ncr:NCU10683`, `VEuPathDB:FungiDB:NCU10683`). Never use a fixed-arity `.split(":")` on this field (that convention belongs to `go_ids` only).
- Allow-listed DR types (only these six are extracted): `VEuPathDB`, `GeneID`, `RefSeq` (first line only, protein-accession field only), `KEGG`, `EnsemblFungi`, `EnsemblBacteria` (both: last non-`-` field, since field position/count varies by organism).
- `PDB` is dropped from v1 (unbounded per-protein multiplicity — see spec).
- No new interning pool — `xrefs` is stored inline per row, following the existing `ec`/`af` precedent (1:1 with the protein, nothing to compress).
- Every interpolated id in a JS URL template goes through `encodeURIComponent`.
- `merge_uniprot_annotations.py`'s new column read uses `a.get("xrefs", "")`, never `a["xrefs"]` (existing annotation TSVs on disk predate this column).

---

## Task 1: NII — extract `xrefs` from UniProt `DR` lines

**Files:**
- Modify: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/bin/extract_dat_annotations.py`
- Test: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/tests/test_extract_dat_annotations.py` (new file — this repo has no `tests/` dir yet, even though `pixi.toml`'s `test` task already points at one)

**Interfaces:**
- Produces: `parse_dat_gz(path: Path)` now yields dicts with one more key, `"xrefs"` (a `|`-joined string, `""` when no allow-listed DR line was found for that record). Output TSV gains a matching `xrefs` column, written last.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_extract_dat_annotations.py`:

```python
import gzip
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "bin"))
from extract_dat_annotations import parse_dat_gz  # noqa: E402

FIXTURE_DAT = """\
ID   A7UWL5_NEUCR            Unreviewed;       301 AA.
AC   A7UWL5;
GN   ORFNames=NCU10683;
OX   NCBI_TaxID=367110;
DE   SubName: Full=NRS/ER;
DR   VEuPathDB; FungiDB:NCU10683; -.
DR   GeneID; 5847462; -.
DR   RefSeq; XP_001728253.1; XM_001728201.2.
DR   RefSeq; XP_999999999.1; XM_999999999.1.
DR   KEGG; ncr:NCU10683; -.
DR   EnsemblFungi; EAA34304; EAA34304; NCU06768.
DR   PDB; 6YWS; EM; 2.74 A; K=1-249.
DR   STRING; 367110.A7UWL5; -.
SQ   SEQUENCE   6 AA;  1 MW;  0000000000000000 CRC64;
     MATTEI
//
ID   B0000000_TEST           Unreviewed;       10 AA.
AC   B0000000;
GN   ORFNames=NCU99999;
OX   NCBI_TaxID=367110;
DE   SubName: Full=Test with Scer-style Ensembl field order;
DR   EnsemblFungi; YML051W_mRNA; YML051W; YML051W.
SQ   SEQUENCE   6 AA;  1 MW;  0000000000000000 CRC64;
     MATTEI
//
"""


def _write_fixture(tmp_path) -> Path:
    p = tmp_path / "fixture.dat.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(FIXTURE_DAT)
    return p


def test_xrefs_packs_allowlisted_db_types_first_colon_only(tmp_path):
    path = _write_fixture(tmp_path)
    records = {r["accession"]: r for r in parse_dat_gz(path)}
    xrefs = records["A7UWL5"]["xrefs"].split("|")
    assert "VEuPathDB:FungiDB:NCU10683" in xrefs
    assert "GeneID:5847462" in xrefs
    assert "KEGG:ncr:NCU10683" in xrefs


def test_xrefs_keeps_only_first_refseq_line(tmp_path):
    path = _write_fixture(tmp_path)
    records = {r["accession"]: r for r in parse_dat_gz(path)}
    xrefs = records["A7UWL5"]["xrefs"].split("|")
    refseq_entries = [x for x in xrefs if x.startswith("RefSeq:")]
    assert refseq_entries == ["RefSeq:XP_001728253.1"]


def test_xrefs_ensemblfungi_takes_last_field_regardless_of_layout(tmp_path):
    path = _write_fixture(tmp_path)
    records = {r["accession"]: r for r in parse_dat_gz(path)}
    a7 = records["A7UWL5"]["xrefs"].split("|")
    b0 = records["B0000000"]["xrefs"].split("|")
    assert "EnsemblFungi:NCU06768" in a7   # Ncra layout: protein;protein;gene
    assert "EnsemblFungi:YML051W" in b0    # Scer layout: transcript;protein;gene


def test_xrefs_drops_pdb_and_non_allowlisted_dbs(tmp_path):
    path = _write_fixture(tmp_path)
    records = {r["accession"]: r for r in parse_dat_gz(path)}
    xrefs = records["A7UWL5"]["xrefs"]
    assert "PDB" not in xrefs
    assert "STRING" not in xrefs
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations && pixi run pytest tests/test_extract_dat_annotations.py -v`
Expected: FAIL — `parse_dat_gz` records have no `"xrefs"` key (or `KeyError`).

- [ ] **Step 3: Implement `xrefs` extraction**

Edit `bin/extract_dat_annotations.py`. Add near the other module-level regexes/constants (after `DE_EC_RE`):

```python
XREF_FIRST_FIELD_DBS = {"VEuPathDB", "GeneID", "KEGG", "RefSeq"}
XREF_LAST_FIELD_DBS = {"EnsemblFungi", "EnsemblBacteria"}
# DR line dispatch, e.g. "DR   GeneID; 5847462; -." -> db="GeneID", rest="5847462; -."
# Kept as one generic split (not a per-DB regex) so adding a 7th allow-listed DB
# later is a one-line set addition, and so a DR line for a DB we don't care
# about (the overwhelming majority -- EMBL/PANTHER/STRING/etc.) costs one
# split+lookup instead of several failed regex .match() attempts (this parser's
# own docstring flags per-line performance as the reason it isn't Biopython;
# Ncra alone has ~250k DR lines).
DR_LINE_RE = re.compile(r"^DR\s+(\w+);\s*(.*?)\.?\s*$")


def _first_xref_field(fields: list[str]) -> str:
    return fields[0].strip() if fields else ""


def _last_xref_field(fields: list[str]) -> str:
    # Field position/count for Ensembl* DR lines is not stable across organisms
    # (Ncra: protein;protein;gene -- Scer: transcript;protein;gene), but the
    # last non-placeholder field is consistently the stable gene ID, which is
    # also what Ensembl's /id/ resolver handles best.
    cleaned = [f.strip() for f in fields if f.strip() and f.strip() != "-"]
    return cleaned[-1] if cleaned else ""
```

Then extend `parse_dat_gz`'s local state (in the outer function, alongside the
existing `alphafold_id = None` etc.):

```python
    xrefs = []
    xref_seen_dbs = set()  # RefSeq: keep only the first line per record
```

and in `reset()`:

```python
        nonlocal xrefs, xref_seen_dbs
        xrefs = []
        xref_seen_dbs = set()
```

and in the yielded dict (alongside `"alphafold_id": alphafold_id or ""`):

```python
                        "xrefs": "|".join(xrefs),
```

Then, inside the `if line.startswith("DR"):` block, **after** the existing
`DR_ALPHAFOLD_RE` check (i.e. it still falls through to here only when none of
GO/Pfam/InterPro/AlphaFoldDB matched — those `continue` on match already):

```python
                m = DR_LINE_RE.match(line)
                if m:
                    db, rest = m.group(1), m.group(2)
                    if db in XREF_FIRST_FIELD_DBS:
                        if db == "RefSeq" and "RefSeq" in xref_seen_dbs:
                            continue
                        val = _first_xref_field(rest.split(";"))
                        if val:
                            xrefs.append(f"{db}:{val}")
                            xref_seen_dbs.add(db)
                        continue
                    if db in XREF_LAST_FIELD_DBS:
                        val = _last_xref_field(rest.split(";"))
                        if val:
                            xrefs.append(f"{db}:{val}")
                        continue
```

Finally, add `"xrefs"` to the `csv.DictWriter` fieldnames list in `main()`:

```python
        w = csv.DictWriter(fh, fieldnames=["accession", "taxon_id", "gene_name", "description", "go_ids", "pfam_ids", "pfam_names", "interpro_ids", "ec_numbers", "alphafold_id", "xrefs"], delimiter="\t")
```

Update the module docstring's "Output columns:" line to mention `xrefs` and
briefly note the allow-list/packing rule (one sentence, pointing at
`docs/superpowers/specs/2026-09-09-uniprot-xref-linkout-design.md` rather than
re-explaining it).

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations && pixi run pytest tests/test_extract_dat_annotations.py -v`
Expected: PASS (all 4 tests)

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
git add bin/extract_dat_annotations.py tests/test_extract_dat_annotations.py
git commit -m "$(cat <<'EOF'
Extract FungiDB/NCBI Gene/RefSeq/KEGG/Ensembl cross-refs into xrefs column

Adds a DR-line dispatch table to the existing single-pass UniProt .dat.gz
parser, packing allow-listed cross-references as DB:id pairs (first colon
only -- KEGG/VEuPathDB ids carry their own colons). RefSeq keeps only the
first line per record (protein accession, not the paired transcript);
EnsemblFungi/EnsemblBacteria take the last field since layout varies by
organism. PDB is deliberately not extracted (unbounded per-protein count).

See docs/superpowers/specs/2026-09-09-uniprot-xref-linkout-design.md.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: NII — pass `xrefs` through `merge_uniprot_annotations.py`

**Files:**
- Modify: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/bin/merge_uniprot_annotations.py`
- Test: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/tests/test_merge_uniprot_annotations.py` (new)

**Interfaces:**
- Consumes: `extract_dat_annotations.py`'s output TSV rows, each optionally carrying an `"xrefs"` field (Task 1).
- Produces: the merged output TSV gains a `uniprot_xrefs` column.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_merge_uniprot_annotations.py`:

```python
import csv
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / "bin" / "merge_uniprot_annotations.py"


def _write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def test_merge_passes_through_xrefs_column(tmp_path):
    annot = tmp_path / "annot.tsv"
    _write_tsv(
        annot,
        [{"accession": "A7UWL5", "gene_name": "", "description": "", "pfam_ids": "",
          "interpro_ids": "", "go_ids": "", "xrefs": "GeneID:5847462|KEGG:ncr:NCU10683"}],
        ["accession", "gene_name", "description", "pfam_ids", "interpro_ids", "go_ids", "xrefs"],
    )
    matrix = tmp_path / "matrix.tsv"
    _write_tsv(matrix, [{"protein_id": "A7UWL5"}], ["protein_id"])
    out = tmp_path / "out.tsv"

    subprocess.run(
        [sys.executable, str(BIN), "--matrix", str(matrix), "--annotations", str(annot), "--output", str(out)],
        check=True,
    )

    with open(out, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["uniprot_xrefs"] == "GeneID:5847462|KEGG:ncr:NCU10683"


def test_merge_lenient_missing_xrefs_column_is_empty_string(tmp_path):
    # An older annotations TSV (predates Task 1) has no xrefs column at all --
    # merge_uniprot_annotations.py must not KeyError on it.
    annot = tmp_path / "annot.tsv"
    _write_tsv(
        annot,
        [{"accession": "A7UWL5", "gene_name": "", "description": "", "pfam_ids": "",
          "interpro_ids": "", "go_ids": ""}],
        ["accession", "gene_name", "description", "pfam_ids", "interpro_ids", "go_ids"],
    )
    matrix = tmp_path / "matrix.tsv"
    _write_tsv(matrix, [{"protein_id": "A7UWL5"}], ["protein_id"])
    out = tmp_path / "out.tsv"

    subprocess.run(
        [sys.executable, str(BIN), "--matrix", str(matrix), "--annotations", str(annot), "--output", str(out)],
        check=True,
    )

    with open(out, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["uniprot_xrefs"] == ""
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations && pixi run pytest tests/test_merge_uniprot_annotations.py -v`
Expected: FAIL — `KeyError: 'uniprot_xrefs'` on the assertion.

- [ ] **Step 3: Implement the passthrough**

Edit `bin/merge_uniprot_annotations.py`. In `main()`, add `"uniprot_xrefs"` to
`new_cols`:

```python
    new_cols = [
        "uniprot_gene_name", "uniprot_description", "uniprot_pfam_ids", "uniprot_pfam_names",
        "uniprot_interpro_ids", "uniprot_go_ids", "uniprot_ec_numbers", "uniprot_alphafold_id",
        "uniprot_xrefs",
    ]
```

and in the per-row loop, alongside the other `.get(...)`-style assignments
(**use `.get`, not bracket access** — per the Global Constraints, existing
annotation TSVs on disk predate this column):

```python
        row["uniprot_xrefs"] = a.get("xrefs", "")
```

Update the module docstring's "Adds new columns ..." line to include
`uniprot_xrefs`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations && pixi run pytest tests/test_merge_uniprot_annotations.py -v`
Expected: PASS (both tests)

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
git add bin/merge_uniprot_annotations.py tests/test_merge_uniprot_annotations.py
git commit -m "$(cat <<'EOF'
Pass xrefs through merge_uniprot_annotations.py as uniprot_xrefs

Uses .get("xrefs", "") rather than bracket access so annotation TSVs
still on disk (predating extract_dat_annotations.py's new xrefs column)
merge without a hard failure.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: nf_NovInvenio — carry `xrefs` through the three payload builders

**Files:**
- Modify: `/bigdata/stajichlab/jstajich/projects/NovInvenio/lib/report_data.py:24-84` (`ROW_FIELDS` + `build_payload`'s `out_rows.append`, `lib/report_data.py:513-538`), `:568-588` + `:667-686` (`CORE_ROW_FIELDS` + `build_core_payload`), `:706-740` + `:875-899` (`LOSSES_ROW_FIELDS` + `build_losses_payload`)
- Test: `/bigdata/stajichlab/jstajich/projects/NovInvenio/tests/test_report_data.py`

**Interfaces:**
- Consumes: `uniprot_xrefs` column on matrix rows (Task 2), read the same way `uniprot_alphafold_id` already is: `row.get('uniprot_xrefs', '') or ''`.
- Produces: all three `payload['fields']` lists gain `'xrefs'` as their last entry; all three `out_rows.append([...])` calls append the matching value last. **No new `_StringTable`** — inline string, same as `'ec'`/`'af'`.

- [ ] **Step 1: Write the failing tests**

In `tests/test_report_data.py`, add a fixture matrix row carrying
`uniprot_xrefs` and three tests (one per payload builder). Insert near the
existing fixtures (after `SUPPORT_MATRIX`, before the `run_dir` fixture):

```python
# Adds a uniprot_xrefs column to exercise the xrefs passthrough (n1 carries
# one, n2 has none -- the empty-string default path).
XREF_MATRIX = """\
protein_id\tsource_proteome\tNcra\tAfum\tSpom\tScer\tgene_name\tproduct_description\tfunction_source\tuniprot_xrefs
n1\tNcra\t1\t1\t0\t0\tada-1\tall development altered-1\tModelOrg_Ncra\tGeneID:5847462|KEGG:ncr:NCU10683
n2\tAfum\t1\t1\t0\t0\t\t\t\t
"""
```

Add the tests (near the other `build_payload`/`build_core_payload`/
`build_losses_payload` tests):

```python
def test_payload_carries_xrefs_field(tmp_path, samples):
    (tmp_path / 'matrix.tsv').write_text(XREF_MATRIX)
    payload = build_payload(tmp_path / 'matrix.tsv', samples, candidates_fa=None,
                             tblastn_path=None)
    idx = payload['fields'].index('xrefs')
    rows = {r[payload['fields'].index('id')]: r for r in payload['rows']}
    assert rows['n1'][idx] == 'GeneID:5847462|KEGG:ncr:NCU10683'
    assert rows['n2'][idx] == ''


def test_core_payload_carries_xrefs_field(tmp_path, samples):
    core_matrix = XREF_MATRIX.replace('n1\tNcra\t1\t1\t0\t0', 'n1\tNcra\t1\t1\t1\t1') \
                             .replace('n2\tAfum\t1\t1\t0\t0', 'n2\tAfum\t1\t1\t1\t1')
    (tmp_path / 'core.tsv').write_text(core_matrix)
    payload = build_core_payload(tmp_path / 'core.tsv', samples, core_min_frac=0.95)
    idx = payload['fields'].index('xrefs')
    rows = {r[payload['fields'].index('id')]: r for r in payload['rows']}
    assert rows['n1'][idx] == 'GeneID:5847462|KEGG:ncr:NCU10683'


def test_losses_payload_carries_xrefs_field(tmp_path, samples):
    losses_matrix = ("protein_id\tsource_proteome\tNcra\tAfum\tSpom\tScer\t"
                      "gene_name\tproduct_description\tfunction_source\tuniprot_xrefs\n"
                      "loss1\tSpom\t0\t0\t1\t1\t\t\t\tGeneID:999\n")
    (tmp_path / 'losses.tsv').write_text(losses_matrix)
    payload = build_losses_payload(tmp_path / 'losses.tsv', samples)
    idx = payload['fields'].index('xrefs')
    rows = {r[payload['fields'].index('id')]: r for r in payload['rows']}
    assert rows['loss1'][idx] == 'GeneID:999'
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && pixi run pytest tests/test_report_data.py -k xrefs -v`
Expected: FAIL — `'xrefs' is not in list` (ValueError from `payload['fields'].index('xrefs')`).

- [ ] **Step 3: Add the `xrefs` field to all three payload builders**

In `lib/report_data.py`, add `'xrefs'` as the last entry of `ROW_FIELDS`
(after `'af'`'s comment block, i.e. right before `'nov'`):

```python
    'xrefs',     # uniprot_xrefs -- NII bin/extract_dat_annotations.py's DR-line
                 # cross-reference extraction (FungiDB/NCBI Gene/RefSeq/KEGG/
                 # Ensembl), packed "DB:id|DB:id" -- see lib/report_common.py's
                 # XREF_LINK_TEMPLATES for the render side. Not interned: 1:1
                 # with the protein like 'ec'/'af', nothing to compress.
```

And append the matching value as the last element of `build_payload`'s
`out_rows.append([...])` list (after `row.get('uniprot_alphafold_id', '') or ''`,
before `is_nov,`):

```python
            row.get('uniprot_xrefs', '') or '',
```

Wait — `'xrefs'` must be added as the *last* field of `ROW_FIELDS` (after
`'start'`) and the matching append must be the *last* element of the list, to
avoid renumbering every field after it in three independent files at once.
Add `'xrefs'` right after `'start'` in `ROW_FIELDS` (`lib/report_data.py:84`,
after the `'start'` comment block) and append
`row.get('uniprot_xrefs', '') or ''` as the last element of `build_payload`'s
`out_rows.append([...])` (`lib/report_data.py:538`, after `start,`).

Do the same for `CORE_ROW_FIELDS` (add `'xrefs'` after `'start'` at
`lib/report_data.py:588`) and `build_core_payload`'s `out_rows.append([...])`
(append after `start,` at `lib/report_data.py:686`).

Do the same for `LOSSES_ROW_FIELDS` (add `'xrefs'` after `'start'` at
`lib/report_data.py:740`) and `build_losses_payload`'s
`out_rows.append([...])` (append after `start,` at `lib/report_data.py:899`).

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && pixi run pytest tests/test_report_data.py -k xrefs -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Run the full existing report_data test suite (regression check)**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && pixi run pytest tests/test_report_data.py -v`
Expected: PASS — appending a field at the end of each positional list/row is
additive and does not renumber any existing field (each report's JS reads
fields by name via `DATA.fields.forEach`, not by fixed position — see Task 4).

- [ ] **Step 6: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
git add lib/report_data.py tests/test_report_data.py
git commit -m "$(cat <<'EOF'
Carry uniprot_xrefs through all three report payload builders

Adds 'xrefs' as the last ROW_FIELDS/CORE_ROW_FIELDS/LOSSES_ROW_FIELDS
entry in each of build_payload/build_core_payload/build_losses_payload,
stored inline (no _StringTable -- these values are 1:1 with the protein,
same as 'ec'/'af', nothing repeats enough to justify interning).

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: nf_NovInvenio — registry-driven renderer in `report_common.py`

**Files:**
- Modify: `/bigdata/stajichlab/jstajich/projects/NovInvenio/lib/report_common.py:539-620` (`externalLinksNode` and the new registry, placed just above it)
- Test: `/bigdata/stajichlab/jstajich/projects/NovInvenio/tests/test_report_js_behaviour.py` (its `MATRIX` fixture) and `/bigdata/stajichlab/jstajich/projects/NovInvenio/tests/js/drive_reports.mjs` (its `check(...)` assertions) — jsdom-driven; self-skips if `jsdom` isn't importable — confirmed not installed in this environment, so this test will not execute here. Write it anyway; it's real coverage wherever jsdom is present, e.g. CI.

**Interfaces:**
- Consumes: `o.xrefs` — a `DB:id|DB:id` string (Task 3's payload field, resolved by each template's `row[F.xrefs]`, wired in Task 5), and `o.proteome` (already existing, for the FungiDB-dedup check below).
- Produces: `externalLinksNode(o)` renders one link per recognized `DB` in `o.xrefs`, in addition to its existing links.

- [ ] **Step 1: Write the failing test**

`tests/test_report_js_behaviour.py` doesn't call into JS directly — it
generates real report HTML pages (via `bin/make_report.py` etc., from the
module's own `CONFIG`/`MATRIX` fixtures) into `fixture_dir`, then runs
`tests/js/drive_reports.mjs` as a **Node script** against those pages inside
jsdom; that script prints `PASS <name>`/`FAIL <name>` lines per assertion,
and `test_report_javascript_behaviour` asserts none of them are `FAIL`.
Extending coverage means: (a) add a `uniprot_xrefs` value to the existing
`MATRIX` fixture, (b) add `check(...)` calls to `drive_reports.mjs`.

`n1`'s row already has `SourceDB=fungidb` (`CONFIG`, `tests/test_report_js_behaviour.py`),
and its existing check at `drive_reports.mjs:108` asserts
`h.some((x) => x.includes('fungidb.org'))` — i.e. "at least one" FungiDB
link. Giving `n1` a `uniprot_xrefs` value makes this row the real,
end-to-end test of Task 4's FungiDB-dedup logic: after this change it must
render **exactly one** `fungidb.org` link, not two.

Edit `tests/test_report_js_behaviour.py`'s `MATRIX` constant: add a
`uniprot_xrefs` column (header) and give `n1` a value, others empty:

```python
MATRIX = (
    "protein_id\tsource_proteome\tNcra\tAfum\tDrome\tSpom\tScer\tgene_name\t"
    "product_description\tfunction_source\tBest_Swissprot\tPfam_Names\t"
    "Pfam_Accessions\tPfam_Evalues\tuniprot_xrefs\n"
    "n1\tNcra\t1\t1\t1\t0\t0\tada-1\tall development altered-1\tModelOrg_Ncra\t\t"
    "bZIP_1\tPF00170.27\t4.5e-09\tVEuPathDB:FungiDB:NCU10683|GeneID:5847462|KEGG:ncr:NCU10683|UnknownDB:xyz\n"
    "n2\tAfum\t1\t1\t1\t0\t0\t\t\t\t\t\t\t\t\n"
    "n3\tDrome\t0\t0\t1\t0\t0\t\t\t\t\t\t\t\t\n"
    "shared\tNcra\t1\t1\t1\t1\t1\t\tconserved thing\tPfam\t"
    "sp|P12345|TEST_YEAST Some protein\tAAA\tPF00004.31\t1e-20\t\n"
)
```

Edit `drive_reports.mjs`: change the existing loose "at least one" n1
FungiDB check at line 108 to an exact count, and add xref-specific checks
right after it (before the `n1: NCBI_TaxID` check, still inside the
`pick(rows, 'n1')` block that already set `h`/`b`):

```js
  const fungidbHrefs = h.filter((x) => x.includes('fungidb.org'));
  check('n1: xrefs-derived FungiDB link replaces genomeDbLink\'s, not both',
        fungidbHrefs.length === 1, fungidbHrefs.join(' '));
  check('n1: the surviving FungiDB link uses the xrefs gene ID, not the UniProt/protein_id',
        fungidbHrefs[0] && fungidbHrefs[0].includes('/gene/NCU10683'), fungidbHrefs[0]);
  check('n1: GeneID xref renders an NCBI Gene link',
        h.some((x) => x.includes('ncbi.nlm.nih.gov/gene/5847462')), h.join(' '));
  check('n1: KEGG xref renders a URL-encoded link (colon -> %3A)',
        h.some((x) => x.includes('ncr%3ANCU10683')), h.join(' '));
  check('n1: an unrecognized xref DB renders no link',
        !h.some((x) => x.includes('xyz')), h.join(' '));
```

- [ ] **Step 2: Run test to verify it fails (or confirm the skip)**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && pixi run pytest tests/test_report_js_behaviour.py -v`
Expected: either FAIL (Task 4 not implemented yet — `n1` still shows the old
single `genomeDbLink`-only FungiDB link with no xrefs-derived links, so the
new `check(...)` lines print `FAIL`, which `test_report_javascript_behaviour`
turns into a python-level assertion failure) if jsdom is installed, or
SKIPPED with a message about `jsdom` not being importable if it isn't
(confirmed not installed in this sandbox — `node -e "require.resolve('jsdom')"`
fails here). Either outcome is expected at this step; a SKIP is not a pass,
so don't treat it as Step 4's green bar.

- [ ] **Step 3: Implement the registry and the render loop**

In `lib/report_common.py`, add the registry just above `externalLinksNode`
(after the `// ---- the whole "External resources" block ----` comment,
before the `function externalLinksNode(o) {` line at `lib/report_common.py:553`):

```python
EXTERNAL_LINKS_JS = r"""
  // ---- generic UniProt DR-cross-reference linkouts ------------------------
  // o.xrefs is "DB:id|DB:id" (NII bin/extract_dat_annotations.py's DR-line
  // extraction, lib/report_data.py's 'xrefs' field) -- split each entry on
  // the FIRST colon only (never a fixed-arity split like goChipsNode's
  // go_ids parsing): KEGG ("ncr:NCU10683") and VEuPathDB ("FungiDB:NCU10683")
  // ids carry their own colons and must not be truncated.
  var VEUPATHDB_HOSTS = {
    FungiDB: "fungidb.org/fungidb",
    PlasmoDB: "plasmodb.org/plasmo",
    ToxoDB: "toxodb.org/toxo"
  };
  var XREF_LINK_TEMPLATES = {
    VEuPathDB: {
      render: function (id) {
        var i = id.indexOf(":");
        var project = i === -1 ? "VEuPathDB" : id.slice(0, i);
        var geneId = i === -1 ? id : id.slice(i + 1);
        var host = VEUPATHDB_HOSTS[project] || ("veupathdb.org/" + project.toLowerCase());
        return {
          label: project + " " + geneId,
          href: "https://" + host + "/app/record/gene/" + encodeURIComponent(geneId),
          title: project + " gene record for " + geneId
        };
      }
    },
    GeneID: {
      render: function (id) {
        return { label: "NCBI Gene " + id,
                 href: "https://www.ncbi.nlm.nih.gov/gene/" + encodeURIComponent(id),
                 title: "NCBI Gene " + id };
      }
    },
    RefSeq: {
      render: function (id) {
        return { label: "RefSeq " + id,
                 href: "https://www.ncbi.nlm.nih.gov/protein/" + encodeURIComponent(id),
                 title: "RefSeq protein " + id };
      }
    },
    KEGG: {
      render: function (id) {
        return { label: "KEGG " + id,
                 href: "https://www.genome.jp/dbget-bin/www_bget?" + encodeURIComponent(id),
                 title: "KEGG entry " + id };
      }
    },
    EnsemblFungi: {
      render: function (id) {
        return { label: "Ensembl Fungi " + id,
                 href: "https://fungi.ensembl.org/id/" + encodeURIComponent(id),
                 title: "Ensembl Fungi gene " + id };
      }
    },
    EnsemblBacteria: {
      render: function (id) {
        return { label: "Ensembl Bacteria " + id,
                 href: "https://bacteria.ensembl.org/id/" + encodeURIComponent(id),
                 title: "Ensembl Bacteria gene " + id };
      }
    }
  };
  // Returns {links: [<a> node, ...], hasVEuPathDB: bool} -- the caller (see
  // externalLinksNode below) needs to know whether a VEuPathDB xref fired so
  // it can suppress genomeDbLink's separate, UniProt-accession-keyed FungiDB
  // link for the same row (that one resolves against o.id, which for a
  // UniProt-sourced protein is the wrong ID space -- see the dedup comment in
  // externalLinksNode below).
  function xrefLinkNodes(xrefsStr) {
    var links = [];
    var hasVEuPathDB = false;
    (xrefsStr ? xrefsStr.split("|") : []).filter(Boolean).forEach(function (entry) {
      var i = entry.indexOf(":");
      if (i === -1) return;
      var db = entry.slice(0, i);
      var id = entry.slice(i + 1);
      var tmpl = XREF_LINK_TEMPLATES[db];
      if (!tmpl || !id) return;
      var r = tmpl.render(id);
      links.push(extLink(r.label, r.href, r.title));
      if (db === "VEuPathDB") hasVEuPathDB = true;
    });
    return { links: links, hasVEuPathDB: hasVEuPathDB };
  }
"""
```

Now edit `externalLinksNode` itself. Two changes:

1. Add the xref links (after the existing `var db = genomeDbLink(...)`
   block, before `// Model-org gene lookup ...`):

```js
    var xrefResult = xrefLinkNodes(o.xrefs);
    // FungiDB dedup: genomeDbLink's "fungidb" branch resolves against
    // geneIdFromProteinId(o.id), which for a UniProt-sourced protein is the
    // UniProt accession, not a real FungiDB gene ID -- that link is already
    // broken for these rows. When the xrefs-derived VEuPathDB link is
    // present, it is the correct one; suppress genomeDbLink's instead of
    // showing both (one dead) side by side.
    var suppressGenomeDbFungiDb = xrefResult.hasVEuPathDB && db &&
      /fungidb\.org/.test(db.href || "");
    if (db && !suppressGenomeDbFungiDb) links.appendChild(db);
    xrefResult.links.forEach(function (a) { links.appendChild(a); });
```

2. Remove the *old*, now-duplicate `if (db) links.appendChild(db);` line
   (the one immediately after `var db = genomeDbLink(...)` at
   `lib/report_common.py:564-565` in the current file) — it's replaced by the
   conditional version above.

Finally, register `EXTERNAL_LINKS_JS` wherever `LINKOUT_HELPERS_JS` is
concatenated into a page's `<script>` block (check `lib/report_template.py`,
`lib/core_report_template.py`, `lib/losses_report_template.py` for the
`LINKOUT_HELPERS_JS` import/concatenation site — add `EXTERNAL_LINKS_JS`
immediately alongside it in each).

- [ ] **Step 4: Run test to verify it passes (or still self-skips)**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && pixi run pytest tests/test_report_js_behaviour.py -v`
Expected: PASS (`test_report_javascript_behaviour` reports zero `FAIL` lines,
including the new xref checks) if jsdom is installed; SKIPPED (not FAILED)
otherwise. If it's FAILED rather than SKIPPED, jsdom-availability detection
broke — stop and check `tests/test_report_js_behaviour.py`'s skip condition
rather than proceeding.

- [ ] **Step 5: Run the full existing JS-behaviour and template test suites (regression check)**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && pixi run pytest tests/test_report_js_behaviour.py tests/test_report_templates.py -v`
Expected: PASS (or consistently SKIPPED for the jsdom-dependent file, same as
before this change — the point is nothing that used to pass now fails).

- [ ] **Step 6: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
git add lib/report_common.py lib/report_template.py lib/core_report_template.py lib/losses_report_template.py tests/test_report_js_behaviour.py tests/js/drive_reports.mjs
git commit -m "$(cat <<'EOF'
Render o.xrefs via a generic DB->URL-template registry

Replaces the implicit "one hardcoded link per DB" ceiling with
XREF_LINK_TEMPLATES + xrefLinkNodes(), covering FungiDB (via VEuPathDB's
Project:GeneID split), NCBI Gene, RefSeq, KEGG, and both Ensembl Fungi/
Bacteria. Suppresses genomeDbLink's separate FungiDB branch when a
UniProt-sourced row's xrefs already carry the correct VEuPathDB link,
since genomeDbLink's own FungiDB link resolves against the UniProt
accession (wrong ID space) for these rows and would otherwise show a
second, dead link next to the correct one.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: nf_NovInvenio — wire `row[F.xrefs]` into all three templates' `o` objects

**Files:**
- Modify: `/bigdata/stajichlab/jstajich/projects/NovInvenio/lib/report_template.py:1085-1094`
- Modify: `/bigdata/stajichlab/jstajich/projects/NovInvenio/lib/core_report_template.py:444-453`
- Modify: `/bigdata/stajichlab/jstajich/projects/NovInvenio/lib/losses_report_template.py:524-533`

**Interfaces:**
- Consumes: `row[F.xrefs]` (Task 3's field, resolved via each template's own `F` name→index map built by `DATA.fields.forEach`).
- Produces: `externalLinksNode`'s `o.xrefs` (Task 4) is populated for real report pages.

- [ ] **Step 1: Edit `report_template.py`**

At `lib/report_template.py:1085-1094`, add `xrefs: row[F.xrefs],` inside the
`externalLinksNode({...})` call, e.g. right after `geneUrl: row[F.gene_url],`:

```js
    detailEl.appendChild(field("External resources", externalLinksNode({
      id: row[F.id],
      gene: row[F.gene],
      sprot: row[F.sprot],
      geneUrl: row[F.gene_url],
      xrefs: row[F.xrefs],
      pfam: row[F.pfam_n],
      fsrcName: row[F.fsrc] >= 0 ? DATA.fsources[row[F.fsrc]] : "",
      seq: seq,
      proteome: row[F.src] >= 0 ? PROTEOMES[row[F.src]] : null
    })));
```

- [ ] **Step 2: Edit `core_report_template.py`**

Same edit at `lib/core_report_template.py:444-453`.

- [ ] **Step 3: Edit `losses_report_template.py`**

Same edit at `lib/losses_report_template.py:524-533`.

- [ ] **Step 4: Run the template + report_data test suites (regression check)**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio && pixi run pytest tests/test_report_templates.py tests/test_report_data.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
git add lib/report_template.py lib/core_report_template.py lib/losses_report_template.py
git commit -m "$(cat <<'EOF'
Wire row xrefs field into all three report templates' externalLinksNode calls

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Test rebuild — regenerate `pezizo_set1`'s annotations and re-sync its reports

This is a manual verification pass against the already-computed
`pezizo_set1` study (no fresh Nextflow run needed — `sync_reports.sh` is
explicitly designed to be safely re-run standalone whenever its inputs
change). No code changes in this task.

**Files touched (regenerated data, not source):**
- `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/pezizo_set1/annotations/*.tsv` (11 files, one per UniProt proteome referenced by this study's config)
- `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/results/pezizo_set1/presence_matrix.uniprot.tsv`, `novelties.html`, `core.html`, `losses.html`

- [ ] **Step 1: Confirm which proteomes this study references**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations && ls studies/fungi/pezizo_set1/annotations/`
Expected: 11 `UP0*.tsv` files (`UP000001261`, `UP000001805` [Ncra], `UP000001861`, `UP000002311`, `UP000002485`, `UP000002530`, `UP000008062`, `UP000010091`, `UP000014254`, `UP000186594`, `UP001658139`).

- [ ] **Step 2: Regenerate all 11 annotation TSVs with the new `xrefs` column**

Run (repeat per proteome, or loop):

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
for tsv in studies/fungi/pezizo_set1/annotations/UP*.tsv; do
  proteome=$(basename "$tsv" .tsv)
  dat="data/uniprot/${proteome}"/*_*.dat.gz
  pixi run python bin/extract_dat_annotations.py --dat-gz $dat --output "$tsv"
done
```

Expected: each command prints `Wrote studies/fungi/pezizo_set1/annotations/UP*.tsv (N proteins)` to stderr, N matching the prior file's row count (re-generating in place must not change protein counts, only add the `xrefs` column).

- [ ] **Step 3: Confirm the Ncra proteome's `xrefs` column looks right**

Run: `grep A7UWL5 studies/fungi/pezizo_set1/annotations/UP000001805.tsv | cut -f11`
(column 11 = `xrefs`, the last column per Task 1's fieldnames order)
Expected: a `|`-joined string containing `VEuPathDB:FungiDB:NCU10683`, `GeneID:5847462`, `RefSeq:XP_001728253.1`, `KEGG:ncr:NCU10683` — **not** `VEuPathDB:NCU10683` (the rev.-1 bug this plan's Task 1 fixed).

- [ ] **Step 4: Re-run `sync_reports.sh` against a local nf_NovInvenio checkout**

Run:

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
NOVINVENIO_ROOT=/bigdata/stajichlab/jstajich/projects/NovInvenio bin/sync_reports.sh fungi/pezizo_set1
```

Expected: stderr shows `== merging UniProt annotation into .../results/pezizo_set1 presence matrices ==`,
then `== regenerating novelties.html (--online) ==`, `== regenerating core.html ==`,
and (if a loss matrix exists for this study) `== regenerating losses.html (--online) ==`,
each completing without a traceback.

- [ ] **Step 5: Visually confirm one working FungiDB link (not two, not dead)**

Open `results/pezizo_set1/novelties.html` in a browser (`file://` path — it's
built `--online` per the script's own note, but still opens fine locally),
find a novelty candidate sourced from Ncra (e.g. search for `NCU10683` or any
Ncra-sourced row in the table), open its detail panel, and check the
"External resources" block:

- Exactly one FungiDB link (label starts with `FungiDB `), and it opens to a
  real FungiDB gene page (`https://fungidb.org/fungidb/app/record/gene/NCU#####`)
  — not a UniProt-accession-shaped dead link.
- NCBI Gene, RefSeq, and (if present for that protein) KEGG links alongside
  it, each opening to the right external record.

If two FungiDB links appear, or the FungiDB link 404s, stop — Task 4's dedup
logic or Task 1's VEuPathDB parsing has a bug; do not proceed to a "deploy"
step (copying into `docs/`, committing, or pushing) until this is clean.

- [ ] **Step 6: Report back, no commit in this task**

This task only regenerates already-gitignored `results/` output and
study-local `annotations/*.tsv` (tracked, but this is a routine
re-extraction, not a new/changed *source*) — nothing here needs a plan-step
commit. If the visual check in Step 5 passes, tell the user what you
confirmed and that `results/pezizo_set1/{novelties,core,losses}.html` are
ready for them to review; **do not** copy into `docs/` or push anything
without their explicit go-ahead (that publish step touches
`sync_reports.sh`'s already-documented `docs/` gallery, which reaches other
people once pushed).
