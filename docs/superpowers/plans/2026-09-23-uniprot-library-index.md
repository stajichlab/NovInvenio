# UniProt Library Index Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Every protein in a study that can be traced to a UniProt record (by accession, RefSeq ID or exact sequence) gets that record's links and annotation in all three reports, from the pipeline alone.

**Architecture:** A one-time `-entry UNIPROT_INDEX` workflow parses a pre-downloaded UniProt library into per-proteome records tables (`.tsv.zst`) plus one SQLite lookup index. A per-run `UNIPROT_LINK` step matches each config proteome against the index and writes `<Short>.uniprot_link.tsv`. That feeds the existing `ANNOTATE_MATRIX --uniprot_xref_files` path in place of `UNIPROT_XREF`. The report payload and the shared JS link builder gain the new fields.

**Tech Stack:** Nextflow DSL2, Python 3 (stdlib `sqlite3`, `hashlib`, `csv`), `zstd` CLI via `lib/compressed_io.py`, pytest, jsdom (via the existing `tests/js/drive_reports.mjs`).

**Spec:** `docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md` (read it before starting; this plan argues from it).

## Global Constraints

- Every `bin/` script: `#!/usr/bin/env python3`, `chmod +x`, argparse flags only (no positionals), and imports `lib/` via `sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))`.
- `lib/` modules have no I/O side effects at import time.
- New Nextflow files use `nextflow.enable.dsl=2`, one process per module file, and a `label` (`low_cpu` / `med_cpu` / `high_cpu`). Workflows declare `take:`/`main:`/`emit:`. The index outputs use `storeDir`. Never use `storeDir` and `publishDir` on the same output.
- Never use `$(dirname "${BASH_SOURCE[0]}")` in any script run by SLURM.
- Compressed output: the records tables are `.tsv.zst` (internal intermediate). Readers accept plain, `.gz` or `.zst` (`lib/compressed_io.open_maybe_compressed`).
- Report pages stay self-contained. Untrusted strings (titles, IDs) are inserted with `textContent` via `el()`. URLs come only from allow-listed templates, with IDs `encodeURIComponent`-escaped.
- Index format version constant: `FORMAT_VERSION = 1`.
- Record columns (exact order): `accession, entry_name, reviewed, protein_existence, taxon_id, gene_name, description, ec_numbers, go_ids, pfam_ids, pfam_names, interpro_ids, alphafold_id, xrefs, pubs, seq_md5, seq_len`.
- Link output columns (exact order): `protein_id, uniprot_accession, uniprot_gene_name, uniprot_description, uniprot_go_ids, uniprot_pfam_ids, uniprot_pfam_names, uniprot_interpro_ids, uniprot_ec_numbers, uniprot_alphafold_id, uniprot_xrefs, uniprot_match, uniprot_match_species, uniprot_reviewed, uniprot_pubs, uniprot_n_matches`.
- Match types (exact strings): `id`, `refseq`, `seq_own`, `seq_other`.
- Library under test: `/bigdata/stajichlab/shared/db/Uniprot/Fungi_2026_03` (release 2026_03, 1526 proteomes, CSV `fungi_proteomes_2026_03.csv`).
- Repo workflow (CLAUDE.md): one GitHub issue per phase, a branch per issue, PRs with `Closes #N`. Confirm with the user before creating issues or pushing.
- All work happens in the worktree `/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/uniprot-library-index`, never the shared `~/projects/NI` checkout.

## Review Focus

1. **Title text that contains `|` or `;`.** Publication titles often contain semicolons, and occasionally pipes. The `pubs` packing must round-trip a title with `;` intact and must never split one reference into two. The test is in Task 1.
2. **A protein ID with a version suffix, or a bare accession with an isoform suffix** (`Q7S6W2-2`, `sp|P12345-3|...`). Accession matching must strip the isoform suffix before lookup, not silently miss. The test is in Task 5.
3. **A config row with no `NCBI_TaxID` whose species name differs from UniProt's** (UniProt names carry strain text: "Neurospora crassa (strain ATCC 24698 / 74-OR23-1A / ...)"). The own-taxid lookup must match on the leading binomial, not the full string, or every protein falls to `seq_other`. The test is in Task 5.
4. **An index directory that exists but is incomplete** (the build was interrupted, so `seq_index.sqlite` is present but `manifest.json` is not). This must be a hard error naming the missing file, never an empty-link run. The test is in Task 4.
5. **A study proteome with lower-case or `*`-terminated sequences** (some FASTA exports end each protein with `*`). The sequence MD5 must be computed on the same normalized form on both sides (upper case, no whitespace, trailing `*` removed), or exact matches are silently lost. The tests are in Tasks 1 and 5.

---

## File Structure

| File | Status | Responsibility |
|---|---|---|
| `lib/uniprot_dat.py` | modify | Parser: add entry_name, reviewed, PE, pubs, seq_md5/len, the wider xref allow-list, skip stats; `normalize_seq()`; `RECORD_COLUMNS` |
| `lib/uniprot_index.py` | create | Index format constant, `UniProtIndex` reader (validate, lookups, record fetch), `IndexFormatError` |
| `bin/uniprot_plan_chunks.py` | create | Group library proteomes into parse batches by compressed size |
| `bin/uniprot_parse_dat.py` | create | Parse one batch → `records/<UP>.records.tsv.zst` + `stats/<UP>.json` |
| `bin/uniprot_build_index.py` | create | Merge records → `seq_index.sqlite` + `manifest.json` |
| `bin/uniprot_link.py` | create | Per-proteome match → `<Short>.uniprot_link.tsv` |
| `modules/uniprot_index.nf` | create | `UNIPROT_PLAN_CHUNKS`, `UNIPROT_PARSE_CHUNK`, `UNIPROT_SEQ_INDEX` processes (one file, since they share the storeDir contract; see Task 7 note) |
| `workflows/uniprot_index.nf` | create | `UNIPROT_INDEX_BUILD` workflow |
| `modules/uniprot_link.nf` | create | `UNIPROT_LINK` process |
| `main.nf` | modify | Named `workflow UNIPROT_INDEX`; replace `UNIPROT_XREF` wiring with `UNIPROT_LINK` |
| `nextflow.config` | modify | `uniprot_index`, `uniprot_library`, `uniprot_library_csv`, `uniprot_index_chunk_gb` params |
| `bin/annotate_presence_matrix.py` | modify | Add the 5 new `uniprot_*` columns |
| `lib/report_data.py` | modify | Add row fields `uacc, umatch, usp, urev, pubs` + `pub_sets` to all 3 payloads |
| `lib/report_common.py` | modify | Own-accession UniProt/AlphaFold links, new xref templates, `seq_other` suppression, publications node |
| `lib/report_template.py`, `lib/core_report_template.py`, `lib/losses_report_template.py` | modify | Pass the new fields to `externalLinksNode`; render Publications |
| `modules/uniprot_xref.nf`, `bin/build_uniprot_refseq_xref.py`, their tests | delete | Replaced by `UNIPROT_LINK` |
| `tests/test_uniprot_dat.py` | modify | New field tests |
| `tests/test_uniprot_index.py`, `tests/test_uniprot_link.py`, `tests/test_uniprot_parse_build.py` | create | Unit tests |
| `tests/data/uniprot_library/` | create | Tiny fixture library (2 proteomes) for tests and `-profile test` |
| `tests/test_report_data.py`, `tests/test_report_js_behaviour.py`, `tests/js/drive_reports.mjs`, `tests/test_annotate_presence_matrix.py` | modify | New fields and link behaviour |
| `CLAUDE.md`, `README.md` | modify | Document the index, the step and the params; remove `UniProtDatGz` docs except the deprecation note |

---

### Task 0: Tickets and setup

**Files:** none (GitHub + environment).

- [ ] **Step 1: Confirm with the user, then create issues.** One issue per phase, each linking the spec:
  - A: "UniProt library index: parser fields + index build (-entry UNIPROT_INDEX)" (Tasks 1–4, 7)
  - B: "UNIPROT_LINK: match every proteome to UniProt (replaces UNIPROT_XREF)" (Tasks 5, 6, 8)
  - C: "Reports: own-accession UniProt/AlphaFold links, wider xrefs, publications" (Tasks 9–11)
  - D: "Build the Fungi_2026_03 index and verify on sordariales_shallow" (Task 12)

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/uniprot-library-index
gh issue create --title "UniProt library index: parser fields + index build (-entry UNIPROT_INDEX)" \
  --body "Implements Components 1–2 of docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md. Plan: docs/superpowers/plans/2026-09-23-uniprot-library-index.md, Tasks 1–4, 7."
# repeat for B, C, D with their component/task references
```

- [ ] **Step 2: Branch.** The worktree is on `uniprot-library-index-spec` (holds the spec and this plan). Create the phase-A branch from it:

```bash
git switch -c <A#>-uniprot-library-index
```

  Later phases branch from the previous phase branch, so each PR stacks. Retarget the base after each merge.

- [ ] **Step 3: Environment.**

```bash
pixi install
pixi run pytest tests/test_uniprot_dat.py -q      # baseline: all pass
```

---

### Task 1: Parser fields in `lib/uniprot_dat.py`

**Files:**
- Modify: `lib/uniprot_dat.py`
- Test: `tests/test_uniprot_dat.py`

**Interfaces:**
- Produces:
  - `RECORD_COLUMNS: tuple[str, ...]` (exact order in Global Constraints)
  - `normalize_seq(seq: str) -> str`
  - `seq_md5(seq: str) -> str` (hex MD5 of `normalize_seq(seq)`)
  - `parse_dat_gz(path, stats: dict | None = None)` now yields dicts with every key in `RECORD_COLUMNS`. When `stats` is given, it is updated with `{'n_records': int, 'n_skipped_no_ac': int, 'n_skipped_no_seq': int}`.
  - `pubs` packing: entries `PMID;DOI;scope;title`, joined by `|`. `scope` is `proteome` or `protein`. Empty PMID/DOI are empty strings. `|` inside a title becomes `/`.
  - `XREF_FIRST_FIELD_DBS` gains: `PDB PANTHER OrthoDB STRING eggNOG Gene3D SUPFAM PROSITE SMART CDD PRINTS PIRSF HAMAP NCBIfam FunFam MEROPS CAZy ESTHER TCDB BRENDA UniPathway`.

- [ ] **Step 1: Write the failing tests.** Add to `tests/test_uniprot_dat.py`:
  - Replace `test_xrefs_drops_pdb_and_non_allowlisted_dbs`: PDB and STRING are now kept, and an unknown DB is still dropped.
  - Append two references and a real SQ block to the first fixture record.

```python
FIXTURE_DAT = """\
ID   A7UWL5_NEUCR            Unreviewed;       6 AA.
AC   A7UWL5;
DE   SubName: Full=NRS/ER;
GN   ORFNames=NCU10683;
OX   NCBI_TaxID=367110;
RN   [1]
RP   NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA].
RX   PubMed=12712197; DOI=10.1038/nature01554;
RA   Galagan J.E.;
RT   "The genome sequence of the filamentous fungus Neurospora crassa.";
RL   Nature 422:859-868(2003).
RN   [2]
RP   FUNCTION; SUBCELLULAR LOCATION.
RX   PubMed=99999999;
RT   "A title; with a semicolon | and a pipe
RT   that wraps.";
RL   J. Test 1:1-2(2020).
DR   VEuPathDB; FungiDB:NCU10683; -.
DR   GeneID; 5847462; -.
DR   RefSeq; XP_001728253.1; XM_001728201.2.
DR   RefSeq; XP_999999999.1; XM_999999999.1.
DR   KEGG; ncr:NCU10683; -.
DR   EnsemblFungi; EAA34304; EAA34304; NCU06768.
DR   PDB; 6YWS; EM; 2.74 A; K=1-249.
DR   STRING; 367110.A7UWL5; -.
DR   PANTHER; PTHR10000; PHOSPHOSERINE PHOSPHATASE; 1.
DR   UnknownDB; zzz; -.
DR   Pfam; PF04321; RmlD_sub_bind; 1.
DR   GO; GO:0016020; C:membrane; IEA:UniProtKB-KW.
PE   3: Inferred from homology;
SQ   SEQUENCE   6 AA;  1 MW;  0000000000000000 CRC64;
     MATTE I
//
ID   B0000000_TEST           Reviewed;         10 AA.
AC   B0000000;
GN   ORFNames=NCU99999;
OX   NCBI_TaxID=367110;
DE   RecName: Full=Test with Scer-style Ensembl field order;
DR   EnsemblFungi; YML051W_mRNA; YML051W; YML051W.
PE   1: Evidence at protein level;
SQ   SEQUENCE   6 AA;  1 MW;  0000000000000000 CRC64;
     MATTEI
//
ID   C0000000_TEST           Unreviewed;       6 AA.
AC   C0000000;
OX   NCBI_TaxID=367110;
DE   SubName: Full=No sequence block;
//
"""
```

```python
import hashlib
from uniprot_dat import RECORD_COLUMNS, normalize_seq, parse_dat_gz, seq_md5  # noqa: E402


def test_record_has_every_column_in_order(tmp_path):
    rec = next(iter(parse_dat_gz(_write_fixture(tmp_path))))
    assert tuple(rec.keys()) == RECORD_COLUMNS


def test_id_line_reviewed_entry_name_and_pe(tmp_path):
    recs = {r["accession"]: r for r in parse_dat_gz(_write_fixture(tmp_path))}
    assert recs["A7UWL5"]["entry_name"] == "A7UWL5_NEUCR"
    assert recs["A7UWL5"]["reviewed"] == "0"
    assert recs["B0000000"]["reviewed"] == "1"
    assert recs["A7UWL5"]["protein_existence"] == "3"
    assert recs["B0000000"]["protein_existence"] == "1"


def test_sequence_md5_and_length(tmp_path):
    recs = {r["accession"]: r for r in parse_dat_gz(_write_fixture(tmp_path))}
    assert recs["A7UWL5"]["seq_md5"] == hashlib.md5(b"MATTEI").hexdigest()
    assert recs["A7UWL5"]["seq_len"] == "6"


def test_normalize_seq_strips_case_whitespace_and_stop():
    assert normalize_seq(" matt\nei* ") == "MATTEI"
    assert seq_md5("mattei*") == hashlib.md5(b"MATTEI").hexdigest()


def test_record_without_sequence_is_skipped_and_counted(tmp_path):
    stats = {}
    recs = {r["accession"] for r in parse_dat_gz(_write_fixture(tmp_path), stats=stats)}
    assert "C0000000" not in recs
    assert stats == {"n_records": 2, "n_skipped_no_ac": 0, "n_skipped_no_seq": 1}


def test_pubs_scope_and_title_with_separators(tmp_path):
    recs = {r["accession"]: r for r in parse_dat_gz(_write_fixture(tmp_path))}
    pubs = recs["A7UWL5"]["pubs"].split("|")
    assert pubs[0] == ("12712197;10.1038/nature01554;proteome;"
                       "The genome sequence of the filamentous fungus Neurospora crassa.")
    pmid, doi, scope, title = pubs[1].split(";", 3)
    assert (pmid, doi, scope) == ("99999999", "", "protein")
    assert title == "A title; with a semicolon / and a pipe that wraps."
    assert recs["B0000000"]["pubs"] == ""


def test_xrefs_keep_new_allowlist_and_drop_unknown(tmp_path):
    recs = {r["accession"]: r for r in parse_dat_gz(_write_fixture(tmp_path))}
    x = recs["A7UWL5"]["xrefs"].split("|")
    assert "PDB:6YWS" in x
    assert "STRING:367110.A7UWL5" in x
    assert "PANTHER:PTHR10000" in x
    assert not any(e.startswith("UnknownDB") for e in x)
```

- [ ] **Step 2: Run the tests to verify they fail.**

Run: `pixi run pytest tests/test_uniprot_dat.py -v`
Expected: FAIL (`ImportError: cannot import name 'RECORD_COLUMNS'`).

- [ ] **Step 3: Implement.** In `lib/uniprot_dat.py`:
  - Add `import hashlib`.
  - Add the regexes, the constants and the two helpers.
  - Extend `parse_dat_gz`. Keep every existing rule unchanged.

```python
ID_RE = re.compile(r"^ID\s+(\S+)\s+(Reviewed|Unreviewed);")
PE_RE = re.compile(r"^PE\s+(\d)")
RX_PMID_RE = re.compile(r"PubMed=(\d+)")
RX_DOI_RE = re.compile(r"DOI=([^;\s]+)")
GENOME_RP = "LARGE SCALE GENOMIC DNA"

RECORD_COLUMNS = (
    "accession", "entry_name", "reviewed", "protein_existence", "taxon_id",
    "gene_name", "description", "ec_numbers", "go_ids", "pfam_ids", "pfam_names",
    "interpro_ids", "alphafold_id", "xrefs", "pubs", "seq_md5", "seq_len",
)

XREF_FIRST_FIELD_DBS = {
    "VEuPathDB", "GeneID", "KEGG", "RefSeq",
    "PDB", "PANTHER", "OrthoDB", "STRING", "eggNOG", "Gene3D", "SUPFAM", "PROSITE",
    "SMART", "CDD", "PRINTS", "PIRSF", "HAMAP", "NCBIfam", "FunFam", "MEROPS",
    "CAZy", "ESTHER", "TCDB", "BRENDA", "UniPathway",
}


def normalize_seq(seq: str) -> str:
    """Upper case, no whitespace, trailing stop '*' removed -- the one form both
    the index and the per-run lookup hash."""
    return "".join(seq.split()).upper().rstrip("*")


def seq_md5(seq: str) -> str:
    return hashlib.md5(normalize_seq(seq).encode("ascii", "replace")).hexdigest()
```

  State to add inside `parse_dat_gz`, reset in `reset()`: `entry_name`, `reviewed`, `pe`, `refs` (list of dicts), `cur_ref` (dict or None), `in_seq` (bool), `seq_chunks` (list).

  Line handling, added before the existing `AC` check:

```python
            if line.startswith("ID"):
                m = ID_RE.match(line)
                if m:
                    entry_name, reviewed = m.group(1), "1" if m.group(2) == "Reviewed" else "0"
                continue
            if in_seq:
                if line.startswith("     "):
                    seq_chunks.append(line)
                    continue
            if line.startswith("RN"):
                cur_ref = {"pmid": "", "doi": "", "rp": "", "rt": []}
                refs.append(cur_ref)
                continue
            if cur_ref is not None and line[:2] in ("RP", "RX", "RT", "RA", "RG", "RL", "RC"):
                body = line[5:].rstrip("\n")
                if line.startswith("RP"):
                    cur_ref["rp"] += " " + body
                elif line.startswith("RX"):
                    m = RX_PMID_RE.search(body)
                    if m:
                        cur_ref["pmid"] = m.group(1)
                    m = RX_DOI_RE.search(body)
                    if m:
                        cur_ref["doi"] = m.group(1)
                elif line.startswith("RT"):
                    cur_ref["rt"].append(body.strip())
                continue
            if line.startswith("PE"):
                m = PE_RE.match(line)
                if m:
                    pe = m.group(1)
                continue
            if line.startswith("SQ"):
                in_seq = True
                continue
```

  At `//`, build `pubs` and the sequence, then yield only when both `accession` and a non-empty sequence exist. Update `stats` otherwise:

```python
def _pack_ref(ref):
    title = " ".join(ref["rt"]).strip().strip('"').rstrip(";").rstrip('"').strip()
    title = title.replace("|", "/")
    scope = "proteome" if GENOME_RP in ref["rp"] else "protein"
    return f'{ref["pmid"]};{ref["doi"]};{scope};{title}'
```

```python
            if line.startswith("//"):
                seq = normalize_seq("".join(seq_chunks))
                if not accession:
                    if stats is not None:
                        stats["n_skipped_no_ac"] += 1
                elif not seq:
                    if stats is not None:
                        stats["n_skipped_no_seq"] += 1
                else:
                    if stats is not None:
                        stats["n_records"] += 1
                    pubs = [_pack_ref(r) for r in refs if r["pmid"] or r["doi"]]
                    yield {
                        "accession": accession,
                        "entry_name": entry_name or "",
                        "reviewed": reviewed or "0",
                        "protein_existence": pe or "",
                        "taxon_id": taxon_id or "",
                        "gene_name": gene_name or "",
                        "description": (rec_description or sub_description or "").strip(),
                        "ec_numbers": "|".join(ec_numbers),
                        "go_ids": "|".join(go_ids),
                        "pfam_ids": "|".join(pfam_ids),
                        "pfam_names": "|".join(pfam_names),
                        "interpro_ids": "|".join(interpro_ids),
                        "alphafold_id": alphafold_id or "",
                        "xrefs": "|".join(xrefs),
                        "pubs": "|".join(pubs),
                        "seq_md5": hashlib.md5(seq.encode("ascii", "replace")).hexdigest(),
                        "seq_len": str(len(seq)),
                    }
                reset()
                continue
```

  At the top of `parse_dat_gz`, initialise the stats: `if stats is not None: stats.update(n_records=0, n_skipped_no_ac=0, n_skipped_no_seq=0)`.
  In the existing `AC`/`GN`/`OX`/`DR`/`DE` handling, set `cur_ref = None` once those line types begin. Reference blocks end at the first non-R line. The `cur_ref` guard above already stops RP/RX/RT capture once `cur_ref` is `None`, so add `cur_ref = None` as the first statement inside each of the `DE`, `OX`, `DR`, `PE`, `SQ` and `CC` branches. For `CC`, add `if line.startswith("CC"): cur_ref = None; continue`.
  Replace `opener = gzip.open ...` with `from compressed_io import open_maybe_compressed` and `with open_maybe_compressed(path) as fh:`, so `.zst` fixtures work too.

- [ ] **Step 4: Run the tests to verify they pass.**

Run: `pixi run pytest tests/test_uniprot_dat.py -v`
Expected: all PASS (existing + 7 new).

- [ ] **Step 5: Check against real records.** Record the result in the commit message:

```bash
pixi run python -c "
import sys; sys.path.insert(0,'lib'); from uniprot_dat import parse_dat_gz
st={}; recs=list(parse_dat_gz('/bigdata/stajichlab/shared/db/Uniprot/Fungi_2026_03/data/UP000001805_367110.dat.gz', stats=st))
print(st); r=[x for x in recs if x['accession']=='Q7S6W2'][0]; print({k:r[k][:80] for k in r})"
```

Expected: `n_records` equals 9759 (the count from the current parser), and Q7S6W2 has `alphafold_id=Q7S6W2` plus a non-empty `seq_md5` and `pubs`.

- [ ] **Step 6: Commit.**

```bash
git add lib/uniprot_dat.py tests/test_uniprot_dat.py
git commit -m "uniprot_dat: entry name, reviewed, PE, publications, sequence MD5, wider xref allow-list (#A)"
```

---

### Task 2: Batch planning (`bin/uniprot_plan_chunks.py`)

**Files:**
- Create: `bin/uniprot_plan_chunks.py`
- Test: `tests/test_uniprot_parse_build.py`
- Create: `tests/data/uniprot_library/` fixture

**Interfaces:**
- Produces:
  - `plan_chunks(rows: list[dict], library: Path, chunk_bytes: int) -> list[list[dict]]`, where each row is a CSV dict with `proteome_id`, `tax_id`, `species_name`, `file_prefix`.
  - The CLI writes a TSV `chunk_id\tproteome_id\tfile_prefix` (chunk ids `chunk_0000`, …) and lists missing files on stderr.
- Fixture library layout, used by every later task:

```
tests/data/uniprot_library/
  README                      # contains a line "Release 2099_01, 01-Jan-2099"
  proteomes.csv               # proteome_id,tax_id,species_name,oscode,n_canonical,n_isoform,file_prefix,lineage
  data/UP000000001_367110.dat.gz   # 2 records: Q7S6W2-like (TrEMBL) + one Swiss-Prot
  data/UP000000002_559292.dat.gz   # 1 record whose sequence is identical to one in UP000000001
  # proteomes.csv also lists UP000000003 with no .dat.gz (missing-file case)
```

- [ ] **Step 1: Create the fixture.** Write `tests/data/uniprot_library/make_fixture.py`, a script that writes the files above from inline record text. The records use the Task 1 fixture format with sequences `MKVLLAQ` (UP1 rec 1, accession `Q7S6W2`, RefSeq `XP_960565.1`), `MSEQTWO` (UP1 rec 2, Swiss-Prot `P00001`) and `MKVLLAQ` (UP2, accession `P00002`, Reviewed). Species names are `"Neurospora crassa (strain ATCC 24698 / 74-OR23-1A)"`, `"Saccharomyces cerevisiae (strain ATCC 204508 / S288c)"` and `"Missing fungus"`. Run it once and commit its outputs together with the script.

- [ ] **Step 2: Write the failing test.**

```python
import csv
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
LIB = REPO / "tests" / "data" / "uniprot_library"
sys.path.insert(0, str(REPO / "bin"))
import uniprot_plan_chunks as upc  # noqa: E402


def _rows():
    with open(LIB / "proteomes.csv") as fh:
        return list(csv.DictReader(fh))


def test_plan_chunks_groups_by_size_and_skips_missing():
    chunks = upc.plan_chunks(_rows(), LIB, chunk_bytes=1)   # 1 byte -> one proteome per chunk
    ids = [[r["proteome_id"] for r in c] for c in chunks]
    assert ids == [["UP000000001"], ["UP000000002"]]
    one = upc.plan_chunks(_rows(), LIB, chunk_bytes=10**9)
    assert len(one) == 1 and len(one[0]) == 2


def test_cli_writes_tsv_and_reports_missing(tmp_path):
    out = tmp_path / "chunks.tsv"
    p = subprocess.run([sys.executable, str(REPO / "bin" / "uniprot_plan_chunks.py"),
                        "--library", str(LIB), "--library-csv", "proteomes.csv",
                        "--chunk-gb", "8", "--output", str(out)],
                       capture_output=True, text=True, check=True)
    assert "UP000000003" in p.stderr
    lines = out.read_text().splitlines()
    assert lines[0] == "chunk_id\tproteome_id\tfile_prefix"
    assert lines[1].startswith("chunk_0000\tUP000000001\t")
```

- [ ] **Step 3: Run it to verify it fails.** Run: `pixi run pytest tests/test_uniprot_parse_build.py -v`. Expected: FAIL (module not found).

- [ ] **Step 4: Implement `bin/uniprot_plan_chunks.py`.**

```python
#!/usr/bin/env python3
"""Group a UniProt library's proteomes into parse batches of about --chunk-gb of
compressed .dat.gz each (toward 1-1.5 h per UNIPROT_PARSE_CHUNK job; the parser
reads ~6.4 MB/s, measured 2026-09-23). Proteomes listed in the CSV with no
.dat.gz are reported on stderr and left out."""
import argparse
import csv
import sys
from pathlib import Path


def plan_chunks(rows, library, chunk_bytes):
    chunks, cur, cur_size = [], [], 0
    for r in rows:
        f = Path(library) / "data" / f"{r['file_prefix']}.dat.gz"
        if not f.exists():
            continue
        size = f.stat().st_size
        if cur and cur_size + size > chunk_bytes:
            chunks.append(cur)
            cur, cur_size = [], 0
        cur.append(r)
        cur_size += size
    if cur:
        chunks.append(cur)
    return chunks


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--library", required=True, type=Path)
    ap.add_argument("--library-csv", required=True, help="CSV file name inside --library")
    ap.add_argument("--chunk-gb", type=float, default=8.0)
    ap.add_argument("--output", required=True, type=Path)
    a = ap.parse_args()
    with open(a.library / a.library_csv) as fh:
        rows = list(csv.DictReader(fh))
    missing = [r["proteome_id"] for r in rows
               if not (a.library / "data" / f"{r['file_prefix']}.dat.gz").exists()]
    for pid in missing:
        print(f"WARNING: no .dat.gz for {pid}; skipped", file=sys.stderr)
    chunks = plan_chunks(rows, a.library, int(a.chunk_gb * 1e9))
    with open(a.output, "w") as out:
        out.write("chunk_id\tproteome_id\tfile_prefix\n")
        for i, c in enumerate(chunks):
            for r in c:
                out.write(f"chunk_{i:04d}\t{r['proteome_id']}\t{r['file_prefix']}\n")
    print(f"{len(chunks)} chunks, {sum(map(len, chunks))} proteomes, "
          f"{len(missing)} missing", file=sys.stderr)


if __name__ == "__main__":
    main()
```

- [ ] **Step 5: Run the tests to verify they pass, then commit.**

```bash
chmod +x bin/uniprot_plan_chunks.py
pixi run pytest tests/test_uniprot_parse_build.py -v    # PASS
git add bin/uniprot_plan_chunks.py tests/test_uniprot_parse_build.py tests/data/uniprot_library
git commit -m "uniprot_plan_chunks: size-based parse batches over a UniProt library (#A)"
```

---

### Task 3: Batch parse (`bin/uniprot_parse_dat.py`)

**Files:**
- Create: `bin/uniprot_parse_dat.py`
- Test: `tests/test_uniprot_parse_build.py` (append)

**Interfaces:**
- Consumes: `parse_dat_gz(path, stats)` and `RECORD_COLUMNS` (Task 1); the chunks TSV (Task 2).
- Produces:
  - `records/<proteome_id>.records.tsv.zst` (header = `RECORD_COLUMNS`);
  - `stats/<proteome_id>.json` = `{"proteome_id", "file_prefix", "dat_sha256", "n_records", "n_skipped_no_ac", "n_skipped_no_seq"}`.

- [ ] **Step 1: Write the failing test.**

```python
import json
from compressed_io import open_maybe_compressed  # after sys.path insert of REPO/'lib'


def test_parse_chunk_writes_records_and_stats(tmp_path):
    chunks = tmp_path / "chunks.tsv"
    chunks.write_text("chunk_id\tproteome_id\tfile_prefix\n"
                      "chunk_0000\tUP000000001\tUP000000001_367110\n")
    subprocess.run([sys.executable, str(REPO / "bin" / "uniprot_parse_dat.py"),
                    "--library", str(LIB), "--chunks", str(chunks),
                    "--chunk-id", "chunk_0000", "--outdir", str(tmp_path)], check=True)
    with open_maybe_compressed(tmp_path / "records" / "UP000000001.records.tsv.zst") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert [r["accession"] for r in rows] == ["Q7S6W2", "P00001"]
    st = json.loads((tmp_path / "stats" / "UP000000001.json").read_text())
    assert st["n_records"] == 2 and len(st["dat_sha256"]) == 64
```

- [ ] **Step 2: Run it to verify it fails.** Expected: FAIL (script missing).

- [ ] **Step 3: Implement.**

```python
#!/usr/bin/env python3
"""Parse one batch (--chunk-id) of a UniProt library's .dat.gz files into
records/<proteome_id>.records.tsv.zst + stats/<proteome_id>.json. See
docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import open_maybe_compressed_write  # noqa: E402
from uniprot_dat import RECORD_COLUMNS, parse_dat_gz  # noqa: E402


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--library", required=True, type=Path)
    ap.add_argument("--chunks", required=True, type=Path)
    ap.add_argument("--chunk-id", required=True)
    ap.add_argument("--outdir", required=True, type=Path)
    a = ap.parse_args()
    (a.outdir / "records").mkdir(parents=True, exist_ok=True)
    (a.outdir / "stats").mkdir(parents=True, exist_ok=True)
    with open(a.chunks) as fh:
        todo = [r for r in csv.DictReader(fh, delimiter="\t") if r["chunk_id"] == a.chunk_id]
    if not todo:
        sys.exit(f"ERROR: no proteomes for {a.chunk_id} in {a.chunks}")
    for r in todo:
        dat = a.library / "data" / f"{r['file_prefix']}.dat.gz"
        stats = {}
        out = a.outdir / "records" / f"{r['proteome_id']}.records.tsv.zst"
        with open_maybe_compressed_write(out) as fh:
            w = csv.DictWriter(fh, fieldnames=RECORD_COLUMNS, delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            for rec in parse_dat_gz(dat, stats=stats):
                w.writerow(rec)
        stats.update(proteome_id=r["proteome_id"], file_prefix=r["file_prefix"],
                     dat_sha256=sha256(dat))
        (a.outdir / "stats" / f"{r['proteome_id']}.json").write_text(json.dumps(stats))
        print(f"{r['proteome_id']}: {stats['n_records']} records "
              f"({stats['n_skipped_no_ac']} no AC, {stats['n_skipped_no_seq']} no SQ)",
              file=sys.stderr)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the tests to verify they pass, then commit.**

```bash
chmod +x bin/uniprot_parse_dat.py
pixi run pytest tests/test_uniprot_parse_build.py -v
git add bin/uniprot_parse_dat.py tests/test_uniprot_parse_build.py
git commit -m "uniprot_parse_dat: per-proteome records tables + stats for one batch (#A)"
```

---

### Task 4: Index build and reader (`bin/uniprot_build_index.py`, `lib/uniprot_index.py`)

**Files:**
- Create: `bin/uniprot_build_index.py`, `lib/uniprot_index.py`
- Test: `tests/test_uniprot_index.py`

**Interfaces:**
- Consumes: the `records/` and `stats/` outputs (Task 3), the library CSV and README.
- Produces, in `lib/uniprot_index.py`:
  - `FORMAT_VERSION = 1`
  - `class IndexFormatError(Exception)`
  - `class UniProtIndex`:
    - `UniProtIndex(index_dir)`: validates `manifest.json`, `seq_index.sqlite` and `records/`, and checks `format_version == FORMAT_VERSION`, else raises `IndexFormatError(msg)` naming the problem;
    - `.by_accession(acc) -> dict | None` with keys `accession, proteome_id, taxid`;
    - `.by_refseq(refseq_id) -> list[str]` (accessions);
    - `.by_md5(md5) -> list[dict]` with keys `accession, proteome_id, taxid, reviewed`;
    - `.taxid_for_species(name) -> int | None`: matches the first two words (the binomial) of `name` against the first two words of `proteome.species_name`, case-insensitive; returns the taxid only when exactly one taxid matches;
    - `.species_name(proteome_id) -> str`;
    - `.records(proteome_id, accessions: set[str]) -> dict[str, dict]`: rows from `records/<proteome_id>.records.tsv.zst`.
  - `build_index(records_dir, stats_dir, library, library_csv, out_dir)` writes `seq_index.sqlite` + `manifest.json` and returns the manifest dict.
- `manifest.json` keys: `format_version, library, release, library_csv, library_csv_sha256, n_proteomes, n_records, missing, proteomes (list of stats dicts), build_date`.

- [ ] **Step 1: Write the failing tests.**

```python
import csv
import json
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
LIB = REPO / "tests" / "data" / "uniprot_library"
sys.path.insert(0, str(REPO / "lib"))
from uniprot_index import FORMAT_VERSION, IndexFormatError, UniProtIndex  # noqa: E402


def build_fixture_index(tmp_path):
    """Plan, parse and build the fixture library into tmp_path/index."""
    work = tmp_path / "work"
    run = lambda *a: subprocess.run([sys.executable, *map(str, a)], check=True)
    run(REPO / "bin" / "uniprot_plan_chunks.py", "--library", LIB,
        "--library-csv", "proteomes.csv", "--chunk-gb", "8", "--output", tmp_path / "chunks.tsv")
    run(REPO / "bin" / "uniprot_parse_dat.py", "--library", LIB, "--chunks",
        tmp_path / "chunks.tsv", "--chunk-id", "chunk_0000", "--outdir", work)
    run(REPO / "bin" / "uniprot_build_index.py", "--records-dir", work / "records",
        "--stats-dir", work / "stats", "--library", LIB, "--library-csv", "proteomes.csv",
        "--output-dir", tmp_path / "index")
    return tmp_path / "index"


def test_manifest_and_tables(tmp_path):
    idx_dir = build_fixture_index(tmp_path)
    m = json.loads((idx_dir / "manifest.json").read_text())
    assert m["format_version"] == FORMAT_VERSION
    assert m["release"] == "2099_01"
    assert m["n_proteomes"] == 2 and m["n_records"] == 3
    assert m["missing"] == ["UP000000003"]
    con = sqlite3.connect(idx_dir / "seq_index.sqlite")
    assert con.execute("select count(*) from seq").fetchone()[0] == 3
    assert con.execute("select accession from refseq where refseq_id='XP_960565.1'").fetchone()[0] == "Q7S6W2"


def test_reader_lookups(tmp_path):
    idx = UniProtIndex(build_fixture_index(tmp_path))
    assert idx.by_accession("Q7S6W2")["proteome_id"] == "UP000000001"
    assert idx.by_refseq("XP_960565.1") == ["Q7S6W2"]
    import hashlib
    hits = idx.by_md5(hashlib.md5(b"MKVLLAQ").hexdigest())
    assert sorted(h["accession"] for h in hits) == ["P00002", "Q7S6W2"]
    assert idx.taxid_for_species("Neurospora crassa") == 367110
    assert idx.taxid_for_species("Neurospora crassa OR74A") == 367110
    assert idx.taxid_for_species("Nonexistent fungus") is None
    recs = idx.records("UP000000001", {"Q7S6W2"})
    assert recs["Q7S6W2"]["alphafold_id"] == "Q7S6W2"


def test_incomplete_index_is_a_hard_error(tmp_path):
    idx_dir = build_fixture_index(tmp_path)
    (idx_dir / "manifest.json").unlink()
    with pytest.raises(IndexFormatError, match="manifest.json"):
        UniProtIndex(idx_dir)


def test_wrong_format_version_is_a_hard_error(tmp_path):
    idx_dir = build_fixture_index(tmp_path)
    m = json.loads((idx_dir / "manifest.json").read_text())
    m["format_version"] = 99
    (idx_dir / "manifest.json").write_text(json.dumps(m))
    with pytest.raises(IndexFormatError, match="format_version 99"):
        UniProtIndex(idx_dir)
```

- [ ] **Step 2: Run them to verify they fail.** Run: `pixi run pytest tests/test_uniprot_index.py -v`. Expected: FAIL (module missing).

- [ ] **Step 3: Implement `lib/uniprot_index.py`.**

```python
"""Reader for a UniProt library index built by bin/uniprot_build_index.py
(docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md)."""
import csv
import json
import sqlite3
from pathlib import Path

from compressed_io import open_maybe_compressed

FORMAT_VERSION = 1


class IndexFormatError(Exception):
    pass


def _binomial(name):
    return " ".join((name or "").lower().replace("(", " ").split()[:2])


class UniProtIndex:
    def __init__(self, index_dir):
        self.dir = Path(index_dir)
        for need in ("manifest.json", "seq_index.sqlite", "records"):
            if not (self.dir / need).exists():
                raise IndexFormatError(f"UniProt index {self.dir} is incomplete: {need} missing")
        self.manifest = json.loads((self.dir / "manifest.json").read_text())
        v = self.manifest.get("format_version")
        if v != FORMAT_VERSION:
            raise IndexFormatError(f"UniProt index {self.dir} has format_version {v}; "
                                   f"this code reads {FORMAT_VERSION}")
        uri = f"file:{self.dir / 'seq_index.sqlite'}?mode=ro"
        self.con = sqlite3.connect(uri, uri=True)
        self._species = {pid: (tax, name) for pid, tax, name in
                         self.con.execute("select proteome_id, taxid, species_name from proteome")}

    def by_accession(self, acc):
        r = self.con.execute("select accession, proteome_id, taxid from acc where accession=?",
                             (acc,)).fetchone()
        return dict(zip(("accession", "proteome_id", "taxid"), r)) if r else None

    def by_refseq(self, refseq_id):
        return [r[0] for r in self.con.execute(
            "select accession from refseq where refseq_id=? order by accession", (refseq_id,))]

    def by_md5(self, md5):
        return [dict(zip(("accession", "proteome_id", "taxid", "reviewed"), r)) for r in
                self.con.execute("select accession, proteome_id, taxid, reviewed from seq "
                                 "where md5=?", (md5,))]

    def taxid_for_species(self, name):
        want = _binomial(name)
        taxids = {tax for tax, sp in self._species.values() if _binomial(sp) == want}
        return taxids.pop() if len(taxids) == 1 else None

    def species_name(self, proteome_id):
        return self._species.get(proteome_id, (None, ""))[1]

    def records(self, proteome_id, accessions):
        out = {}
        path = self.dir / "records" / f"{proteome_id}.records.tsv.zst"
        with open_maybe_compressed(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["accession"] in accessions:
                    out[row["accession"]] = row
        return out
```

- [ ] **Step 4: Implement `bin/uniprot_build_index.py`.**

```python
#!/usr/bin/env python3
"""Merge records/<UP>.records.tsv.zst + stats/<UP>.json into seq_index.sqlite +
manifest.json (UniProt library index, format 1)."""
import argparse
import csv
import datetime
import hashlib
import json
import re
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import open_maybe_compressed  # noqa: E402
from uniprot_index import FORMAT_VERSION  # noqa: E402

RELEASE_RE = re.compile(r"Release\s+(\d{4}_\d{2})")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--records-dir", required=True, type=Path)
    ap.add_argument("--stats-dir", required=True, type=Path)
    ap.add_argument("--library", required=True, type=Path)
    ap.add_argument("--library-csv", required=True)
    ap.add_argument("--output-dir", required=True, type=Path)
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = a.library / a.library_csv
    with open(csv_path) as fh:
        lib_rows = list(csv.DictReader(fh))
    readme = (a.library / "README").read_text(errors="replace") if (a.library / "README").exists() else ""
    m = RELEASE_RE.search(readme)
    db = a.output_dir / "seq_index.sqlite"
    if db.exists():
        db.unlink()
    con = sqlite3.connect(db)
    con.executescript("""
        create table seq(md5 text, accession text, proteome_id text, taxid integer, reviewed integer);
        create table acc(accession text primary key, proteome_id text, taxid integer);
        create table refseq(refseq_id text, accession text);
        create table proteome(proteome_id text primary key, taxid integer, species_name text,
                              file_prefix text);
    """)
    stats = {p.stem: json.loads(p.read_text()) for p in a.stats_dir.glob("*.json")}
    n_records = 0
    for row in lib_rows:
        pid = row["proteome_id"]
        if pid not in stats:
            continue
        con.execute("insert into proteome values (?,?,?,?)",
                    (pid, int(row["tax_id"]), row["species_name"], row["file_prefix"]))
        with open_maybe_compressed(a.records_dir / f"{pid}.records.tsv.zst") as fh:
            batch_seq, batch_acc, batch_ref = [], [], []
            for r in csv.DictReader(fh, delimiter="\t"):
                tax = int(r["taxon_id"] or row["tax_id"])
                batch_seq.append((r["seq_md5"], r["accession"], pid, tax, int(r["reviewed"])))
                batch_acc.append((r["accession"], pid, tax))
                for x in r["xrefs"].split("|"):
                    if x.startswith("RefSeq:"):
                        batch_ref.append((x.split(":", 1)[1], r["accession"]))
            con.executemany("insert into seq values (?,?,?,?,?)", batch_seq)
            con.executemany("insert or ignore into acc values (?,?,?)", batch_acc)
            con.executemany("insert into refseq values (?,?)", batch_ref)
            n_records += len(batch_seq)
        con.commit()
    con.executescript("create index seq_md5 on seq(md5); create index refseq_id on refseq(refseq_id);")
    con.commit()
    con.close()
    missing = [r["proteome_id"] for r in lib_rows if r["proteome_id"] not in stats]
    manifest = {
        "format_version": FORMAT_VERSION,
        "library": str(a.library),
        "release": m.group(1) if m else "",
        "library_csv": a.library_csv,
        "library_csv_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
        "n_proteomes": len(stats),
        "n_records": n_records,
        "missing": missing,
        "proteomes": [stats[k] for k in sorted(stats)],
        "build_date": datetime.date.today().isoformat(),
    }
    (a.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=1))
    if a.records_dir.resolve() != (a.output_dir / "records").resolve():
        (a.output_dir / "records").mkdir(exist_ok=True)
        for p in a.records_dir.glob("*.records.tsv.zst"):
            target = a.output_dir / "records" / p.name
            if not target.exists():
                target.write_bytes(p.read_bytes())
    print(f"index: {len(stats)} proteomes, {n_records} records, {len(missing)} missing",
          file=sys.stderr)


if __name__ == "__main__":
    main()
```

  Note: `manifest.json` is written last and the records are copied after it. `UniProtIndex` checks all three, so an interrupted build still fails the completeness check.

- [ ] **Step 5: Run the tests to verify they pass, then commit.**

```bash
chmod +x bin/uniprot_build_index.py
pixi run pytest tests/test_uniprot_index.py tests/test_uniprot_parse_build.py -v
git add lib/uniprot_index.py bin/uniprot_build_index.py tests/test_uniprot_index.py
git commit -m "UniProt library index: SQLite build + reader with format/completeness checks (#A)"
```

---

### Task 5: Per-proteome matching (`bin/uniprot_link.py`)

**Files:**
- Create: `bin/uniprot_link.py`
- Test: `tests/test_uniprot_link.py`

**Interfaces:**
- Consumes: `UniProtIndex` (Task 4), `seq_md5` (Task 1), and `build_fixture_index` from `tests/test_uniprot_index.py` (import it).
- Produces:
  - `choose(cands: list[dict], own_taxid: int | None, restrict: str | None) -> tuple[dict | None, str]`: returns the chosen candidate and the match type (`seq_own` / `seq_other`), or `(None, '')`;
  - `accession_from_id(pid: str) -> str | None`;
  - `LINK_COLUMNS` (exact order in Global Constraints);
  - CLI: `--index DIR --protein-fasta FA --short S --species NAME [--taxid N] [--restrict-proteome UPID] --output TSV`.

- [ ] **Step 1: Write the failing tests.**

```python
import csv
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "bin"))
sys.path.insert(0, str(REPO / "lib"))
sys.path.insert(0, str(REPO / "tests"))
import uniprot_link as ul  # noqa: E402
from test_uniprot_index import build_fixture_index  # noqa: E402


def test_accession_from_id():
    assert ul.accession_from_id("tr|Q7S6W2|Q7S6W2_NEUCR") == "Q7S6W2"
    assert ul.accession_from_id("sp|P12345-3|X_YEAST") == "P12345"
    assert ul.accession_from_id("Q7S6W2-2") == "Q7S6W2"
    assert ul.accession_from_id("NCU05603T0") is None
    assert ul.accession_from_id("XP_960565.1") is None


def test_choose_order():
    own = {"accession": "B2", "proteome_id": "UP1", "taxid": 1, "reviewed": 0}
    other_rev = {"accession": "A1", "proteome_id": "UP2", "taxid": 2, "reviewed": 1}
    other_unrev = {"accession": "A0", "proteome_id": "UP3", "taxid": 3, "reviewed": 0}
    assert ul.choose([other_rev, own], 1, None) == (own, "seq_own")
    assert ul.choose([other_unrev, other_rev], 1, None) == (other_rev, "seq_other")
    assert ul.choose([other_unrev, other_rev], None, None) == (other_rev, "seq_other")
    assert ul.choose([own, other_rev], 1, "UP2") == (other_rev, "seq_own")
    assert ul.choose([], 1, None) == (None, "")


def _run(tmp_path, fasta, *extra):
    idx = build_fixture_index(tmp_path)
    fa = tmp_path / "p.fa"
    fa.write_text(fasta)
    out = tmp_path / "out.tsv"
    subprocess.run([sys.executable, str(REPO / "bin" / "uniprot_link.py"), "--index", str(idx),
                    "--protein-fasta", str(fa), "--short", "Ncra",
                    "--species", "Neurospora crassa", *extra, "--output", str(out)],
                   check=True, capture_output=True, text=True)
    with open(out) as fh:
        return {r["protein_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def test_match_rules_end_to_end(tmp_path):
    rows = _run(tmp_path,
                ">tr|Q7S6W2|Q7S6W2_NEUCR\nMKVLLAQ\n"   # id
                ">XP_960565.1\nXXXX\n"                 # refseq
                ">NCU_gene1\nmkvllaq*\n"               # seq, normalized; own species wins
                ">NCU_gene2\nMSEQTWO\n"                # seq, own species (only candidate)
                ">NCU_gene3\nNOMATCH\n")               # none
    assert rows["tr|Q7S6W2|Q7S6W2_NEUCR"]["uniprot_match"] == "id"
    assert rows["XP_960565.1"]["uniprot_match"] == "refseq"
    assert rows["XP_960565.1"]["uniprot_accession"] == "Q7S6W2"
    assert rows["NCU_gene1"]["uniprot_match"] == "seq_own"
    assert rows["NCU_gene1"]["uniprot_accession"] == "Q7S6W2"
    assert rows["NCU_gene1"]["uniprot_n_matches"] == "2"
    assert rows["NCU_gene2"]["uniprot_accession"] == "P00001"
    assert "NCU_gene3" not in rows
    assert list(next(iter(rows.values())).keys()) == list(ul.LINK_COLUMNS)


def test_species_without_own_proteome_falls_to_seq_other(tmp_path):
    idx = build_fixture_index(tmp_path)
    fa = tmp_path / "p.fa"
    fa.write_text(">g1\nMKVLLAQ\n")
    out = tmp_path / "o.tsv"
    p = subprocess.run([sys.executable, str(REPO / "bin" / "uniprot_link.py"), "--index", str(idx),
                        "--protein-fasta", str(fa), "--short", "Xxx", "--species", "Other fungus",
                        "--output", str(out)], check=True, capture_output=True, text=True)
    row = next(csv.DictReader(open(out), delimiter="\t"))
    assert row["uniprot_match"] == "seq_other"
    assert row["uniprot_accession"] == "P00002"          # reviewed wins the tie
    assert "Saccharomyces cerevisiae" in row["uniprot_match_species"]
    assert "seq_other=1" in p.stderr and "none=0" in p.stderr


def test_unusable_index_is_a_hard_error(tmp_path):
    fa = tmp_path / "p.fa"
    fa.write_text(">g1\nMK\n")
    p = subprocess.run([sys.executable, str(REPO / "bin" / "uniprot_link.py"), "--index",
                        str(tmp_path / "nope"), "--protein-fasta", str(fa), "--short", "X",
                        "--species", "X y", "--output", str(tmp_path / "o.tsv")],
                       capture_output=True, text=True)
    assert p.returncode != 0 and "incomplete" in p.stderr
```

- [ ] **Step 2: Run them to verify they fail.** Run: `pixi run pytest tests/test_uniprot_link.py -v`. Expected: FAIL (module missing).

- [ ] **Step 3: Implement `bin/uniprot_link.py`.**

```python
#!/usr/bin/env python3
"""Match one proteome to UniProt records via a library index: protein ID as
UniProt accession, then RefSeq ID, then exact sequence (own species first, then
any). Writes <Short>.uniprot_link.tsv keyed by protein_id. See
docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md."""
import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import open_maybe_compressed  # noqa: E402
from uniprot_dat import seq_md5  # noqa: E402
from uniprot_index import IndexFormatError, UniProtIndex  # noqa: E402

ACC_RE = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
LINK_COLUMNS = (
    "protein_id", "uniprot_accession", "uniprot_gene_name", "uniprot_description",
    "uniprot_go_ids", "uniprot_pfam_ids", "uniprot_pfam_names", "uniprot_interpro_ids",
    "uniprot_ec_numbers", "uniprot_alphafold_id", "uniprot_xrefs", "uniprot_match",
    "uniprot_match_species", "uniprot_reviewed", "uniprot_pubs", "uniprot_n_matches",
)
REC_TO_COL = {"gene_name": "uniprot_gene_name", "description": "uniprot_description",
              "go_ids": "uniprot_go_ids", "pfam_ids": "uniprot_pfam_ids",
              "pfam_names": "uniprot_pfam_names", "interpro_ids": "uniprot_interpro_ids",
              "ec_numbers": "uniprot_ec_numbers", "alphafold_id": "uniprot_alphafold_id",
              "xrefs": "uniprot_xrefs", "reviewed": "uniprot_reviewed", "pubs": "uniprot_pubs"}


def accession_from_id(pid):
    m = re.match(r"^(?:sp|tr)\|([^|]+)\|", pid)
    token = m.group(1) if m else pid.split()[0]
    token = token.split("-")[0]
    return token if ACC_RE.match(token) else None


def choose(cands, own_taxid, restrict):
    if not cands:
        return None, ""
    def own(c):
        return c["proteome_id"] == restrict if restrict else (
            own_taxid is not None and c["taxid"] == own_taxid)
    ranked = sorted(cands, key=lambda c: (not own(c), -int(c["reviewed"]), c["accession"]))
    best = ranked[0]
    return best, ("seq_own" if own(best) else "seq_other")


def read_fasta(path):
    pid, chunks = None, []
    with open_maybe_compressed(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if pid is not None:
                    yield pid, "".join(chunks)
                pid, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if pid is not None:
        yield pid, "".join(chunks)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--index", required=True, type=Path)
    ap.add_argument("--protein-fasta", required=True, type=Path)
    ap.add_argument("--short", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--taxid", type=int, default=None)
    ap.add_argument("--restrict-proteome", default=None, dest="restrict")
    ap.add_argument("--output", required=True, type=Path)
    a = ap.parse_args()
    try:
        idx = UniProtIndex(a.index)
    except IndexFormatError as e:
        sys.exit(f"ERROR: {e}")
    own_taxid = a.taxid if a.taxid is not None else idx.taxid_for_species(a.species)

    chosen = {}          # pid -> (candidate dict with accession/proteome_id, match, n)
    counts = Counter()
    n_total = 0
    for pid, seq in read_fasta(a.protein_fasta):
        n_total += 1
        acc = accession_from_id(pid)
        hit = idx.by_accession(acc) if acc else None
        if hit:
            chosen[pid] = (hit, "id", 1)
        else:
            refs = idx.by_refseq(pid)
            if refs:
                hit = idx.by_accession(refs[0])
                if hit:
                    chosen[pid] = (hit, "refseq", len(refs))
            if pid not in chosen:
                cands = idx.by_md5(seq_md5(seq))
                best, kind = choose(cands, own_taxid, a.restrict)
                if best:
                    chosen[pid] = (best, kind, len(cands))
        counts[chosen[pid][1] if pid in chosen else "none"] += 1

    by_proteome = defaultdict(set)
    for hit, _, _ in chosen.values():
        by_proteome[hit["proteome_id"]].add(hit["accession"])
    recs = {}
    for prot, accs in by_proteome.items():
        for acc, r in idx.records(prot, accs).items():
            recs[(prot, acc)] = r

    with open(a.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=LINK_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for pid, (hit, kind, n) in chosen.items():
            r = recs.get((hit["proteome_id"], hit["accession"]))
            if r is None:
                continue
            row = {c: "" for c in LINK_COLUMNS}
            row.update(protein_id=pid, uniprot_accession=hit["accession"], uniprot_match=kind,
                       uniprot_match_species=idx.species_name(hit["proteome_id"]),
                       uniprot_n_matches=str(n))
            for k, col in REC_TO_COL.items():
                row[col] = r[k]
            w.writerow(row)
    summary = " ".join(f"{k}={counts[k]}" for k in ("id", "refseq", "seq_own", "seq_other", "none"))
    print(f"{a.short}: {n_total} proteins; {summary} (own taxid {own_taxid})", file=sys.stderr)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the tests to verify they pass, then commit.**

```bash
chmod +x bin/uniprot_link.py
pixi run pytest tests/test_uniprot_link.py -v
git add bin/uniprot_link.py tests/test_uniprot_link.py
git commit -m "uniprot_link: match a proteome to UniProt by id, RefSeq, or exact sequence (#B)"
```

---

### Task 6: Annotation columns (`bin/annotate_presence_matrix.py`)

**Files:**
- Modify: `bin/annotate_presence_matrix.py:87-90` (`UNIPROT_XREF_COLS`) and the `load_uniprot_xrefs` docstring
- Test: `tests/test_annotate_presence_matrix.py`

**Interfaces:**
- Consumes: `LINK_COLUMNS` (Task 5).
- Produces: `presence_matrix.function.tsv` with the 15 `uniprot_*` columns.

- [ ] **Step 1: Write the failing test.** Add to `tests/test_annotate_presence_matrix.py`. Follow the file's existing pattern of running the script on a tiny matrix, and add a `--uniprot_xref_files` input TSV with header `LINK_COLUMNS` and one row for `n1` with `uniprot_match=seq_other`, `uniprot_pubs=1;;protein;T`:

```python
def test_uniprot_link_columns_pass_through(tmp_path):
    matrix = tmp_path / "m.tsv"
    matrix.write_text("protein_id\tsource_proteome\tA\nn1\tA\t1\nn2\tA\t1\n")
    link = tmp_path / "A.uniprot_link.tsv"
    cols = ["protein_id", "uniprot_accession", "uniprot_gene_name", "uniprot_description",
            "uniprot_go_ids", "uniprot_pfam_ids", "uniprot_pfam_names", "uniprot_interpro_ids",
            "uniprot_ec_numbers", "uniprot_alphafold_id", "uniprot_xrefs", "uniprot_match",
            "uniprot_match_species", "uniprot_reviewed", "uniprot_pubs", "uniprot_n_matches"]
    vals = ["n1", "Q7S6W2", "", "", "", "", "", "", "", "Q7S6W2", "", "seq_other",
            "Saccharomyces cerevisiae", "0", "1;;protein;T", "2"]
    link.write_text("\t".join(cols) + "\n" + "\t".join(vals) + "\n")
    out = tmp_path / "out.tsv"
    subprocess.run([sys.executable, str(REPO / "bin" / "annotate_presence_matrix.py"),
                    "--matrix", str(matrix), "--output", str(out),
                    "--uniprot_xref_files", str(link)], check=True)
    rows = {r["protein_id"]: r for r in csv.DictReader(open(out), delimiter="\t")}
    assert rows["n1"]["uniprot_match"] == "seq_other"
    assert rows["n1"]["uniprot_pubs"] == "1;;protein;T"
    assert rows["n2"]["uniprot_match"] == ""
```

  (The file already imports `csv`, `subprocess`, `sys` and defines `REPO`. The script requires only `--matrix` and `--output`, checked at `bin/annotate_presence_matrix.py:113,131`.)

- [ ] **Step 2: Run it to verify it fails.** Expected: FAIL (`KeyError: 'uniprot_match'`).

- [ ] **Step 3: Implement.**

```python
UNIPROT_XREF_COLS = ['uniprot_accession', 'uniprot_gene_name', 'uniprot_description',
                      'uniprot_go_ids', 'uniprot_pfam_ids', 'uniprot_pfam_names',
                      'uniprot_interpro_ids', 'uniprot_ec_numbers', 'uniprot_alphafold_id',
                      'uniprot_xrefs', 'uniprot_match', 'uniprot_match_species',
                      'uniprot_reviewed', 'uniprot_pubs', 'uniprot_n_matches']
```

  Update the `load_uniprot_xrefs` docstring: the inputs are now `bin/uniprot_link.py` TSVs, and protein IDs are unique per run because each protein comes from one proteome FASTA.

- [ ] **Step 4: Run the tests to verify they pass, then commit.**

```bash
pixi run pytest tests/test_annotate_presence_matrix.py -v
git add bin/annotate_presence_matrix.py tests/test_annotate_presence_matrix.py
git commit -m "annotate_presence_matrix: carry uniprot_match/species/reviewed/pubs/n_matches (#B)"
```

---

### Task 7: Index workflow (`-entry UNIPROT_INDEX`)

**Files:**
- Create: `modules/uniprot_index.nf`, `workflows/uniprot_index.nf`
- Modify: `main.nf` (add a named workflow), `nextflow.config` (params)
- Test: a `-profile test` stub run against the fixture library

**Interfaces:**
- Consumes: `bin/uniprot_plan_chunks.py`, `bin/uniprot_parse_dat.py`, `bin/uniprot_build_index.py`.
- Produces: the index directory at `params.uniprot_index`: `manifest.json`, `seq_index.sqlite`, `records/`.

Note: the three processes live in one module file because they are only ever used together. They pass one contract (the index directory). This follows the CLAUDE.md rule "one process per file" in spirit for tools, while keeping the index build readable. If a reviewer prefers the strict rule, split the file into `modules/uniprot_plan_chunks.nf`, `modules/uniprot_parse_chunk.nf` and `modules/uniprot_seq_index.nf` with no other change.

- [ ] **Step 1: Params in `nextflow.config`** (inside `params { }`, next to `pfam_hmm`):

```groovy
    uniprot_index        = null   // UniProt library index dir (-entry UNIPROT_INDEX output); enables UNIPROT_LINK
    uniprot_library      = null   // pre-downloaded UniProt library (dir with data/<prefix>.dat.gz), for -entry UNIPROT_INDEX
    uniprot_library_csv  = null   // proteome CSV file name inside uniprot_library (proteome_id,tax_id,species_name,file_prefix,...)
    uniprot_index_chunk_gb = 8    // compressed .dat.gz per UNIPROT_PARSE_CHUNK job (~1-1.5 h at 6.4 MB/s)
```

- [ ] **Step 2: `modules/uniprot_index.nf`.**

```groovy
nextflow.enable.dsl=2

// UniProt library index build (docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md).
process UNIPROT_PLAN_CHUNKS {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    input:
    val(library)
    output:
    path("chunks.tsv"), emit: chunks
    script:
    """
    uniprot_plan_chunks.py --library ${library} --library-csv ${params.uniprot_library_csv} \\
        --chunk-gb ${params.uniprot_index_chunk_gb} --output chunks.tsv
    """
}

process UNIPROT_PARSE_CHUNK {
    label 'low_cpu'
    tag "${chunk_id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    input:
    tuple val(chunk_id), path(chunks_tsv), val(library)
    output:
    path("records/*.records.tsv.zst"), emit: records
    path("stats/*.json"),              emit: stats
    script:
    """
    uniprot_parse_dat.py --library ${library} --chunks ${chunks_tsv} \\
        --chunk-id ${chunk_id} --outdir .
    """
}

process UNIPROT_SEQ_INDEX {
    label 'low_cpu'
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    storeDir { params.uniprot_index }
    input:
    path(records, stageAs: 'in_records/*')
    path(stats,   stageAs: 'in_stats/*')
    val(library)
    output:
    path("manifest.json")
    path("seq_index.sqlite")
    path("records")
    script:
    """
    mkdir -p records
    cp in_records/* records/
    uniprot_build_index.py --records-dir records --stats-dir in_stats \\
        --library ${library} --library-csv ${params.uniprot_library_csv} --output-dir .
    """
}
```

- [ ] **Step 3: `workflows/uniprot_index.nf`.**

```groovy
nextflow.enable.dsl=2

include { UNIPROT_PLAN_CHUNKS; UNIPROT_PARSE_CHUNK; UNIPROT_SEQ_INDEX } from '../modules/uniprot_index'

workflow UNIPROT_INDEX_BUILD {
    take:
    library   // val: absolute path to the UniProt library directory

    main:
    UNIPROT_PLAN_CHUNKS(library)
    chunk_ids = UNIPROT_PLAN_CHUNKS.out.chunks
        .splitCsv(header: true, sep: '\t')
        .map { it.chunk_id }
        .unique()
    UNIPROT_PARSE_CHUNK(chunk_ids.combine(UNIPROT_PLAN_CHUNKS.out.chunks).combine(library))
    UNIPROT_SEQ_INDEX(UNIPROT_PARSE_CHUNK.out.records.flatten().collect(),
                      UNIPROT_PARSE_CHUNK.out.stats.flatten().collect(), library)

    emit:
    manifest = UNIPROT_SEQ_INDEX.out[0]
}
```

- [ ] **Step 4: `main.nf`.** Add the include and a named workflow at the end of the file:

```groovy
include { UNIPROT_INDEX_BUILD } from './workflows/uniprot_index'

// Run with: nextflow run main.nf -entry UNIPROT_INDEX --uniprot_library <dir>
//   --uniprot_library_csv <csv> --uniprot_index <out dir>
workflow UNIPROT_INDEX {
    if (!params.uniprot_library || !params.uniprot_library_csv || !params.uniprot_index) {
        error "UNIPROT_INDEX needs --uniprot_library, --uniprot_library_csv and --uniprot_index"
    }
    UNIPROT_INDEX_BUILD(Channel.value(file(params.uniprot_library).toAbsolutePath().toString()))
}
```

- [ ] **Step 5: Verify on the fixture library.** Run it outside the repo, per CLAUDE.md:

```bash
OUT=$SCRATCH/uniprot_index_test
nextflow run main.nf -entry UNIPROT_INDEX -profile test \
    --uniprot_library $PWD/tests/data/uniprot_library --uniprot_library_csv proteomes.csv \
    --uniprot_index $OUT/index --outdir $OUT/results -work-dir $OUT/work
pixi run python -c "import sys; sys.path.insert(0,'lib'); from uniprot_index import UniProtIndex; \
i=UniProtIndex('$OUT/index'); print(i.manifest['n_records'], i.manifest['missing'])"
```

Expected: `3 ['UP000000003']`. A second identical run completes with `UNIPROT_SEQ_INDEX` skipped (storeDir).

  If `-profile test` requires `--config`/`--data_dir` for the default workflow's param checks, pass `--config tests/data/test.csv --data_dir tests/data`. The named entry does not use them.

- [ ] **Step 6: Commit.**

```bash
git add modules/uniprot_index.nf workflows/uniprot_index.nf main.nf nextflow.config
git commit -m "-entry UNIPROT_INDEX: build a UniProt library index (plan, parse, SQLite) (#A)"
```

---

### Task 8: Per-run `UNIPROT_LINK` wiring (replaces `UNIPROT_XREF`)

**Files:**
- Create: `modules/uniprot_link.nf`
- Modify: `main.nf:17` (include), `main.nf:99-116` (channel), `main.nf:141-147` (call)
- Delete: `modules/uniprot_xref.nf`, `bin/build_uniprot_refseq_xref.py`, `tests/test_build_uniprot_refseq_xref.py` (if present; check with `ls tests | grep refseq_xref`)
- Test: `-profile test` run with `--uniprot_index` pointing at the Task 7 fixture index

**Interfaces:**
- Consumes: `bin/uniprot_link.py` (Task 5); `ANNOTATE`'s existing `uniprot_xref_files` input (unchanged name).
- Produces: `uniprot_xref_files` = collected `*.uniprot_link.tsv`.

- [ ] **Step 1: `modules/uniprot_link.nf`.**

```groovy
nextflow.enable.dsl=2

// Match one proteome to UniProt records via the library index (bin/uniprot_link.py).
// Replaces UNIPROT_XREF (issue #92). Runs for every config proteome when --uniprot_index is set.
process UNIPROT_LINK {
    label 'low_cpu'
    tag "${meta.id}"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    tuple val(meta), path(protein_fa)
    val(index_dir)

    output:
    path("${meta.id}.uniprot_link.tsv"), emit: tsv

    script:
    def taxid    = meta.taxid ? "--taxid ${meta.taxid}" : ''
    def restrict = meta.uniprot_restrict ? "--restrict-proteome ${meta.uniprot_restrict}" : ''
    """
    uniprot_link.py --index ${index_dir} --protein-fasta ${protein_fa} \\
        --short ${meta.id} --species "${meta.species}" ${taxid} ${restrict} \\
        --output ${meta.id}.uniprot_link.tsv
    """
}
```

- [ ] **Step 2: `main.nf`.** Replace the `uniprot_xref_ch` block (lines 99–116) and the `UNIPROT_XREF` call (141–147):

```groovy
include { UNIPROT_LINK } from './modules/uniprot_link'
```

```groovy
    // UniProt linking (docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md):
    // every config proteome, when --uniprot_index is set. UniProtDatGz is a deprecated
    // alias (one release): its basename's UP id restricts the own-species proteome.
    uniprot_link_ch = Channel
        .fromPath(params.config)
        .splitCsv(header: true)
        .map { row ->
            def dat = row.UniProtDatGz?.trim()
            if (dat) log.warn "UniProtDatGz is deprecated (${row.Short}); use --uniprot_index. " +
                              "Treated as --restrict-proteome ${dat.tokenize('_')[0]}"
            [ [id: row.Short, species: row.Species, taxid: row.NCBI_TaxID?.trim() ?: null,
               uniprot_restrict: dat ? dat.tokenize('_')[0] : null],
              resolve_fa(row.Protein, ['pep', 'proteins']) ]
        }
```

```groovy
    if (params.uniprot_index) {
        UNIPROT_LINK(uniprot_link_ch, file(params.uniprot_index).toAbsolutePath().toString())
        uniprot_xref_files = UNIPROT_LINK.out.tsv.collect().ifEmpty([])
    } else {
        log.info "--uniprot_index not set: UniProt linking skipped"
        uniprot_xref_files = Channel.value([])
    }
```

  Leave the `ANNOTATE(..., uniprot_xref_files)` and `LOSS_ANNOTATE(..., uniprot_xref_files)` call sites unchanged.

- [ ] **Step 3: Delete the replaced files.**

```bash
git rm modules/uniprot_xref.nf bin/build_uniprot_refseq_xref.py
ls tests | grep -i refseq_xref && git rm tests/test_build_uniprot_refseq_xref.py
grep -rn "UNIPROT_XREF\|build_uniprot_refseq_xref" --include=*.nf --include=*.py --include=*.md . | grep -v docs/superpowers
```

  Expected: the final grep prints only CLAUDE.md/README lines, which Task 11 updates.

- [ ] **Step 4: Verify with `-profile test`.** The test config's proteins must be matchable by the fixture index. Add to `tests/data/` a protein whose sequence is `MKVLLAQ` in one test proteome's FASTA, or run `uniprot_link.py` output inspection only.

```bash
OUT=$SCRATCH/uniprot_link_test
nextflow run main.nf -profile test --uniprot_index $SCRATCH/uniprot_index_test/index \
    --outdir $OUT/results -work-dir $OUT/work
head -3 $OUT/results/*/presence_matrix.function.tsv | cut -f1,2 ; \
grep -c "uniprot_match" <(head -1 $OUT/results/*/presence_matrix.function.tsv)
nextflow run main.nf -profile test --outdir $OUT/results_noidx -work-dir $OUT/work2   # still succeeds, logs "skipped"
```

  Expected: the header contains `uniprot_match`, and the run without the index succeeds.

- [ ] **Step 5: Commit.**

```bash
git add modules/uniprot_link.nf main.nf
git commit -m "UNIPROT_LINK replaces UNIPROT_XREF; UniProtDatGz deprecated to a restrict alias (#B)"
```

---

### Task 9: Payload fields (`lib/report_data.py`)

**Files:**
- Modify: `lib/report_data.py`: `ROW_FIELDS` (~line 24), `CORE_ROW_FIELDS` (~643), `LOSSES_ROW_FIELDS` (~783), the three row builders (~580–606, ~740–765, ~955–980) and the three payload dicts (`'go_sets'` neighbours)
- Test: `tests/test_report_data.py`

**Interfaces:**
- Produces:
  - Row fields appended at the end of all three field lists: `'uacc'` (uniprot_accession), `'umatch'`, `'usp'` (uniprot_match_species, only when umatch == 'seq_other', else ''), `'urev'` (int 0/1, or -1 when no match), `'pubs'` (an index into `payload['pub_sets']`, or -1);
  - `payload['pub_sets']`: list of the unique `uniprot_pubs` strings (`_StringTable`, like `go_sets`).

- [ ] **Step 1: Write the failing test.** Add to `tests/test_report_data.py`. Use the file's existing helper that builds a payload from an in-memory matrix. Mirror whichever test already asserts `payload['go_sets']`, and add the uniprot columns to its matrix row:

```python
def test_uniprot_link_fields_reach_all_three_payloads(tmp_path):
    # matrix row n1 carries uniprot_accession=Q7S6W2, uniprot_match=seq_other,
    # uniprot_match_species="Saccharomyces cerevisiae", uniprot_reviewed=1,
    # uniprot_pubs="1;;protein;T"; row n2 has no uniprot columns set.
    for payload in _build_all_three(tmp_path):          # helper: novelty, core, losses payloads
        F = {n: i for i, n in enumerate(payload['fields'])}
        rows = {r[F['id']]: r for r in payload['rows']}
        n1 = rows['n1']
        assert n1[F['uacc']] == 'Q7S6W2'
        assert n1[F['umatch']] == 'seq_other'
        assert n1[F['usp']] == 'Saccharomyces cerevisiae'
        assert n1[F['urev']] == 1
        assert payload['pub_sets'][n1[F['pubs']]] == '1;;protein;T'
        n2 = rows['n2']
        assert n2[F['umatch']] == '' and n2[F['urev']] == -1 and n2[F['pubs']] == -1
```

  Write `_build_all_three(tmp_path)` in the test file. It writes a matrix/config pair that is valid for all three builders (reuse the fixtures already in the file for novelty and core; the losses payload needs `--query-group OUT`-shaped rows, see the existing losses tests). It returns `[build_payload(...), build_core_payload(...), build_losses_payload(...)]`. Copy the exact call signatures from the existing tests in the same file.

- [ ] **Step 2: Run it to verify it fails.** Run: `pixi run pytest tests/test_report_data.py -v -k uniprot_link`. Expected: FAIL (`KeyError: 'uacc'`).

- [ ] **Step 3: Implement.** In each of the three field lists, append:

```python
    'uacc',      # uniprot_accession (UNIPROT_LINK), or ''
    'umatch',    # uniprot_match: id | refseq | seq_own | seq_other | ''
    'usp',       # uniprot_match_species when umatch == 'seq_other', else ''
    'urev',      # uniprot_reviewed as int (1 Swiss-Prot, 0 TrEMBL), -1 when no match
    'pubs',      # index into payload['pub_sets'], or -1
```

  In each builder, create `pub_sets = _StringTable()` next to `go_sets`, and append to each row list, in the same order:

```python
            row.get('uniprot_accession', '') or '',
            row.get('uniprot_match', '') or '',
            (row.get('uniprot_match_species', '') or '') if row.get('uniprot_match') == 'seq_other' else '',
            int(row['uniprot_reviewed']) if (row.get('uniprot_reviewed') or '') != '' else -1,
            pub_sets.intern(row.get('uniprot_pubs', '') or ''),
```

  `_StringTable.intern('')` returns -1 (checked, `lib/report_data.py:111-113`), so an empty `uniprot_pubs` gives `pubs == -1` with no special case. Add `'pub_sets': pub_sets.table,` to each payload dict.

- [ ] **Step 4: Run the tests to verify they pass, then commit.**

```bash
pixi run pytest tests/test_report_data.py -v
git add lib/report_data.py tests/test_report_data.py
git commit -m "report_data: UniProt link fields (uacc/umatch/usp/urev/pubs) in all three payloads (#C)"
```

---

### Task 10: Card links and publications (`lib/report_common.py` + templates)

**Files:**
- Modify: `lib/report_common.py`: `XREF_LINK_TEMPLATES` (~line 700–775), `xrefLinkNodes` (~785), `externalLinksNode` (~801); add `publicationsNode`
- Modify: `lib/report_template.py` (~1124), `lib/core_report_template.py` (~451), `lib/losses_report_template.py` (~532): pass the new fields and append a Publications field
- Test: `tests/test_report_js_behaviour.py` (fixture matrix), `tests/js/drive_reports.mjs` (checks)

**Interfaces:**
- Consumes: the payload fields from Task 9.
- Produces: in JS, `externalLinksNode(o)` accepts the extra keys `uacc`, `af`, `umatch`, `usp`, and `publicationsNode(pubsStr) -> Node | null`.

- [ ] **Step 1: Write the failing jsdom checks.**
  In `tests/test_report_js_behaviour.py`, add a candidate row `n4` to `MATRIX` and `>n4\nMKVLLA...` to `_fasta()`, present in all three ingroup species like `n2`. Its values:
  - `Best_Swissprot` empty;
  - `uniprot_accession=Q7S6W2`, `uniprot_alphafold_id=Q7S6W2`, `uniprot_match=seq_other`;
  - `uniprot_match_species=Saccharomyces cerevisiae (strain S288c)`;
  - `uniprot_xrefs=VEuPathDB:FungiDB:NCU05603|PANTHER:PTHR10000|PDB:6YWS|Evil:<img src=x onerror=alert(1)>`;
  - `uniprot_pubs=12712197;10.1038/nature01554;proteome;Genome paper|99999999;;protein;<script>alert(1)</script> title`.

  Add the new columns to the `MATRIX` header, with empty cells on the existing rows. In `tests/js/drive_reports.mjs`, after the `shared` checks, add:

```js
  pick(allRows, 'n4').dispatchEvent(ev(w, 'click'));
  h = hrefs(detail());
  check('n4: own UniProt accession links to UniProt with no SwissProt hit',
        h.some((x) => x.includes('uniprot.org/uniprotkb/Q7S6W2')), h.join(' '));
  check('n4: own accession gives an AlphaFold link',
        h.some((x) => x.includes('alphafold.ebi.ac.uk/entry/Q7S6W2')));
  check('n4: seq_other suppresses gene-database links',
        !h.some((x) => x.includes('fungidb.org')), h.join(' '));
  check('n4: seq_other keeps family/structure links',
        h.some((x) => x.includes('pantherdb.org')) && h.some((x) => x.includes('rcsb.org/structure/6YWS')),
        h.join(' '));
  check('n4: seq_other match is labelled with the other species',
        /Identical sequence in Saccharomyces cerevisiae/.test(detail().textContent));
  check('n4: unknown xref DB renders no link and no markup',
        !detail().querySelector('img') && !h.some((x) => x.includes('onerror')));
  check('n4: protein-scope publication links to PubMed',
        h.some((x) => x.includes('pubmed.ncbi.nlm.nih.gov/99999999')), h.join(' '));
  check('n4: genome paper collapsed to one line',
        /Genome paper/.test(detail().textContent));
  check('n4: publication title is text, not markup',
        !detail().querySelector('script') && /<script>alert\(1\)<\/script> title/.test(detail().textContent));
```

- [ ] **Step 2: Run to verify it fails.** Run: `pixi run pytest tests/test_report_js_behaviour.py -v`. Expected: FAIL lines for n4. If jsdom is missing, set `NOVINVENIO_JSDOM` per CLAUDE.md "Testing template JS changes" first.

- [ ] **Step 3: Implement in `lib/report_common.py`.**
  Add templates to `XREF_LINK_TEMPLATES`, in the same shape as `GeneID`:

```js
    PDB:      { render: function (id) { return { label: "PDB " + id, href: "https://www.rcsb.org/structure/" + encodeURIComponent(id), title: "PDB structure " + id }; } },
    PANTHER:  { render: function (id) { return { label: "PANTHER " + id, href: "https://www.pantherdb.org/panther/family.do?clsAccession=" + encodeURIComponent(id), title: "PANTHER family " + id }; } },
    OrthoDB:  { render: function (id) { return { label: "OrthoDB " + id, href: "https://www.orthodb.org/?query=" + encodeURIComponent(id), title: "OrthoDB " + id }; } },
    STRING:   { render: function (id) { return { label: "STRING " + id, href: "https://string-db.org/network/" + encodeURIComponent(id), title: "STRING network " + id }; } },
    eggNOG:   { render: function (id) { return { label: "eggNOG " + id, href: "http://eggnog5.embl.de/#/app/results?target_nogs=" + encodeURIComponent(id), title: "eggNOG " + id }; } },
    Gene3D:   { render: function (id) { return { label: "Gene3D " + id, href: "https://www.cathdb.info/version/latest/superfamily/" + encodeURIComponent(id.replace(/^G3DSA:/, "")), title: "CATH/Gene3D " + id }; } },
    SUPFAM:   { render: function (id) { return { label: "SUPFAM " + id, href: "https://supfam.org/SUPERFAMILY/cgi-bin/scop.cgi?ipid=" + encodeURIComponent(id), title: "SUPERFAMILY " + id }; } },
    PROSITE:  { render: function (id) { return { label: "PROSITE " + id, href: "https://prosite.expasy.org/" + encodeURIComponent(id), title: "PROSITE " + id }; } },
    SMART:    { render: function (id) { return { label: "SMART " + id, href: "https://smart.embl.de/smart/do_annotation.pl?DOMAIN=" + encodeURIComponent(id), title: "SMART " + id }; } },
    CDD:      { render: function (id) { return { label: "CDD " + id, href: "https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=" + encodeURIComponent(id), title: "NCBI CDD " + id }; } },
    PRINTS:   { render: function (id) { return { label: "PRINTS " + id, href: "https://www.ebi.ac.uk/interpro/entry/prints/" + encodeURIComponent(id), title: "PRINTS " + id }; } },
    PIRSF:    { render: function (id) { return { label: "PIRSF " + id, href: "https://www.ebi.ac.uk/interpro/entry/pirsf/" + encodeURIComponent(id), title: "PIRSF " + id }; } },
    HAMAP:    { render: function (id) { return { label: "HAMAP " + id, href: "https://hamap.expasy.org/rule/" + encodeURIComponent(id), title: "HAMAP rule " + id }; } },
    NCBIfam:  { render: function (id) { return { label: "NCBIfam " + id, href: "https://www.ebi.ac.uk/interpro/entry/ncbifam/" + encodeURIComponent(id), title: "NCBIfam " + id }; } },
    FunFam:   { render: function (id) { return { label: "FunFam " + id, href: "https://www.cathdb.info/search?q=" + encodeURIComponent(id), title: "CATH FunFam " + id }; } },
    MEROPS:   { render: function (id) { return { label: "MEROPS " + id, href: "https://www.ebi.ac.uk/merops/cgi-bin/pepsum?id=" + encodeURIComponent(id), title: "MEROPS " + id }; } },
    CAZy:     { render: function (id) { return { label: "CAZy " + id, href: "http://www.cazy.org/" + encodeURIComponent(id) + ".html", title: "CAZy family " + id }; } },
    ESTHER:   { render: function (id) { return { label: "ESTHER " + id, href: "https://bioweb.supagro.inrae.fr/ESTHER/gene_locus?name=" + encodeURIComponent(id), title: "ESTHER " + id }; } },
    TCDB:     { render: function (id) { return { label: "TCDB " + id, href: "https://www.tcdb.org/search/result.php?tc=" + encodeURIComponent(id), title: "TCDB " + id }; } },
    BRENDA:   { render: function (id) { return { label: "BRENDA " + id, href: "https://www.brenda-enzymes.org/enzyme.php?ecno=" + encodeURIComponent(id), title: "BRENDA EC " + id }; } },
    UniPathway: { render: function (id) { return { label: "UniPathway " + id, href: "https://www.uniprot.org/unipathway/" + encodeURIComponent(id), title: "UniPathway " + id }; } }
```

  Gene-database set and suppression, in `xrefLinkNodes(xrefsStr, suppressGeneDb)`:

```js
  var GENE_DB_XREFS = { VEuPathDB: 1, GeneID: 1, KEGG: 1, RefSeq: 1, EnsemblFungi: 1, EnsemblBacteria: 1 };
  function xrefLinkNodes(xrefsStr, suppressGeneDb) {
    var links = [];
    var hasVEuPathDB = false;
    (xrefsStr ? xrefsStr.split("|") : []).filter(Boolean).forEach(function (entry) {
      var i = entry.indexOf(":");
      if (i === -1) return;
      var db = entry.slice(0, i);
      var id = entry.slice(i + 1);
      var tmpl = XREF_LINK_TEMPLATES[db];
      if (!tmpl || !id) return;
      if (suppressGeneDb && GENE_DB_XREFS[db]) return;
      var r = tmpl.render(id);
      links.push(extLink(r.label, r.href, r.title));
      if (db === "VEuPathDB") hasVEuPathDB = true;
    });
    return { links: links, hasVEuPathDB: hasVEuPathDB };
  }
```

  In `externalLinksNode(o)`, replace the `var acc = uniprotAcc(o.sprot); if (acc) {...}` block and the xref call:

```js
    var own = o.uacc || "";
    var hitAcc = uniprotAcc(o.sprot);
    var acc = own || hitAcc;
    if (acc) {
      links.appendChild(extLink("UniProt " + acc,
        "https://www.uniprot.org/uniprotkb/" + encodeURIComponent(acc) + "/entry"));
      links.appendChild(extLink("AlphaFold",
        "https://alphafold.ebi.ac.uk/entry/" + encodeURIComponent(o.af || acc),
        "Predicted structure for " + (o.af || acc)));
    }
    if (own && hitAcc && hitAcc !== own) {
      links.appendChild(extLink("Best SwissProt hit: " + hitAcc,
        "https://www.uniprot.org/uniprotkb/" + encodeURIComponent(hitAcc) + "/entry"));
    }
    var otherSpecies = o.umatch === "seq_other";
    var db = otherSpecies ? null : genomeDbLink(o.proteome && o.proteome.source_db, o.id);
    var xrefResult = xrefLinkNodes(o.xrefs, otherSpecies);
```

  After `box.appendChild(links);`, add:

```js
    if (otherSpecies && o.usp) {
      box.appendChild(el("p", "links-note",
        "Identical sequence in " + o.usp + " (" + own + "); gene-database links omitted."));
    }
```

  Add `publicationsNode`, next to `alphafoldLinkNode`:

```js
  // uniprot_pubs: "PMID;DOI;scope;title|..." (bin/uniprot_link.py). Genome-sequencing
  // papers (scope=proteome) collapse to one line; up to 5 protein-scope papers listed.
  function publicationsNode(pubsStr) {
    if (!pubsStr) return null;
    var box = el("div", "pubs");
    var genome = null, n = 0;
    pubsStr.split("|").forEach(function (entry) {
      var p = entry.split(";");
      var pmid = p[0] || "", doi = p[1] || "", scope = p[2] || "", title = p.slice(3).join(";");
      var href = pmid ? "https://pubmed.ncbi.nlm.nih.gov/" + encodeURIComponent(pmid) + "/"
                      : (doi ? "https://doi.org/" + encodeURI(doi) : "");
      if (!href) return;
      if (scope === "proteome") { if (!genome) genome = { href: href, title: title }; return; }
      if (n >= 5) return;
      n += 1;
      var row = el("div", "pub");
      row.appendChild(extLink(title || (pmid ? "PubMed " + pmid : doi), href,
                              pmid ? "PubMed " + pmid : "DOI " + doi));
      box.appendChild(row);
    });
    if (genome) {
      var g = el("div", "pub");
      g.appendChild(extLink("Genome paper: " + (genome.title || ""), genome.href, "Proteome-wide reference"));
      box.appendChild(g);
    }
    return box.childNodes.length ? box : null;
  }
```

  `extLink(label, href, title)` builds its label through `el("a", null, label)` (`lib/report_common.py:355`), which sets `textContent`, so titles render as text.

- [ ] **Step 4: Templates.** In each of the three `externalLinksNode({...})` calls, add:

```js
      uacc: row[F.uacc],
      af: row[F.af],
      umatch: row[F.umatch],
      usp: row[F.usp],
```

  After the External resources field in each template, add:

```js
    var pubsNode = row[F.pubs] >= 0 ? publicationsNode(DATA.pub_sets[row[F.pubs]]) : null;
    if (pubsNode) detailEl.appendChild(field("Publications", pubsNode));
```

  For `losses_report_template.py`, keep its "(outgroup protein)" label pattern: `field("Publications (outgroup protein)", pubsNode)`.

- [ ] **Step 5: Run the tests to verify they pass.**

```bash
pixi run pytest tests/test_report_js_behaviour.py tests/test_report_templates.py tests/test_skins.py -v
```

  Expected: PASS. This includes `test_page_javascript_parses` and the no-hardcoded-colour test. The CSS for `.pubs`/`.pub` must use existing tokens only.

- [ ] **Step 6: Commit.**

```bash
git add lib/report_common.py lib/report_template.py lib/core_report_template.py \
        lib/losses_report_template.py tests/test_report_js_behaviour.py tests/js/drive_reports.mjs
git commit -m "Reports: UniProt/AlphaFold from own accession, 21 more xref DBs, publications, seq_other rule (#C)"
```

---

### Task 11: Documentation

**Files:**
- Modify: `CLAUDE.md` (Repository Layout, ANNOTATE section item 5, Key Parameters table, Analysis Config CSV Format `UniProtDatGz` paragraph, Modifying `lib/`), `README.md` (UniProt section)
- Modify: `.living/decisions.md` (append), then regenerate `.living/INDEX.md`

- [ ] **Step 1: CLAUDE.md edits.**
  - Layout: replace `uniprot_xref.nf` with `uniprot_link.nf` and `uniprot_index.nf`. Add `bin/uniprot_plan_chunks.py`, `bin/uniprot_parse_dat.py`, `bin/uniprot_build_index.py`, `bin/uniprot_link.py` and `lib/uniprot_index.py`, one line each, describing what they do. Remove `build_uniprot_refseq_xref.py`.
  - ANNOTATE item 5: rewrite the "UniProt cross-reference lookup (issue #92)" paragraph. Describe `UNIPROT_LINK`, the four match types, the own-species-first rule and the index requirement, and link the spec.
  - Key Parameters: add rows for `--uniprot_index`, `--uniprot_library`, `--uniprot_library_csv` and `--uniprot_index_chunk_gb`, with the defaults from Task 7.
  - Config CSV: replace the `UniProtDatGz` paragraph with a deprecation note: it restricts the own-species match to the proteome named by its `UP…` prefix, and it will be removed in the next release.
  - Add a "Building the UniProt index" command block with the Task 12 Step 1 command.

- [ ] **Step 2: `.living/decisions.md`.** Append:

```markdown
## 2026-09-23 — UniProt linking moves into the pipeline (library index + UNIPROT_LINK)

**Context**: UniProt annotation reached reports only through NII's post-run merge (UniProt-sourced studies) or UNIPROT_XREF (RefSeq-only, hand-set UniProtDatGz). Q7S6W2 lost its links when a study dir lacked annotations/.
**Decision**: One-time -entry UNIPROT_INDEX over a pre-downloaded library; per-run UNIPROT_LINK matching by id, RefSeq, exact sequence (own species first, then any).
**Alternatives**: parse each species' .dat per run (no cross-species fallback); UniProt REST per protein (network, rate limits, not reproducible).
**Rationale**: spec docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md.
```

```bash
/usr/bin/python3.12 /rhome/jstajich/.claude/plugins/marketplaces/mycelium/skills/core/scripts/generate_index.py --living-dir .living
```

- [ ] **Step 3: Commit.**

```bash
git add CLAUDE.md README.md .living/decisions.md .living/INDEX.md
git commit -m "Docs: UniProt library index + UNIPROT_LINK; UniProtDatGz deprecated (#C)"
```

---

### Task 12: Build the real index and accept on sordariales_shallow

**Files:** none in the repo. Outputs go to the shared library and to `$SCRATCH`/NII. Findings are recorded in the PR body and `.living/learnings.md`.

- [ ] **Step 1: Build the `Fungi_2026_03` index** (SLURM head job from the worktree; workers on preempt as in NII's `conf/diamond_sensitivity_bench.config`):

```bash
nextflow run /bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/uniprot-library-index/main.nf \
    -entry UNIPROT_INDEX -profile slurm \
    -c /bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/uniprot-library-index/conf/ucr_hpcc_slurm.config \
    --uniprot_library /bigdata/stajichlab/shared/db/Uniprot/Fungi_2026_03 \
    --uniprot_library_csv fungi_proteomes_2026_03.csv \
    --uniprot_index /bigdata/stajichlab/shared/db/Uniprot/Fungi_2026_03/novinvenio_index/v1 \
    --outdir $SCRATCH/uniprot_index_build -work-dir $SCRATCH/uniprot_index_work
```

  Acceptance 1: `manifest.json` has `n_proteomes == 1526` (or lists the missing ones). Record `n_records`, the `seq_index.sqlite` size and the total `UNIPROT_PARSE_CHUNK` realtime from the trace.

- [ ] **Step 2: Link the sordariales proteomes** without a pipeline rerun (cheap, standalone):

```bash
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
IDX=/bigdata/stajichlab/shared/db/Uniprot/Fungi_2026_03/novinvenio_index/v1
OUT=$SCRATCH/uniprot_accept; mkdir -p $OUT
tail -n +2 $NII/studies/fungi/sordariales_shallow/config.csv | while IFS=, read -r grp sp st prot dna gff short rest; do
  fa=$(find -L $NII/studies/fungi/sordariales_shallow/data_dir -name "$prot" | head -1)
  pixi run python bin/uniprot_link.py --index $IDX --protein-fasta "$fa" --short "$short" \
      --species "$sp" --output $OUT/$short.uniprot_link.tsv
done 2>&1 | tee $OUT/coverage.txt
grep Q7S6W2 $OUT/Neurcras.uniprot_link.tsv | cut -f1,2,10,12
```

  Check the config column order with `head -1 config.csv` first, and adjust the `read` field list to match.
  Expected: the Q7S6W2 row has `uniprot_match=id`, `uniprot_alphafold_id=Q7S6W2`, and xrefs that include `VEuPathDB:FungiDB:NCU05603` and PANTHER/OrthoDB/STRING entries if its record has them. Compare with `zgrep -A200 "^AC   Q7S6W2" .../UP000001805_367110.dat.gz | grep "^DR"`.

- [ ] **Step 3: Regenerate both sordariales novelty pages** from the existing results with the new columns (outputs in `$SCRATCH`, never in either repo):

```bash
for s in sordariales_shallow sordariales_shallow_cluster; do
  R=$NII/results/$s
  pixi run python bin/annotate_presence_matrix.py --matrix $R/presence_matrix.tsv \
      --output $OUT/$s.presence_matrix.function.tsv --uniprot_xref_files $OUT/*.uniprot_link.tsv
  pixi run python bin/make_report.py --matrix $OUT/$s.presence_matrix.function.tsv \
      --config $NII/studies/fungi/sordariales_shallow/config.csv \
      --tblastn_summary $R/tblastn_summary.tsv --novelties $R/novelties.*.tsv \
      --candidates_fa $R/candidates.fa --output $OUT/$s.novelties.html
done
```

  Check the required flags of `annotate_presence_matrix.py` with `--help`, and pass the same Pfam/SwissProt inputs the original run used, if any. The sordariales runs had none.
  Acceptance 2: open `$OUT/sordariales_shallow_cluster.novelties.html`, select Q7S6W2, and confirm the UniProt, AlphaFold and FungiDB links and the Publications field. With jsdom, extend a throwaway copy of `drive_reports.mjs` pointed at this file. Delete it afterwards.

- [ ] **Step 4: Coverage on a non-UniProt study** (acceptance 3). Pick one config whose protein IDs are BFD/funannotate gene IDs, for example a `configs/<focal>_<batch>.csv` built by `bin/build_targeted_configs.py`. Run Step 2's loop for it and record the per-species `seq_own/seq_other/none` counts from stderr. Report the measured numbers. Do not assume a value.

- [ ] **Step 5: Record and open PRs.**
  - Append the measured build cost and coverage to `.living/learnings.md`.
  - Open the phase PRs (A, then B, then C) with `Closes #N`, each linking the spec and plan. Confirm with the user before the first push.
  - Put the Task 12 numbers in PR D's body, or in the issue D comment if D has no code.
