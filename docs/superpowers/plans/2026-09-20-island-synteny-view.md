# Island Synteny View — Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Render, for each selected accessory island, a presence/absence grid — rows = strains (collapsed into distinct haplotypes), columns = member families in locus order — as one self-contained HTML page, so deletion breakpoints in a gain/loss island become visible.

**Architecture:** Three pure-Python layers with no I/O at import time, matching the existing report pattern. `lib/pfam_classes.py` classifies Pfam domain names into display classes. `lib/island_synteny.py` turns the pipeline's TSVs into a compact JSON payload (island selection, haplotype collapsing, bitpacking). `lib/island_synteny_template.py` holds the HTML page, assembled from the shared fragments in `lib/report_common.py` and `lib/skins.py`. `bin/pangenome_island_synteny.py` is the CLI; `modules/pangenome/island_synteny.nf` wires it into `PANGENOME_PROFILE`.

**Tech Stack:** Python 3.12 (pixi), pytest, Nextflow DSL2, vanilla JS + `<canvas>` (no dependencies — the page must open from `file://`).

**Spec:** `docs/superpowers/specs/2026-09-19-pangenome-gainloss-visualization-design.md` (on `main`)

**Issue:** #116 — **Branch:** `116-island-synteny-view` — **Worktree:** `/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/116-island-synteny-view`

## Global Constraints

- **Self-contained page.** No CDN scripts, external stylesheets, fonts, or fetches. The page must open from `file://` offline. Inline everything.
- **Untrusted strings.** Protein IDs, family IDs and Pfam descriptions come from FASTA headers and HMM output. Insert them with `textContent` (the `el()` helper), never `innerHTML`. Escape `</` in the JSON payload so an annotation cannot close the `<script>` block early.
- **Colour lives in `lib/skins.py`, never in a template.** A raw hex in a template is a test failure (`tests/test_report_templates.py::test_report_has_no_hardcoded_theme_colours`). Emit `skins_css()`.
- **Compressed inputs.** Read every pipeline TSV through `lib/compressed_io.py`'s `open_maybe_compressed()`. As of PR #113 `family_positions.tsv` publishes as `.zst`. `presence_matrix.rescued.tsv` is still plain — the helper handles both.
- **Payload is single-run.** Gene-family IDs are not stable across pipeline re-runs (73% ID overlap between two runs of the same study; a family's ID is its mmseqs representative). The payload must carry the run's own project name, and nothing may join it to another run by family ID.
- **Degrade gracefully.** Zero islands after filtering produces a valid page stating so and exit code 0 — never a crash. Mirrors `bin/pangenome_detect_trans_modules.py`'s zero-trans-edge behaviour.
- **No results in this repo.** Any manual run must pass `--outdir` outside this checkout.
- **Python style.** `bin/` scripts: `#!/usr/bin/env python3`, `chmod +x`, argparse flags only (never positional — Nextflow reorders tokens), import shared logic via `sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))`.

## Input data shapes (verified against the real `genus_vs_ureesii` run)

`report_tables/islands_with_domains.tsv` — 13 columns:
`n_strains, example_strain, island_size, member_families, n_supporting_pairs, classifications, pfam_domains, locus_id, locus_contig, locus_start, locus_end, n_members_with_coordinates, n_contigs_in_locus`

- `member_families` — comma-separated **family IDs**. A family ID is its mmseqs representative's sequence ID, so it looks like `UTAH_20380X10|CE3D476_008562-T1`. It is an opaque key, not a strain indicator.
- `pfam_domains` — comma-separated domain **names** (e.g. `AAA,DEAD,Peptidase_M41`), island-level, not per-family.
- `locus_id` — `strain:contig:start-end`.

`presence_matrix.rescued.tsv` — `family` column then one column per strain; cells are `present`/`absent`. Load with `lib.pangenome_matrix.PresenceMatrix.from_tsv(path)`; use `.is_present(family, strain)` and `.strains`.

`family_positions.tsv.zst` — `Short, family, contig, rank`. Gives locus order within a strain.

---

### Task 1: Pfam domain class map

**Files:**
- Create: `lib/pfam_classes.py`
- Test: `tests/test_pfam_classes.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `classify_domain(name: str) -> str` returning one of `"nlr"`, `"secondary_metabolite"`, `"transporter"`, `"other"`, `"unannotated"`; `dominant_class(names: list[str]) -> str`; `CLASS_ORDER: tuple[str, ...]` (precedence, most specific first); `CLASS_LABELS: dict[str, str]` (display names).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_pfam_classes.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pfam_classes import CLASS_ORDER, classify_domain, dominant_class


def test_nlr_domains_classify_as_nlr():
    assert classify_domain("NACHT") == "nlr"
    assert classify_domain("HET") == "nlr"
    assert classify_domain("Ank_2") == "nlr"


def test_secondary_metabolite_domains():
    assert classify_domain("ketoacyl-synt") == "secondary_metabolite"
    assert classify_domain("AMP-binding") == "secondary_metabolite"
    assert classify_domain("Thioesterase") == "secondary_metabolite"


def test_transporter_domains():
    assert classify_domain("MFS_1") == "transporter"
    assert classify_domain("ABC_tran") == "transporter"


def test_matching_is_case_insensitive():
    # Pfam name casing is not stable across releases.
    assert classify_domain("nacht") == "nlr"
    assert classify_domain("Ketoacyl-Synt") == "secondary_metabolite"


def test_unrecognised_domain_is_other():
    assert classify_domain("DUF1653") == "other"


def test_empty_or_missing_domain_is_unannotated():
    assert classify_domain("") == "unannotated"
    assert classify_domain("-") == "unannotated"


def test_dominant_class_is_most_frequent():
    assert dominant_class(["MFS_1", "ABC_tran", "NACHT"]) == "transporter"


def test_dominant_class_breaks_ties_by_precedence():
    # One each: precedence decides, and it must not depend on input order.
    assert dominant_class(["MFS_1", "NACHT"]) == "nlr"
    assert dominant_class(["NACHT", "MFS_1"]) == "nlr"
    assert CLASS_ORDER.index("nlr") < CLASS_ORDER.index("transporter")


def test_dominant_class_of_nothing_is_unannotated():
    assert dominant_class([]) == "unannotated"
    assert dominant_class(["-"]) == "unannotated"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/116-island-synteny-view && pixi run -e default python -m pytest tests/test_pfam_classes.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'pfam_classes'`

- [ ] **Step 3: Write minimal implementation**

```python
# lib/pfam_classes.py
"""Map a Pfam domain NAME to a display class for the island synteny view.

This is data, not logic: one curated keyword -> class table, matched
case-insensitively as a substring of the domain name, because Pfam names
carry version/variant suffixes (Ank_2, Ank_4, ...) and casing is not stable
across Pfam releases.

A family whose domains fall in two classes takes the first match in
CLASS_ORDER, which is a documented precedence, not an accident of dict
iteration order. The order is most-specific-biology first: an island
carrying both an NLR domain and a transporter is far more interesting as an
NLR island, and transporter domains are the more promiscuous of the two.
"""
from __future__ import annotations

UNANNOTATED = "unannotated"
OTHER = "other"

# Precedence, most specific first. OTHER and UNANNOTATED are terminal.
CLASS_ORDER: tuple[str, ...] = (
    "nlr",
    "secondary_metabolite",
    "transporter",
    OTHER,
    UNANNOTATED,
)

CLASS_LABELS: dict[str, str] = {
    "nlr": "NLR / incompatibility",
    "secondary_metabolite": "Secondary metabolite",
    "transporter": "Transporter",
    OTHER: "Other",
    UNANNOTATED: "No Pfam hit",
}

# keyword (lowercase, substring-matched) -> class
_KEYWORDS: dict[str, str] = {
    "nacht": "nlr",
    "het": "nlr",
    "ank": "nlr",
    "tpr": "nlr",
    "ketoacyl-synt": "secondary_metabolite",
    "amp-binding": "secondary_metabolite",
    "pp-binding": "secondary_metabolite",
    "thioesterase": "secondary_metabolite",
    "acp_": "secondary_metabolite",
    "ps-dh": "secondary_metabolite",
    "mfs_": "transporter",
    "sugar_tr": "transporter",
    "abc_tran": "transporter",
    "mfs_1": "transporter",
}

_MISSING = {"", "-", "na", "none"}


def classify_domain(name: str) -> str:
    """Class for one Pfam domain name. Missing/sentinel -> UNANNOTATED."""
    if name is None:
        return UNANNOTATED
    key = name.strip().lower()
    if key in _MISSING:
        return UNANNOTATED
    matches = {cls for kw, cls in _KEYWORDS.items() if kw in key}
    if not matches:
        return OTHER
    for cls in CLASS_ORDER:
        if cls in matches:
            return cls
    return OTHER


def dominant_class(names: list[str]) -> str:
    """Most frequent class among `names`, ties broken by CLASS_ORDER.

    OTHER and UNANNOTATED never win over a real class: they are what you get
    when there is nothing better, so they are only returned if no name in the
    list classified into a real class.
    """
    classes = [classify_domain(n) for n in (names or [])]
    real = [c for c in classes if c not in (OTHER, UNANNOTATED)]
    pool = real or classes
    if not pool:
        return UNANNOTATED
    counts: dict[str, int] = {}
    for c in pool:
        counts[c] = counts.get(c, 0) + 1
    best = max(counts.values())
    tied = [c for c, n in counts.items() if n == best]
    for cls in CLASS_ORDER:
        if cls in tied:
            return cls
    return UNANNOTATED
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run -e default python -m pytest tests/test_pfam_classes.py -q`
Expected: PASS (9 tests)

- [ ] **Step 5: Lint and commit**

```bash
pixi run lint
git add lib/pfam_classes.py tests/test_pfam_classes.py
git commit -m "lib: add Pfam domain class map for the island synteny view (#116)"
```

---

### Task 2: Island selection

**Files:**
- Create: `lib/island_synteny.py`
- Test: `tests/test_island_synteny.py`

**Interfaces:**
- Consumes: nothing from Task 1.
- Produces: `select_islands(rows: list[dict], min_strains: int = 2, top_n: int = 50) -> list[dict]`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_island_synteny.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from island_synteny import select_islands


def island(locus_id, size, n_strains):
    return {"locus_id": locus_id, "island_size": str(size),
            "n_strains": str(n_strains), "member_families": "a,b",
            "pfam_domains": "-"}


def test_single_strain_islands_are_excluded_by_default():
    # 62% of located islands on a real run are single-strain; their grid is
    # one filled row and N-1 empty ones, so there is no breakpoint to show.
    rows = [island("S1:c1:1-9000", 69, 1), island("S2:c2:1-400", 5, 3)]
    kept = select_islands(rows)
    assert [r["locus_id"] for r in kept] == ["S2:c2:1-400"]


def test_islands_are_ranked_by_size_descending():
    rows = [island("small", 3, 2), island("big", 40, 2), island("mid", 9, 2)]
    assert [r["locus_id"] for r in select_islands(rows)] == ["big", "mid", "small"]


def test_top_n_truncates_after_filtering_not_before():
    rows = [island("private", 99, 1), island("keep1", 8, 2), island("keep2", 7, 2)]
    kept = select_islands(rows, top_n=2)
    assert [r["locus_id"] for r in kept] == ["keep1", "keep2"]


def test_min_strains_is_configurable():
    rows = [island("two", 5, 2), island("four", 4, 4)]
    assert [r["locus_id"] for r in select_islands(rows, min_strains=4)] == ["four"]


def test_rows_without_a_locus_are_excluded():
    rows = [island("-", 40, 5), island("S1:c1:1-400", 5, 5)]
    assert [r["locus_id"] for r in select_islands(rows)] == ["S1:c1:1-400"]


def test_malformed_counts_do_not_crash():
    rows = [island("bad", "", ""), island("ok", 5, 2)]
    assert [r["locus_id"] for r in select_islands(rows)] == ["ok"]


def test_no_islands_returns_empty_list():
    assert select_islands([]) == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e default python -m pytest tests/test_island_synteny.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'island_synteny'`

- [ ] **Step 3: Write minimal implementation**

```python
# lib/island_synteny.py
"""Build the island synteny view's JSON payload.

For one accessory island: rows = strains, columns = member families in locus
order, filled cell = family present in that strain. Deletion breakpoints show
up as vertical edges.

No I/O at import time; every function here is pure and takes already-parsed
rows, so the whole payload is unit-testable without touching a pipeline run.
"""
from __future__ import annotations


def _as_int(value, default: int = 0) -> int:
    """Parse a TSV cell to int, tolerating '' and the '-' missing sentinel."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def select_islands(rows: list[dict], min_strains: int = 2,
                   top_n: int = 50) -> list[dict]:
    """Islands worth drawing, largest first.

    Two filters, in this order:

    1. A locus must be resolved (`locus_id != '-'`), since the view is
       locus-ordered and an unplaced island has no column order.
    2. The island must be carried by at least `min_strains` strains. On the
       real genus_vs_ureesii run 17,218 of 27,836 located islands (62%) are
       single-strain, and all 20 of the largest were -- an unfiltered size
       ranking is therefore entirely strain-private content, and a
       single-strain island's grid is one filled row with no breakpoint in
       it. This is the same filter PR #110 added to report.md's top-islands
       table, for the same reason.

    Truncation to `top_n` happens AFTER filtering, so the filter never costs
    slots in the output.
    """
    located = [r for r in rows if r.get("locus_id", "-") not in ("-", "", None)]
    kept = [r for r in located if _as_int(r.get("n_strains")) >= min_strains]
    kept.sort(key=lambda r: -_as_int(r.get("island_size")))
    return kept[:top_n]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run -e default python -m pytest tests/test_island_synteny.py -q`
Expected: PASS (7 tests)

- [ ] **Step 5: Commit**

```bash
git add lib/island_synteny.py tests/test_island_synteny.py
git commit -m "lib: island selection for the synteny view, excluding strain-private islands (#116)"
```

---

### Task 3: Haplotype collapsing

**Files:**
- Modify: `lib/island_synteny.py`
- Test: `tests/test_island_synteny.py`

**Interfaces:**
- Consumes: nothing from Tasks 1-2.
- Produces: `collapse_haplotypes(presence_rows: dict[str, list[bool]]) -> list[dict]`, each dict `{"pattern": str, "count": int, "strains": list[str]}`, sorted by pattern string ascending.

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_island_synteny.py
from island_synteny import collapse_haplotypes


def test_identical_rows_collapse_into_one_haplotype():
    # 530 strains x a 20-56 gene island is unreadable as 530 rows, but the
    # number of DISTINCT presence patterns is typically tens.
    rows = {"s1": [True, False], "s2": [True, False], "s3": [False, True]}
    haps = collapse_haplotypes(rows)
    assert len(haps) == 2
    by_pattern = {h["pattern"]: h for h in haps}
    assert by_pattern["10"]["count"] == 2
    assert by_pattern["10"]["strains"] == ["s1", "s2"]
    assert by_pattern["01"]["count"] == 1


def test_haplotypes_sort_lexicographically_on_the_pattern():
    # Lexicographic order on the 0/1 string groups related haplotypes, which
    # is what makes breakpoints line up visually down the grid.
    rows = {"a": [True, True], "b": [False, False], "c": [False, True]}
    assert [h["pattern"] for h in collapse_haplotypes(rows)] == ["00", "01", "11"]


def test_strains_within_a_haplotype_are_sorted():
    rows = {"zeta": [True], "alpha": [True]}
    assert collapse_haplotypes(rows)[0]["strains"] == ["alpha", "zeta"]


def test_counts_sum_to_the_number_of_strains():
    rows = {f"s{i}": [i % 2 == 0] for i in range(10)}
    assert sum(h["count"] for h in collapse_haplotypes(rows)) == 10


def test_no_strains_gives_no_haplotypes():
    assert collapse_haplotypes({}) == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e default python -m pytest tests/test_island_synteny.py -q -k haplotype`
Expected: FAIL — `ImportError: cannot import name 'collapse_haplotypes'`

- [ ] **Step 3: Write minimal implementation**

```python
# append to lib/island_synteny.py

def collapse_haplotypes(presence_rows: dict[str, list[bool]]) -> list[dict]:
    """Collapse identical strain rows into distinct presence patterns.

    Returns one entry per distinct pattern: {"pattern", "count", "strains"},
    sorted lexicographically on the 0/1 pattern string. That sort is
    equivalent to a Hamming-distance leaf ordering for this data and is what
    makes deletion breakpoints align into visible vertical edges.

    Row collapsing is the load-bearing decision of this view: without it, a
    530-strain study draws 530 rows for a 20-56 column island and nothing is
    legible. With it, the row count is the number of haplotypes, typically
    tens.
    """
    by_pattern: dict[str, list[str]] = {}
    for strain, vector in presence_rows.items():
        pattern = "".join("1" if v else "0" for v in vector)
        by_pattern.setdefault(pattern, []).append(strain)
    return [
        {"pattern": p, "count": len(s), "strains": sorted(s)}
        for p, s in sorted(by_pattern.items())
    ]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run -e default python -m pytest tests/test_island_synteny.py -q`
Expected: PASS (12 tests)

- [ ] **Step 5: Commit**

```bash
git add lib/island_synteny.py tests/test_island_synteny.py
git commit -m "lib: collapse identical strain rows into haplotypes (#116)"
```

---

### Task 4: Locus ordering of member families

**Files:**
- Modify: `lib/island_synteny.py`
- Test: `tests/test_island_synteny.py`

**Interfaces:**
- Consumes: nothing from Tasks 1-3.
- Produces: `order_families_by_locus(families: list[str], positions: dict[tuple[str, str], int], strain: str) -> list[str]`, where `positions` maps `(strain, family) -> rank`.

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_island_synteny.py
from island_synteny import order_families_by_locus


def test_families_are_ordered_by_rank_in_the_example_strain():
    positions = {("S1", "famC"): 3, ("S1", "famA"): 1, ("S1", "famB"): 2}
    assert order_families_by_locus(["famA", "famB", "famC"], positions, "S1") == \
        ["famA", "famB", "famC"]


def test_input_order_does_not_matter():
    positions = {("S1", "famC"): 3, ("S1", "famA"): 1, ("S1", "famB"): 2}
    assert order_families_by_locus(["famC", "famA", "famB"], positions, "S1") == \
        ["famA", "famB", "famC"]


def test_families_without_a_position_sort_last_in_stable_name_order():
    # A member with no coordinate in this strain still deserves a column --
    # dropping it would silently shrink the island.
    positions = {("S1", "famB"): 2}
    assert order_families_by_locus(["famZ", "famB", "famA"], positions, "S1") == \
        ["famB", "famA", "famZ"]


def test_positions_from_another_strain_are_ignored():
    positions = {("S2", "famA"): 1, ("S1", "famB"): 5}
    assert order_families_by_locus(["famA", "famB"], positions, "S1") == \
        ["famB", "famA"]


def test_empty_family_list_gives_empty_order():
    assert order_families_by_locus([], {}, "S1") == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e default python -m pytest tests/test_island_synteny.py -q -k locus`
Expected: FAIL — `ImportError: cannot import name 'order_families_by_locus'`

- [ ] **Step 3: Write minimal implementation**

```python
# append to lib/island_synteny.py

def order_families_by_locus(families: list[str],
                            positions: dict[tuple[str, str], int],
                            strain: str) -> list[str]:
    """Member families in locus order within `strain`.

    Gene ORDER is what makes this view meaningful -- the question it answers
    is "where does the block break?" -- so the columns must follow the
    island's own layout, taken from the example strain's gene ranks in
    family_positions.tsv.

    A family with no position in this strain keeps a column, sorted last by
    name. Dropping it would silently shrink the island and change the
    breakpoint the reader sees.
    """
    ranked = [(positions[(strain, f)], f) for f in families
              if (strain, f) in positions]
    unranked = sorted(f for f in families if (strain, f) not in positions)
    return [f for _, f in sorted(ranked)] + unranked
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run -e default python -m pytest tests/test_island_synteny.py -q`
Expected: PASS (17 tests)

- [ ] **Step 5: Commit**

```bash
git add lib/island_synteny.py tests/test_island_synteny.py
git commit -m "lib: order island member families by locus rank (#116)"
```

---

### Task 5: Payload assembly

**Files:**
- Modify: `lib/island_synteny.py`
- Test: `tests/test_island_synteny.py`

**Interfaces:**
- Consumes: `classify_domain`/`dominant_class` (Task 1), `select_islands` (Task 2), `collapse_haplotypes` (Task 3), `order_families_by_locus` (Task 4).
- Produces: `build_payload(island_rows, matrix, positions, project, min_strains=2, top_n=50) -> dict` with keys `project`, `generated`, `islands`, `classes`, `n_islands_total`, `n_islands_excluded`. `matrix` is a `lib.pangenome_matrix.PresenceMatrix`.

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_island_synteny.py
from island_synteny import build_payload
from pangenome_matrix import PresenceMatrix


def make_matrix(calls):
    """calls: {family: {strain: bool}} -> a PresenceMatrix."""
    m = PresenceMatrix()
    for family, per_strain in calls.items():
        for strain, present in per_strain.items():
            m.set_call(family, strain, "present" if present else "absent")
    return m


def test_payload_carries_one_island_with_ordered_columns_and_haplotypes():
    rows = [{"locus_id": "S1:c1:1-400", "island_size": "2", "n_strains": "2",
             "member_families": "famB,famA", "pfam_domains": "NACHT,MFS_1",
             "example_strain": "S1"}]
    matrix = make_matrix({
        "famA": {"S1": True, "S2": True, "S3": False},
        "famB": {"S1": True, "S2": False, "S3": False},
    })
    positions = {("S1", "famA"): 1, ("S1", "famB"): 2}
    payload = build_payload(rows, matrix, positions, project="demo")

    assert payload["project"] == "demo"
    assert len(payload["islands"]) == 1
    island = payload["islands"][0]
    assert island["families"] == ["famA", "famB"]
    assert island["locus_id"] == "S1:c1:1-400"
    patterns = {h["pattern"]: h["count"] for h in island["haplotypes"]}
    assert patterns == {"11": 1, "10": 1, "00": 1}


def test_payload_reports_how_many_islands_were_excluded():
    rows = [
        {"locus_id": "S1:c1:1-9", "island_size": "9", "n_strains": "1",
         "member_families": "famA", "pfam_domains": "-", "example_strain": "S1"},
        {"locus_id": "S2:c2:1-9", "island_size": "2", "n_strains": "2",
         "member_families": "famA", "pfam_domains": "-", "example_strain": "S2"},
    ]
    matrix = make_matrix({"famA": {"S1": True, "S2": True}})
    payload = build_payload(rows, matrix, {}, project="demo")
    assert payload["n_islands_total"] == 2
    assert payload["n_islands_excluded"] == 1
    assert len(payload["islands"]) == 1


def test_each_family_carries_its_pfam_class_for_the_glyph_strip():
    rows = [{"locus_id": "S1:c1:1-400", "island_size": "1", "n_strains": "2",
             "member_families": "famA", "pfam_domains": "NACHT",
             "example_strain": "S1"}]
    matrix = make_matrix({"famA": {"S1": True, "S2": True}})
    payload = build_payload(rows, matrix, {}, project="demo")
    assert payload["islands"][0]["dominant_class"] == "nlr"
    assert "nlr" in payload["classes"]


def test_zero_islands_produces_a_valid_empty_payload():
    # The single-species-ingroup case. A valid page saying so, not a crash.
    payload = build_payload([], make_matrix({}), {}, project="demo")
    assert payload["islands"] == []
    assert payload["n_islands_total"] == 0
    assert payload["project"] == "demo"


def test_families_absent_from_the_matrix_are_scored_absent_not_dropped():
    rows = [{"locus_id": "S1:c1:1-400", "island_size": "2", "n_strains": "2",
             "member_families": "famA,ghost", "pfam_domains": "-",
             "example_strain": "S1"}]
    matrix = make_matrix({"famA": {"S1": True, "S2": True}})
    payload = build_payload(rows, matrix, {}, project="demo")
    assert payload["islands"][0]["families"] == ["famA", "ghost"]
    for hap in payload["islands"][0]["haplotypes"]:
        assert len(hap["pattern"]) == 2
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e default python -m pytest tests/test_island_synteny.py -q -k payload`
Expected: FAIL — `ImportError: cannot import name 'build_payload'`

- [ ] **Step 3: Write minimal implementation**

```python
# add near the top of lib/island_synteny.py, after the existing imports
from datetime import datetime, timezone

from pfam_classes import CLASS_LABELS, dominant_class

# ... then append:

def build_payload(island_rows: list[dict], matrix, positions: dict,
                  project: str, min_strains: int = 2,
                  top_n: int = 50) -> dict:
    """The JSON payload embedded in the island synteny page.

    `matrix` is a lib.pangenome_matrix.PresenceMatrix; `positions` maps
    (strain, family) -> gene rank.

    The payload records `project` because family IDs are NOT stable across
    pipeline re-runs -- a family's ID is its mmseqs cluster representative,
    and mmseqs does not pick the same one twice (measured: 73% ID overlap
    between two runs of the same study). The page is only ever valid for the
    run that produced it, and must never be joined to another run's output by
    family ID.
    """
    selected = select_islands(island_rows, min_strains=min_strains, top_n=top_n)
    located = [r for r in island_rows
               if r.get("locus_id", "-") not in ("-", "", None)]
    strains = sorted(matrix.strains)

    islands = []
    for row in selected:
        families = [f for f in row.get("member_families", "").split(",") if f]
        example = row.get("example_strain", "")
        ordered = order_families_by_locus(families, positions, example)
        presence_rows = {
            s: [matrix.is_present(f, s) for f in ordered] for s in strains
        }
        domains = [d for d in row.get("pfam_domains", "").split(",") if d]
        islands.append({
            "locus_id": row.get("locus_id", "-"),
            "locus_contig": row.get("locus_contig", ""),
            "locus_start": _as_int(row.get("locus_start"), -1),
            "locus_end": _as_int(row.get("locus_end"), -1),
            "example_strain": example,
            "size": len(ordered),
            "n_strains": _as_int(row.get("n_strains")),
            "families": ordered,
            "domains": domains,
            "dominant_class": dominant_class(domains),
            "haplotypes": collapse_haplotypes(presence_rows),
        })

    return {
        "project": project,
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "strains": strains,
        "classes": CLASS_LABELS,
        "islands": islands,
        "n_islands_total": len(located),
        "n_islands_excluded": len(located) - len(selected),
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run -e default python -m pytest tests/test_island_synteny.py -q`
Expected: PASS (22 tests)

- [ ] **Step 5: Lint and commit**

```bash
pixi run lint
git add lib/island_synteny.py tests/test_island_synteny.py
git commit -m "lib: assemble the island synteny payload (#116)"
```

---

### Task 6: The HTML page

**Files:**
- Create: `lib/island_synteny_template.py`
- Modify: `tests/test_report_templates.py:38-44` (add the page to `PAGES`)
- Test: `tests/test_report_templates.py`

**Interfaces:**
- Consumes: `lib/report_common.py`'s `BASE_PAGE_CSS`, `SKIN_VARS_CSS`, `LOGO_CSS`, `FOOTER_CSS`, `FOOTER_HTML`, `FAVICON_LINK_HTML`, `LOGO_IMG_HTML`, `EL_HELPER_JS`; `lib/skins.py`'s `skin_picker_html()`, `skin_boot_js()`, `skin_picker_js()`.
- Produces: `ISLAND_SYNTENY_TEMPLATE: str`, with substitution markers `__PROJECT_TITLE__` and `/*__PAYLOAD__*/`.

Follow `lib/core_report_template.py` as the model — it is the smallest of the existing single-table pages and uses the same fragment imports. Read it first.

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_report_templates.py, and add
#   'island_synteny': ISLAND_SYNTENY_TEMPLATE
# to the PAGES dict at line 38, with the matching import at the top.

def test_island_synteny_page_has_substitution_markers():
    from island_synteny_template import ISLAND_SYNTENY_TEMPLATE
    assert '__PROJECT_TITLE__' in ISLAND_SYNTENY_TEMPLATE
    assert '/*__PAYLOAD__*/' in ISLAND_SYNTENY_TEMPLATE


def test_island_synteny_page_is_self_contained():
    from island_synteny_template import ISLAND_SYNTENY_TEMPLATE
    # The page is copied off the cluster and opened from file://.
    for forbidden in ('http://', 'https://cdn', '<link rel="stylesheet" href="http'):
        assert forbidden not in ISLAND_SYNTENY_TEMPLATE.replace(
            'https://github.com/stajichlab/NovInvenio', '')


def test_island_synteny_page_emits_skin_tokens():
    from island_synteny_template import ISLAND_SYNTENY_TEMPLATE
    assert '--bg' in ISLAND_SYNTENY_TEMPLATE
    assert 'data-theme' in ISLAND_SYNTENY_TEMPLATE
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e default python -m pytest tests/test_report_templates.py -q -k island_synteny`
Expected: FAIL — `ModuleNotFoundError: No module named 'island_synteny_template'`

- [ ] **Step 3: Write the template**

Build `lib/island_synteny_template.py` on the `lib/core_report_template.py` skeleton. Requirements, all from the spec:

- A **sidebar** listing the selected islands: locus id, size (families), distinct-haplotype count, carrying-strain count, and a dominant-Pfam-class chip. Sortable and filterable. Selecting one redraws a single canvas.
- A **canvas grid**: rows = haplotypes (with a `xN` count badge), columns = families in locus order, filled cell = present. Read colours through `getComputedStyle` in a `palette()` helper so a skin change repaints — wire `window.onSkinChange` exactly as `lib/report_template.py` does.
- A **glyph strip above the grid**: one coloured tick per column keyed to the family's Pfam class; hover reveals the full domain string and family ID. Annotation never goes in the cells.
- A **row-sort `<select>`**: `pattern` (default), `count`, `strain name`. **Do not offer a phylogeny option** — this pipeline has no strain tree, and offering it would imply a capability that does not exist (see `todo/pangenome-phylogeny-aware-gain-loss.md`).
- An **empty state**: when `payload.islands` is empty, show a message naming `n_islands_total` and `n_islands_excluded` rather than an empty canvas.
- A **note line** stating how many islands were excluded and why, mirroring `report.md`'s wording.
- All text inserted with the `el()` helper (`EL_HELPER_JS`), never `innerHTML`.

- [ ] **Step 4: Run the template tests**

Run: `pixi run -e default python -m pytest tests/test_report_templates.py -q`
Expected: PASS, including `test_page_javascript_parses[island_synteny]` (the `node --check` pass) and `test_report_has_no_hardcoded_theme_colours`.

If node is unavailable the JS check skips — verify by hand before committing:

```bash
pixi run -e default python -c "
from lib.island_synteny_template import ISLAND_SYNTENY_TEMPLATE
import re
m = re.findall(r'<script>\n(.*?)</script>', ISLAND_SYNTENY_TEMPLATE, re.S)
open('/tmp/isv.js','w').write(m[-1])
"
node --check /tmp/isv.js && rm /tmp/isv.js
```

- [ ] **Step 5: Commit**

```bash
git add lib/island_synteny_template.py tests/test_report_templates.py
git commit -m "lib: island synteny HTML page (#116)"
```

---

### Task 7: The CLI

**Files:**
- Create: `bin/pangenome_island_synteny.py`
- Test: `tests/test_pangenome_island_synteny.py`

**Interfaces:**
- Consumes: `build_payload` (Task 5), `ISLAND_SYNTENY_TEMPLATE` (Task 6), `PresenceMatrix.from_tsv`, `open_maybe_compressed`.
- Produces: a CLI with flags `--islands_with_domains`, `--presence_matrix`, `--family_positions`, `--project`, `--min_strains`, `--top_islands`, `--output`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_pangenome_island_synteny.py
import json
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).parent.parent / "bin"


def write_inputs(tmp_path):
    islands = tmp_path / "islands_with_domains.tsv"
    islands.write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\t"
        "n_supporting_pairs\tclassifications\tpfam_domains\tlocus_id\t"
        "locus_contig\tlocus_start\tlocus_end\tn_members_with_coordinates\t"
        "n_contigs_in_locus\n"
        "2\tS1\t2\tfamA,famB\t1\tunexplained_physical\tNACHT\t"
        "S1:c1:100-400\tc1\t100\t400\t2\t1\n"
        "1\tS2\t9\tfamC\t1\tunexplained_physical\t-\t"
        "S2:c2:1-900\tc2\t1\t900\t1\t1\n"
    )
    matrix = tmp_path / "presence_matrix.tsv"
    matrix.write_text(
        "family\tS1\tS2\nfamA\tpresent\tpresent\nfamB\tpresent\tabsent\n"
        "famC\tabsent\tpresent\n"
    )
    positions = tmp_path / "family_positions.tsv"
    positions.write_text("Short\tfamily\tcontig\trank\nS1\tfamA\tc1\t1\n"
                         "S1\tfamB\tc1\t2\n")
    return islands, matrix, positions


def run_cli(tmp_path, *extra):
    islands, matrix, positions = write_inputs(tmp_path)
    out = tmp_path / "island_synteny.html"
    subprocess.run(
        [sys.executable, str(BIN / "pangenome_island_synteny.py"),
         "--islands_with_domains", str(islands),
         "--presence_matrix", str(matrix),
         "--family_positions", str(positions),
         "--project", "demo", "--output", str(out), *extra],
        check=True,
    )
    return out


def payload_of(html_path):
    html = html_path.read_text()
    start = html.index('<script type="application/json"')
    start = html.index(">", start) + 1
    return json.loads(html[start:html.index("</script>", start)])


def test_cli_writes_a_self_contained_page(tmp_path):
    out = run_cli(tmp_path)
    html = out.read_text()
    assert html.lstrip().startswith("<!DOCTYPE html")
    assert "demo" in html


def test_cli_excludes_the_single_strain_island(tmp_path):
    payload = payload_of(run_cli(tmp_path))
    assert [i["locus_id"] for i in payload["islands"]] == ["S1:c1:100-400"]
    assert payload["n_islands_excluded"] == 1


def test_cli_orders_columns_by_locus_rank(tmp_path):
    payload = payload_of(run_cli(tmp_path))
    assert payload["islands"][0]["families"] == ["famA", "famB"]


def test_cli_min_strains_is_configurable(tmp_path):
    payload = payload_of(run_cli(tmp_path, "--min_strains", "3"))
    assert payload["islands"] == []


def test_cli_reads_zst_inputs(tmp_path):
    # Nextflow publishes family_positions.tsv as .zst since PR #113, and
    # stages it as a SYMLINK -- zstd returns an EMPTY stream for a symlink
    # without -f, which open_maybe_compressed handles.
    import shutil
    if shutil.which("zstd") is None:
        import pytest
        pytest.skip("zstd not available")
    islands, matrix, positions = write_inputs(tmp_path)
    subprocess.run(["zstd", "-q", "-f", str(positions),
                    "-o", str(positions) + ".zst"], check=True)
    link = tmp_path / "staged_positions.tsv.zst"
    link.symlink_to(Path(str(positions) + ".zst"))
    out = tmp_path / "out.html"
    subprocess.run(
        [sys.executable, str(BIN / "pangenome_island_synteny.py"),
         "--islands_with_domains", str(islands),
         "--presence_matrix", str(matrix),
         "--family_positions", str(link),
         "--project", "demo", "--output", str(out)], check=True)
    assert payload_of(out)["islands"][0]["families"] == ["famA", "famB"]


def test_cli_with_no_islands_still_writes_a_page(tmp_path):
    islands = tmp_path / "empty.tsv"
    islands.write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\t"
        "n_supporting_pairs\tclassifications\tpfam_domains\tlocus_id\t"
        "locus_contig\tlocus_start\tlocus_end\tn_members_with_coordinates\t"
        "n_contigs_in_locus\n")
    matrix = tmp_path / "m.tsv"
    matrix.write_text("family\tS1\n")
    positions = tmp_path / "p.tsv"
    positions.write_text("Short\tfamily\tcontig\trank\n")
    out = tmp_path / "empty.html"
    proc = subprocess.run(
        [sys.executable, str(BIN / "pangenome_island_synteny.py"),
         "--islands_with_domains", str(islands),
         "--presence_matrix", str(matrix),
         "--family_positions", str(positions),
         "--project", "demo", "--output", str(out)],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert payload_of(out)["islands"] == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e default python -m pytest tests/test_pangenome_island_synteny.py -q`
Expected: FAIL — the script does not exist (`No such file or directory`)

- [ ] **Step 3: Write the CLI**

```python
#!/usr/bin/env python3
"""Render the island synteny / presence-absence view (issue #116).

One self-contained HTML page: for each selected accessory island, a
presence/absence grid with rows = strains collapsed into distinct
haplotypes and columns = member families in locus order.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import open_maybe_compressed  # noqa: E402
from island_synteny import build_payload  # noqa: E402
from island_synteny_template import ISLAND_SYNTENY_TEMPLATE  # noqa: E402
from pangenome_matrix import PresenceMatrix  # noqa: E402


def load_positions(path: str) -> dict[tuple[str, str], int]:
    """(strain, family) -> gene rank, from family_positions.tsv[.zst]."""
    positions: dict[tuple[str, str], int] = {}
    with open_maybe_compressed(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                positions[(row["Short"], row["family"])] = int(row["rank"])
            except (KeyError, ValueError):
                continue
    return positions


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--islands_with_domains", required=True)
    ap.add_argument("--presence_matrix", required=True)
    ap.add_argument("--family_positions", required=True)
    ap.add_argument("--project", required=True)
    ap.add_argument("--min_strains", type=int, default=2)
    ap.add_argument("--top_islands", type=int, default=50)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    with open_maybe_compressed(args.islands_with_domains) as fh:
        island_rows = list(csv.DictReader(fh, delimiter="\t"))
    matrix = PresenceMatrix.from_tsv(args.presence_matrix)
    positions = load_positions(args.family_positions)

    payload = build_payload(
        island_rows, matrix, positions, project=args.project,
        min_strains=args.min_strains, top_n=args.top_islands,
    )

    # Escape `</` so a Pfam description or family ID cannot close the
    # <script> block early -- these strings come from HMM output and FASTA
    # headers, which this pipeline does not control.
    blob = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
    html = (ISLAND_SYNTENY_TEMPLATE
            .replace("__PROJECT_TITLE__", args.project)
            .replace("/*__PAYLOAD__*/", blob))
    Path(args.output).write_text(html)

    print(f"pangenome_island_synteny: {len(payload['islands'])} islands drawn, "
          f"{payload['n_islands_excluded']} excluded (< {args.min_strains} "
          f"strains), wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run tests and make the script executable**

```bash
chmod +x bin/pangenome_island_synteny.py
pixi run -e default python -m pytest tests/test_pangenome_island_synteny.py -q
```
Expected: PASS (6 tests)

- [ ] **Step 5: Lint and commit**

```bash
pixi run lint
git add bin/pangenome_island_synteny.py tests/test_pangenome_island_synteny.py
git commit -m "bin: island synteny view CLI (#116)"
```

---

### Task 8: Nextflow wiring

**Files:**
- Create: `modules/pangenome/island_synteny.nf`
- Modify: `workflows/pangenome_profile.nf` (include line ~49; call site after `REPORT_TABLES`)
- Modify: `nextflow.config` (params block, near `pangenome_top_islands_min_strains`)
- Modify: `pangenome.nf` (`--help` text)

**Interfaces:**
- Consumes: `bin/pangenome_island_synteny.py` (Task 7); `REPORT_TABLES.out` islands table; the rescued presence matrix; `FAMILY_POSITIONS.out.positions`.
- Produces: process `ISLAND_SYNTENY`, output `island_synteny.html`, published to `${params.outdir}/${Helpers.projectName(params)}/pangenome`.

- [ ] **Step 1: Write the module**

```groovy
// ISLAND_SYNTENY -- per-island presence/absence grid (issue #116), rows =
// strains collapsed into distinct haplotypes, columns = member families in
// locus order, so a deletion breakpoint shows up as a vertical edge.
//
// Runs inside the --pangenome_island_pfam_hmm block because the islands
// table it consumes (REPORT_TABLES' islands_with_domains.tsv) only exists
// there. Degrades to a valid "no islands" page rather than failing when a
// study has none.
process ISLAND_SYNTENY {
    label 'low_cpu'
    tag "island_synteny"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(islands_with_domains)
    path(presence_matrix)
    path(family_positions)

    output:
    path("island_synteny.html"), emit: page

    script:
    """
    pangenome_island_synteny.py \
        --islands_with_domains ${islands_with_domains} \
        --presence_matrix ${presence_matrix} \
        --family_positions ${family_positions} \
        --project ${Helpers.projectName(params)} \
        --min_strains ${params.pangenome_top_islands_min_strains} \
        --top_islands ${params.pangenome_viz_top_islands} \
        --output island_synteny.html
    """
}
```

- [ ] **Step 2: Add the parameter**

In `nextflow.config`, next to `pangenome_top_islands_min_strains`:

```groovy
    // Islands drawn in the island synteny view (issue #116). Selection is
    // top-N by size AFTER the pangenome_top_islands_min_strains filter, so
    // the filter never costs slots. 50 keeps the embedded payload small --
    // 530 strains x a 56-gene island is ~30 kB raw per island, so ~100
    // islands still fits under 1 MB.
    pangenome_viz_top_islands = 50
```

In `pangenome.nf`'s `--help` block, after `--pangenome_top_islands_min_strains`:

```
      --pangenome_viz_top_islands      Islands drawn in island_synteny.html
                                       (default: 50), selected by size after
                                       the min-strains filter above.
```

- [ ] **Step 3: Wire the call site**

In `workflows/pangenome_profile.nf`, add to the includes near line 49:

```groovy
include { ISLAND_SYNTENY }                                          from '../modules/pangenome/island_synteny'
```

and after the `REPORT_TABLES(...)` call:

```groovy
        ISLAND_SYNTENY(
            REPORT_TABLES.out.islands_with_domains,
            rescued_matrix,
            FAMILY_POSITIONS.out.positions,
        )
```

Check `REPORT_TABLES`'s `output:` block for the real emit name of the islands table; if it is not already named, add `emit: islands_with_domains` to it rather than referencing it positionally.

- [ ] **Step 4: Verify the workflow parses**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/116-island-synteny-view
pixi run nextflow run pangenome.nf --help
```
Expected: the help text prints, including `--pangenome_viz_top_islands`, with no Groovy parse error.

- [ ] **Step 5: Commit**

```bash
git add modules/pangenome/island_synteny.nf workflows/pangenome_profile.nf nextflow.config pangenome.nf
git commit -m "pangenome: wire ISLAND_SYNTENY into the profile workflow (#116)"
```

---

### Task 9: Real smoke run and CHANGES

**Files:**
- Modify: `CHANGES.md`

**Interfaces:** none — this task validates Tasks 1-8 end to end.

- [ ] **Step 1: Run the full test suite and lint**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/116-island-synteny-view
pixi run -e default python -m pytest tests/ -q
pixi run lint
```
Expected: all pass except `test_novelty_discovery_integration.py`, which fails with `MMSEQS_FAMILY_CLUSTER` exit 132 (SIGILL) on non-AVX2 nodes — a known pre-existing issue, unrelated.

- [ ] **Step 2: Submit a 5-strain smoke run**

Do NOT `nohup &` this from the session — a background run dies with the session's SLURM job, and `/scratch` is node-local and vanishes with it. Submit it, with every path on `/bigdata`:

```bash
B=/bigdata/stajichlab/jstajich/pangenome_smoke_116
mkdir -p $B
cp /bigdata/stajichlab/jstajich/pangenome_smoke_2026-09-19/config_smoke.csv $B/
cat > $B/run.sh <<'EOF'
#!/bin/bash
#SBATCH -p short -c 8 --mem=24G -t 02:00:00 -C milan
#SBATCH -J pg_smoke_116 -o /bigdata/stajichlab/jstajich/pangenome_smoke_116/slurm-%j.out -e /bigdata/stajichlab/jstajich/pangenome_smoke_116/slurm-%j.err
set -euo pipefail
cd /bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/116-island-synteny-view
pixi run nextflow run pangenome.nf \
  --pangenome_samplesheet /bigdata/stajichlab/jstajich/pangenome_smoke_116/config_smoke.csv \
  --pangenome_data_dir /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/data_dir \
  --pangenome_cluster_backend diamond \
  --pangenome_project smoke_116 \
  --outdir /bigdata/stajichlab/jstajich/pangenome_smoke_116/out \
  --pangenome_island_pfam_hmm /bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm \
  --max_cpus 8 \
  -work-dir /bigdata/stajichlab/jstajich/pangenome_smoke_116/work \
  -profile local
EOF
sbatch $B/run.sh
```

- [ ] **Step 3: Verify the page**

```bash
B=/bigdata/stajichlab/jstajich/pangenome_smoke_116
grep -E "SUCCESS|FAILED|ISLAND_SYNTENY" $B/slurm-*.out | tail -5
P=$B/out/output/pangenome/island_synteny.html
ls -l $P
node --check <(python3 -c "
import re,sys
h=open('$P').read()
print(re.findall(r'<script>\n(.*?)</script>', h, re.S)[-1])
")
```
Expected: `[SUCCESS] ... failed=0`; the page exists; the JS parses. Open it in a browser and confirm a grid draws, the sidebar lists islands, the skin picker repaints the canvas, and the row-sort select reorders rows.

- [ ] **Step 4: Update CHANGES.md**

Add a `### New: island synteny view` entry under `## Unreleased`, naming the new files, the two parameters, the single-strain exclusion and its 62% justification, and the single-run-validity caveat about family IDs.

- [ ] **Step 5: Commit and open the PR**

```bash
git add CHANGES.md
git commit -m "docs: record the island synteny view in CHANGES (#116)"
git push -u origin 116-island-synteny-view
gh pr create --base main --title "pangenome: island synteny / presence-absence view (Phase 1)" --body "Closes #116 ..."
```

---

## Self-Review

**Spec coverage (Phase 1 sections):** row collapsing → Task 3; lexicographic sort → Task 3; alternate row sorts and the deliberate absence of a phylogeny sort → Task 6; domain glyph strip → Tasks 1, 5, 6; island navigation sidebar → Task 6; selection with the min-strains filter → Tasks 2, 5; Pfam class map + precedence → Task 1; single palette/legend → Tasks 1, 6; annotation lookup → Task 5; self-contained HTML, no per-island images → Tasks 6, 7; compressed inputs → Task 7; single-run validity → Task 5; graceful zero-islands → Tasks 5, 7; Nextflow integration → Task 8. **Out of scope, per the spec:** rearrangement detection; View B.

**Gap accepted deliberately:** the spec mentions bitpacking presence arrays. Tasks 3 and 5 emit the pattern as a `"0110"` string instead. For 50 islands at 530 strains that is well inside the 10 MB target, and a string is directly usable as the collapse key and the sort key. If a real run's payload exceeds ~5 MB, bitpack in a follow-up — do not pre-optimize here.

**Type consistency:** `select_islands` / `collapse_haplotypes` / `order_families_by_locus` / `build_payload` are used in Tasks 5 and 7 with exactly the signatures defined in Tasks 2-5. `dominant_class` returns a key of `CLASS_LABELS`, which is what Task 5 puts in `payload["classes"]` and Task 6 renders. `PresenceMatrix.is_present(family, strain)` and `.strains` match the real API in `lib/pangenome_matrix.py`.
