# Island Locus View and Clinker Synteny Panel Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the island presence grid as the default view of `island_synteny.html` with an exemplar-anchored locus view (flank anchors, five per-strain cell states, row classes, a shared-breakpoint track), and add a clinker synteny figure per drawn locus.

**Architecture:** Part A adds a pure library `lib/island_locus.py` (loci, exemplar, columns, cell states, row classes, breakpoints), a CLI `bin/pangenome_island_loci.py` that streams `family_positions` three times and writes `island_loci.json`, a new process `ISLAND_LOCI`, and a locus-view fragment `lib/island_locus_template.py` that `island_synteny.html` shows by default. Part B adds per-strain regions and clinker strain choice to the same library, GenBank slicing (`lib/genbank_slice.py`, `ISLAND_GBK_SLICE`), a batched clinker runner with sequence slimming (`ISLAND_CLINKER`), a "Synteny (clinker)" panel in the page, and the NovInvenio_Investigations publish change.

**Tech Stack:** Python 3.12 (pixi), pytest, Biopython, Nextflow 26 DSL2 (strict parser), vanilla JS + `<canvas>` (the page must open from `file://`), node for JS unit checks, jsdom and headless Chromium for behaviour checks, gamcil/clinker 0.0.32 from PyPI, blastn (spot checks only).

**Spec:** `docs/superpowers/specs/2026-09-24-island-locus-view-design.md` (on `main` at `cc06785`, sections 1-8). The spec is the authority; this plan adds Rulings only where the spec is silent or open.

## Global Constraints

Values copied from the spec:

- Flank genes per side `F = 5`; fallback `F_min = 3` when no carrier strain has `F` genes on both sides; the locus then uses `F_min` flank columns and the page marks it "short flanks". If no strain meets `F_min`, use the strain with the most flank genes and mark the locus "exemplar at contig end".
- Exemplar: the strain that carries the locus's largest variant and has at least `F` genes on the same contig on both sides; ties: highest N50 (`assembly_quality_vs_content.tsv`), then strain name.
- Gene window `k = 10` (same as `pair_class_k`).
- Loci: an island joins the locus of a larger island when >= 50% of its families are in that island (the M4 criterion).
- Cell states: in place (a copy on contig X at rank r and a copy of at least one other column family on X within `k` ranks), elsewhere, rescue (only copy is a `genome_only` TBLASTN hit; drawn hatched), absent, contig break (absent or elsewhere, and the nearest in-place column copy is within `k` ranks of that contig's end). Positions are compared only within one strain's own assembly.
- Flanks intact: some left-flank and some right-flank family "in place" on the same contig, within `locus length + 2k` ranks of each other.
- Row classes: full locus (flanks intact, all locus columns in place); empty site (flanks intact, >= 80% of locus columns absent); partial (flanks intact, some present, some absent); uninformative (flanks not intact).
- Rows: strains collapsed into identical state vectors, grouped by row class, then by species.
- Breakpoint track: per boundary between adjacent columns, the number of flank-intact strains whose state changes between "in place" and "absent", split by species.
- Default ranking: at least 10 flank-intact strains in the empty-site class and at least 2 in the full-locus class, ranked by the smaller of the two counts. Alternative sorts: strain count, size, text search. `--top_loci` default 50.
- Payload: state vectors as 3-bit codes per column, stored once per collapsed pattern.
- Page: sidebar (size, variant count, empty-site / full / partial / uninformative counts, Pfam class chip); main (title = exemplar locus strain:contig:start-end; note with the exemplar rule and `F`/`k`; breakpoint track; column header with labels, class strip and anchor marks; grid; legend for the five states); popups (column: family, bin, domains, exemplar location; cell: strain(s), state, the strain's own contig:rank, and why the state was assigned); diagnostics banner kept; one line under it when `assembly_quality_confound` triggered.
- Clinker: gamcil/clinker **0.0.32 from PyPI**, added to the pixi environment as a PyPI dependency and to the container. Bioconda's `clinker` 1.33 is a different, unrelated RNA-seq tool; do not use it. `clinker <locus>/*.gbk -gf <locus>.groups.csv -p <locus>.html`; the tier-1 family is the group.
- Clinker strains: at most `--clinker_max_strains` (default 12; at most 66 pairs per locus): the exemplar; per row class (full, partial, empty) and species the flank-intact strain with the highest N50, ties by name; then more full or partial strains by N50 until the cap. Uninformative strains are never chosen.
- Clinker region per strain: on the contig where the flanks are intact, from the lowest to the highest rank among the strain's in-place column copies, then `F` genes further on each side (`F_min` for a short-flank locus), clipped at the contig end. For an empty-site strain this is the left-flank block plus the right-flank block.
- GenBank: one `.gbk` per (locus, strain); each CDS carries `/locus_tag` = the gene ID, `/translation` from the protein FASTA, and `/note="family=<tier-1 family>"`. `ISLAND_GBK_SLICE` is one task for all drawn loci.
- After clinker runs, remove `sequence` and `translation` from every gene object, in clusters and in links.
- Page panel "Synteny (clinker)" below the grid loads `clinker/<locus_id>.html` in an iframe, lists the strains shown, the reason each was chosen and each region; a locus without a file shows "no synteny figure for this locus"; with `--pangenome_clinker false` the panel says the step was not run. `--pangenome_clinker true|false`, default true.
- `clinker/` files are release assets (class 3), like `island_synteny.html`; `bin/sync_pangenome_report.py` (NII) copies `clinker/` with the page.

Repository rules that apply to every task (from `CLAUDE.md` and the user's global instructions):

- The page stays self-contained (no CDN, fetch or webfont), inserts untrusted strings with `textContent` only, and takes every colour from `lib/skins.py` tokens (`--series-1` = protein-search presence, `--series-2` = TBLASTN hit). No hex colour in a template fragment.
- No analysis output inside this repo. Every real-data command writes under `$SCRATCH` (node-local; use `${SCRATCH:?}`).
- `bin/` scripts: `#!/usr/bin/env python3`, `chmod +x`, argparse flags only. `lib/`: no I/O at import time. Read inputs through `lib/compressed_io.open_maybe_compressed()` (`.gz`/`.zst`).
- `pixi.toml` keeps `[workspace]`.
- Do not change `lib/Helpers.groovy` or the `process {}` block of `nextflow.config`: open PR branch `fix-194-189-project-name-presence-matrix-mem` changes the `projectName` fallback and `BUILD_PRESENCE_MATRIX` memory there. This plan only adds lines at the end of the `params {}` block and does not depend on that PR.
- Workflow per `CLAUDE.md` "Implementation workflow": one GitHub issue and one branch per part, commits reference the issue, one PR per part. Confirm with the user before creating an issue, pushing or opening a PR.
- Commands run from the worktree root. Tests: `pixi run python -m pytest -q <files>`.

## Review Focus

Inputs the spec implies but does not name, most likely first. Each has a test in the task that owns the code.

1. **A family fills two columns** (a tandem paralog in the exemplar's block). A strain with one lone copy of that family must not count as "in place" because the same copy sits in both columns. Test: `test_a_tandem_paralog_is_not_its_own_neighbour` (Task 4).
2. **Two strains of one locus share a protein ID** (for example NCBI-style `g1.t1` models). clinker's `-gf` file is keyed by label, so shared labels would mis-group genes. Expected: that locus's labels become `Short|protein_id`. Test: `test_protein_id_clash_across_strains_prefixes_labels` (Task 18).
3. **clinker is missing or fails for one locus** (the published container `0.5.0` has no clinker). Expected: the other loci still get pages, the failed one gets none, and the page says "No synteny figure for this locus." Tests: `tests/test_pangenome_island_clinker.py` (Task 19) and the jsdom `clinker panel:` checks (Task 21).
4. **Compressed study inputs**: `.fa.gz` genome and protein FASTAs, `.zst` position tables. Expected: read directly, no manual decompression. Tests: `test_gzipped_genome_and_proteins_are_read` (Task 18), `test_reads_zst_family_positions` (Task 7).
5. **No species data** (no `--config`). Expected: every strain in one "unknown species" group, one indel bar per boundary, no crash. Test: `test_without_config_all_strains_are_one_unknown_species_group` (Task 7).

## Rulings

Decisions on points the spec leaves open or silent. "Ruling Rn" in code comments refers to this list.

- **R1 Page.** `island_synteny.html` keeps its name, its header, tiles, skins and the whole island presence view. When the payload has `loci`, the page opens on the locus view; a two-button switch ("Loci (flank-anchored)" / "Islands (presence)") shows the old view. Without `loci` the page is unchanged. Reason: the user asked to keep the current layout as the base; the spec's open question 3 is answered "keep both, locus view default".
- **R2 Ranking.** The pipeline draws the top `--pangenome_top_loci` loci by `--pangenome_locus_rank` (`informative` default; `strains`, `size` also accepted). Loci that do not qualify for the informative rule (score -1) fill the remaining slots by carrier count. The sidebar re-sorts the drawn loci by informative polymorphism, strain count, size or locus ID, and has a text search.
- **R3 Carriers.** `islands_with_domains.tsv` does not list which strains carry an island. A strain carries the root (largest) variant when one contig holds a copy of every member family within `len(set(members)) - 1 + k` ranks.
- **R4 N50 gaps.** `assembly_quality_vs_content.tsv` covers the ingroup only (169 of 529 strains on the Coccidioides run). Strains with an N50 rank before strains without; without an N50, the tie-break is the gene count of the locus contig, then the name.
- **R5 Candidate pool.** States are computed for the top `--pangenome_locus_candidates` (200) loci by carrier proxy (sum of variant `n_strains`, capped at the strain count; by root size for `rank_by size`), then ranked. Measured: 1301 loci, 1031 with proxy >= 2.
- **R6 Present but unplaced.** A matrix call `present`/`genome_only` with no `family_positions` row is drawn "elsewhere"; its base state for row classes stays "absent" (the prototype rule). Measured: 0 such cells in the feasibility loci.
- **R7 Contig break.** Applied only to a column that lies outside the column range of the strain's in-place columns, when the nearest in-place copy is within `k` ranks of either end of its contig (orientation is Phase 2, so both ends count). `rescue, elsewhere` cells are not converted. A contig's ends are its minimum and maximum rank, because `rank` is a strain-wide running index (`bin/pangenome_build_family_positions.py`).
- **R8 Breakpoint track.** Indel bars count flank-intact strains whose state changes between in place (including rescue in place) and absent, one bar per species side by side in a fixed, labelled order, all in `--text-secondary` (position encodes species; hue does not). A last bar per boundary counts strains of any class whose state changes between in place and contig break, in `--warn` (the user's choice: contig breaks as a separate colour).
- **R9 Neighbour and class rules.** The in-place neighbour must be a copy of a different family, not only a different column (Review Focus 1). Row classes use the base states (rescue counts as its base state; contig break as absent). "Flanks intact" requires both flank copies to be in place, as the spec says; the regression test uses any copy pair, as the prototype did (the two rules differed on 1 of 26,450 strain-locus pairs on the top 50 loci by strain count).
- **R10 Clinker batching.** `ISLAND_CLINKER` runs `--pangenome_clinker_batch` (50) loci per task, not one task per locus. Reason: the user's HPCC rule sizes jobs toward about 1 hour; one locus measured 56 s at 4 cores, so 50 loci are one job of under an hour. Batch 1 restores the spec's one task per locus.
- **R11 Region anchors.** A strain's in-place copies count toward its region only when they sit on the flank-pair contig within `n_locus + 2k` ranks of the pair, so a stray in-place paralog elsewhere on the contig cannot stretch the region. The exemplar without a flank pair uses its carrier block.
- **R12 Locus keys.** Drawn loci get keys `L001`, `L002`, ... in drawn order; the clinker file is `clinker/<key>.html` (locus IDs contain `:`). The page only loads keys matching `^L[0-9]+$`.
- **R13 Popup detail.** Each collapsed row stores the per-column contig, rank and neighbour of its first strain, for the 300 largest rows per locus; the popup names that strain. The bp position appears only for the exemplar (existing `family_locations`). Measured: `island_loci.json` 1.94 MB, of which the state codes are about 83 KB; the page is about 2.7 MB. The spec's "well under 0.8 MB" estimate covered the codes only.
- **R14 GenBank.** Each gene is written as a `gene` plus `CDS` feature (silences clinker's "Could not find parent gene" warning: 145 warnings in the spike, 0 with the gene feature). TBLASTN rescue hits have no gene model and are not drawn. A strain whose DNA, protein or GFF3 file is missing is skipped with a warning. Rank rebuild uses the exact rule of `bin/pangenome_build_family_positions.py`; a mismatch with the region's anchor families stops the task.
- **R15 clinker failure.** A failed clinker run skips that locus with a warning; the page reports no figure for it.
- **R16 Slimming.** On by default; `--pangenome_clinker_slim false` keeps clinker's page unchanged.

## Facts measured while writing this plan (2026-09-26)

All on NII `studies/fungi/coccidioides_pangenome/results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome/` (529 strains: 169 *C. immitis* IN, 360 *C. posadasii* OUT). Planning code lived only in the session scratchpad; nothing was committed.

- The feasibility script reproduces every spec Feasibility number (both sort orders, `--empty_frac 0.8` and `1.0`), about 68 s and 90-180 MB each. This plan's library, in prototype mode, matches all of them exactly (Task 8).
- `family_positions.tsv.zst`: 5,074,953 rows, 236,315 (strain, contig) pairs; ranks within each contig are contiguous (0 gaps), so a rank difference equals a gene count.
- Locus grouping: 11,880 located islands form 1301 loci; 10,579 islands join a larger one (the spec's M4 number).
- `bin/pangenome_island_loci.py` on this run: 1 min 41 s wall (1 min 52 s on a second run), 1.36 GB peak RSS, `island_loci.json` 1.94 MB, 50 of 50 drawn loci have informative score >= 0; tiers 48 full, 1 short flanks, 1 contig end.
- Clinker chain on 3 loci x 12 strains: GenBank step 35 s, 60 MB; clinker 56 s and 128 MB for one locus at 4 cores; the page 5.35 MB, slimmed 1.01 MB; headless Chromium draws the same 237 genes, 12 clusters and 1178 link paths for both.
- `pixi` 0.71.3 installs clinker 0.0.32 from `[pypi-dependencies]`.
- Sequence spot check (Task 14's tool): **0 of 6 empty-site calls on L001-L003 confirmed**; the "empty" strains carry the locus DNA at 97.0-99.9% identity. At L001 the difference is one merged gene model versus two split models (see Task 14).

## File structure

| File | Part | Responsibility |
|---|---|---|
| `lib/island_locus.py` | A, B | Pure locus-view logic: loci, carriers, exemplar, columns, cell states, row classes, breakpoints, collapsed rows, ranking; Part B adds clinker strains and regions |
| `bin/pangenome_island_loci.py` | A, B | Streams inputs, runs the library, writes `island_loci.json` (and `island_regions.tsv`) |
| `lib/island_locus_template.py` | A, B | Locus view CSS/HTML/JS fragments; Part B adds the clinker panel |
| `lib/island_synteny_template.py` | A | Inserts the fragments; wraps the island view in `#island-view` |
| `bin/pangenome_island_synteny.py` | A, B | Embeds loci, confound flag, clinker metadata into the page payload |
| `modules/pangenome/island_loci.nf` | A, B | `ISLAND_LOCI` |
| `lib/genbank_slice.py` | B | Rank rebuild, GFF3 CDS parsing, GenBank records |
| `bin/pangenome_island_gbk_slice.py` | B | `ISLAND_GBK_SLICE` script |
| `lib/clinker_html.py` | B | Removes embedded sequences from a clinker page |
| `bin/pangenome_island_clinker.py` | B | Runs clinker per locus, slims, skips failures |
| `modules/pangenome/island_gbk_slice.nf`, `modules/pangenome/island_clinker.nf` | B | The two new processes |
| `workflows/pangenome_profile.nf`, `modules/pangenome/island_synteny.nf`, `nextflow.config`, `pangenome.nf` | A, B | Wiring, params, help, validation |
| `pixi.toml`, `pixi.lock`, `Dockerfile`, `novinvenio.def` | B | clinker 0.0.32 from PyPI |
| `tests/test_island_locus*.py`, `tests/test_pangenome_island_loci.py`, `tests/test_genbank_slice.py`, `tests/test_pangenome_island_gbk_slice.py`, `tests/test_clinker_html.py`, `tests/test_pangenome_island_clinker.py`, `tests/test_clinker_render.py` | A, B | New tests |
| `tests/test_pangenome_island_synteny.py`, `tests/test_report_js_behaviour.py`, `tests/js/drive_reports.mjs` | A, B | Extended tests |
| NII: `bin/sync_pangenome_report.py`, `lib/pangenome_site.py`, `.gitignore`, `bin/publish_report_release.sh`, `.github/workflows/static.yml`, `tests/test_sync_pangenome_report.py` | B | Publish `clinker/` (other repo) |

## Before Task 1 (setup for Part A)

```bash
git -C /bigdata/stajichlab/jstajich/projects/NovInvenio fetch -q origin
git -C /bigdata/stajichlab/jstajich/projects/NovInvenio worktree add \
  /bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/locus-view-a -b locus-view-a origin/main
cd /bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/locus-view-a
pixi install
```

Ask the user to approve one issue for Part A (`gh issue create --title "Island locus view (spec 2026-09-24 sections 1-7)" --body "See docs/superpowers/specs/2026-09-24-island-locus-view-design.md and docs/superpowers/plans/2026-09-26-island-locus-view-clinker.md, Part A."`). Then `export ISSUE_A=<the new issue number>` in the shell that runs the commits.

---

## Part A: locus computation and page (spec sections 1-7)

Part A ends in working software: `island_synteny.html` opens on the locus view for every pangenome run with `--pangenome_island_pfam_hmm`, and the island view is one click away.

### Task 1: Locus grouping (spec section 1)

**Files:**
- Create: `lib/island_locus.py` (module header, constants, section 1)
- Create: `tests/test_island_locus.py`

**Interfaces:**
- Consumes: rows of `report_tables/islands_with_domains.tsv` as `dict[str, str]` (keys `locus_id`, `member_families`, `n_strains`, `example_strain`, `locus_contig`).
- Produces: `Locus(root: dict, variants: list[dict])` with `.members -> list[str]`, `.locus_id -> str`, `.size -> int`, `.carrier_proxy(n_strains_total: int) -> int`; `group_loci(island_rows: list[dict], containment: float = 0.5) -> list[Locus]`; `candidate_loci(loci, n_strains_total: int, rank_by: str, n_candidates: int, min_strains: int = 2) -> list[Locus]`; `island_members(row) -> list[str]`; the state codes `ABSENT='0' IN_PLACE='1' ELSEWHERE='2' RESCUE_IN_PLACE='3' RESCUE_ELSEWHERE='4' CONTIG_BREAK='5'`, `IN_PLACE_CODES`, `STATE_LABELS`, `ROW_CLASSES = ('full', 'partial', 'empty', 'uninformative')`, and the `DEFAULT_*` / `INFORMATIVE_MIN_*` / `MAX_DETAIL_ROWS` constants used by later tasks.


- [ ] **Step 1: Write the failing test**

Create `tests/test_island_locus.py`:

```python
"""Unit tests for lib/island_locus.py on synthetic strains (spec validation
plan item 1). Positions are (contig, rank); a contig's ranks run from its
min to its max rank (see the module docstring). Each section imports the
names it adds, so the file grows task by task."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from island_locus import Locus, candidate_loci, group_loci  # noqa: E402


def isl(locus_id, members, n_strains=2, example="S1", contig="c1"):
    return {"locus_id": locus_id, "member_families": ",".join(members),
            "n_strains": str(n_strains), "example_strain": example, "locus_contig": contig}


# ---- group_loci -------------------------------------------------------------

def test_island_with_half_its_families_in_a_larger_island_joins_its_locus():
    loci = group_loci([isl("big", ["a", "b", "c", "d"]), isl("small", ["a", "x"])])
    assert [loc.locus_id for loc in loci] == ["big"]
    assert [v["locus_id"] for v in loci[0].variants] == ["big", "small"]


def test_island_below_containment_starts_its_own_locus():
    loci = group_loci([isl("big", ["a", "b", "c", "d"]), isl("small", ["a", "x", "y"])])
    assert sorted(loc.locus_id for loc in loci) == ["big", "small"]


def test_equal_size_islands_never_merge():
    loci = group_loci([isl("one", ["a", "b"]), isl("two", ["a", "b"])])
    assert len(loci) == 2


def test_grouping_is_transitive_to_the_root():
    rows = [isl("root", ["a", "b", "c", "d", "e", "f"]), isl("mid", ["a", "b", "c", "z"]),
            isl("leaf", ["z", "q"])]
    loci = group_loci(rows)
    assert [loc.locus_id for loc in loci] == ["root"]
    assert [v["locus_id"] for v in loci[0].variants] == ["root", "mid", "leaf"]


def test_tie_between_two_larger_islands_goes_to_the_earlier_one():
    rows = [isl("first", ["a", "b", "c"]), isl("second", ["x", "y", "z"]), isl("pair", ["a", "x"])]
    loci = group_loci(rows)
    by_id = {loc.locus_id: [v["locus_id"] for v in loc.variants] for loc in loci}
    assert by_id == {"first": ["first", "pair"], "second": ["second"]}


def test_unplaced_islands_are_not_grouped():
    assert group_loci([isl("-", ["a", "b"])]) == []


def test_carrier_proxy_is_capped_at_the_strain_count():
    loc = Locus(root={"members": ["a"], "n_strains": "8"},
                variants=[{"n_strains": "8"}, {"n_strains": "5"}])
    assert loc.carrier_proxy(10) == 10


def test_candidates_filter_on_min_strains_and_rank_by_proxy():
    loci = group_loci([isl("L1", ["a", "b"], 1), isl("L2", ["c", "d"], 9), isl("L3", ["e", "f", "g"], 3)])
    assert [loc.locus_id for loc in candidate_loci(loci, 20, "informative", 5, 2)] == ["L2", "L3"]
    assert [loc.locus_id for loc in candidate_loci(loci, 20, "size", 5, 2)] == ["L3", "L2"]
```


- [ ] **Step 2: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus.py
```

Expected: FAIL: `ModuleNotFoundError: No module named 'island_locus'`


- [ ] **Step 3: Make the change**

Create `lib/island_locus.py`:

```python
"""Island locus view: exemplar-anchored columns, per-strain cell states,
row classes and shared indel breakpoints.

Spec: docs/superpowers/specs/2026-09-24-island-locus-view-design.md.

Pure functions, no I/O at import time. A position is a (contig, rank) pair
from family_positions.tsv. A rank is a strain-wide running index
(bin/pangenome_build_family_positions.py sorts every strain's genes by
(contig, start) and enumerates them), so two ranks compare only on the same
contig, and a contig's first and last gene are its minimum and maximum rank.

"Ruling Rn" refers to the Rulings section of
docs/superpowers/plans/2026-09-26-island-locus-view-clinker.md.
"""
from __future__ import annotations

import collections
from dataclasses import dataclass

# ---- cell state codes (one character per column in a row's `codes`) -------
ABSENT = "0"
IN_PLACE = "1"
ELSEWHERE = "2"
RESCUE_IN_PLACE = "3"
RESCUE_ELSEWHERE = "4"
CONTIG_BREAK = "5"
STATE_LABELS = {
    ABSENT: "absent",
    IN_PLACE: "in place",
    ELSEWHERE: "elsewhere",
    RESCUE_IN_PLACE: "rescue (in place)",
    RESCUE_ELSEWHERE: "rescue (elsewhere)",
    CONTIG_BREAK: "contig break",
}
IN_PLACE_CODES = frozenset({IN_PLACE, RESCUE_IN_PLACE})

ROW_CLASSES = ("full", "partial", "empty", "uninformative")

DEFAULT_FLANK = 5
DEFAULT_FLANK_MIN = 3
DEFAULT_K = 10
DEFAULT_EMPTY_FRAC = 0.8
DEFAULT_CONTAINMENT = 0.5
INFORMATIVE_MIN_EMPTY = 10
INFORMATIVE_MIN_FULL = 2
MAX_DETAIL_ROWS = 300

Positions = dict  # {(strain, family): [(contig, rank), ...]}
ContigSpans = dict  # {(strain, contig): (min_rank, max_rank)}


def _as_int(value, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def island_members(row: dict) -> list[str]:
    return [f for f in row.get("member_families", "").split(",") if f]


# ---- 1. loci ---------------------------------------------------------------
@dataclass
class Locus:
    """A group of located islands: the largest (`root`) plus every island
    that joined it (spec section 1). `variants[0]` is the root."""
    root: dict
    variants: list

    @property
    def members(self) -> list[str]:
        return self.root["members"]

    @property
    def locus_id(self) -> str:
        return self.root.get("locus_id", "-")

    @property
    def size(self) -> int:
        return len(self.root["members"])

    def carrier_proxy(self, n_strains_total: int) -> int:
        """Sum of the variants' n_strains, capped at the strain count. A
        cheap pre-filter only: one strain can carry two variants."""
        return min(n_strains_total, sum(_as_int(v.get("n_strains")) for v in self.variants))


def group_loci(island_rows: list[dict], containment: float = DEFAULT_CONTAINMENT) -> list[Locus]:
    """Group located islands into loci (spec section 1, the M4 criterion).

    Islands are visited largest first (stable on input order, as the
    feasibility script does). An island joins the locus of a strictly larger,
    earlier island when >= `containment` of its families are in that island.
    With several such islands it joins the one sharing the most families;
    a tie goes to the earlier (larger) island. It inherits that island's
    locus, so grouping is transitive to the root.
    """
    rows = []
    for r in island_rows:
        if r.get("locus_id", "-") in ("-", "", None):
            continue
        members = island_members(r)
        if members:
            rows.append(dict(r, members=members))
    by_size = sorted(rows, key=lambda r: -len(r["members"]))
    sets = [set(r["members"]) for r in by_size]
    fam_index: dict[str, list[int]] = collections.defaultdict(list)
    locus_of: list[Locus] = []
    loci: list[Locus] = []
    for i, s in enumerate(sets):
        overlap = collections.Counter(
            j for f in s for j in fam_index[f] if len(sets[j]) > len(s))
        best = None
        for j, c in overlap.items():
            if c / len(s) < containment:
                continue
            if best is None or c > overlap[best] or (c == overlap[best] and j < best):
                best = j
        if best is None:
            loc = Locus(root=by_size[i], variants=[by_size[i]])
            loci.append(loc)
        else:
            loc = locus_of[best]
            loc.variants.append(by_size[i])
        locus_of.append(loc)
        for f in s:
            fam_index[f].append(i)
    return loci


def candidate_loci(loci: list[Locus], n_strains_total: int, rank_by: str,
                   n_candidates: int, min_strains: int = 2) -> list[Locus]:
    """Loci whose states are computed. `rank_by == "size"` pre-selects by root
    size; any other rank pre-selects by carrier_proxy()."""
    kept = [loc for loc in loci if loc.carrier_proxy(n_strains_total) >= min_strains]
    if rank_by == "size":
        kept.sort(key=lambda loc: (-loc.size, -loc.carrier_proxy(n_strains_total), loc.locus_id))
    else:
        kept.sort(key=lambda loc: (-loc.carrier_proxy(n_strains_total), -loc.size, loc.locus_id))
    return kept[:n_candidates]
```


- [ ] **Step 4: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus.py
```

Expected: 8 passed


- [ ] **Step 5: Commit**

```bash
git add lib/island_locus.py tests/test_island_locus.py
git commit -m "island locus view: group islands into loci (spec section 1) (#${ISSUE_A})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 2: Exemplar choice (spec section 2)

**Files:**
- Modify: `lib/island_locus.py` (append section 2)
- Modify: `tests/test_island_locus.py` (append)

**Interfaces:**
- Consumes: `positions: dict[(strain, family)] -> list[(contig, rank)]` (family_positions), `spans: dict[(strain, contig)] -> (min_rank, max_rank)`, `n50: dict[str, int]`.
- Produces: `Placement(strain, contig, lo, hi, left_avail, right_avail, contig_genes)` (frozen); `carrier_placements(members, positions, strains, spans, k=10) -> list[Placement]`; `quality_key(strain, n50, contig_genes) -> tuple`; `choose_exemplar(placements, n50, flank=5, flank_min=3) -> tuple[Placement, str] | None` with tier `'full'`, `'short_flanks'` or `'contig_end'`.


- [ ] **Step 1: Write the failing test**

Append to the end of `tests/test_island_locus.py`:

```python


from island_locus import Placement, carrier_placements, choose_exemplar, quality_key  # noqa: E402


# ---- exemplar -----------------------------------------------------------------

SPANS = {("S1", "c1"): (0, 19), ("S2", "c1"): (0, 19), ("S3", "c2"): (100, 104)}


def test_carrier_needs_every_member_within_the_span_window():
    pos = {("S1", "a"): [("c1", 8)], ("S1", "b"): [("c1", 9)],
           ("S2", "a"): [("c1", 1)], ("S2", "b"): [("c1", 19)]}
    places = carrier_placements(["a", "b"], pos, ["S1", "S2"], SPANS, k=10)
    assert [(p.strain, p.lo, p.hi, p.left_avail, p.right_avail) for p in places] == [("S1", 8, 9, 8, 10)]


def test_strain_missing_a_member_is_not_a_carrier():
    pos = {("S1", "a"): [("c1", 8)]}
    assert carrier_placements(["a", "b"], pos, ["S1"], SPANS) == []


def test_tightest_window_wins_for_a_multicopy_member():
    pos = {("S1", "a"): [("c1", 2), ("c1", 12)], ("S1", "b"): [("c1", 13)]}
    (p,) = carrier_placements(["a", "b"], pos, ["S1"], SPANS)
    assert (p.lo, p.hi) == (12, 13)


def place(strain, left, right, genes=20):
    return Placement(strain, "c1", left, left + 1, left, right, genes)


def test_exemplar_prefers_full_flanks_then_n50_then_name():
    picks = [place("B", 6, 6), place("A", 6, 6), place("C", 2, 9)]
    assert choose_exemplar(picks, {"A": 10, "B": 50}) == (picks[0], "full")
    assert choose_exemplar(picks, {"A": 50, "B": 50})[0].strain == "A"


def test_exemplar_falls_back_to_short_flanks():
    picks = [place("A", 4, 3), place("B", 9, 1)]
    assert choose_exemplar(picks, {}, flank=5, flank_min=3) == (picks[0], "short_flanks")


def test_exemplar_at_contig_end_takes_the_most_flank_genes():
    picks = [place("A", 0, 2), place("B", 1, 2)]
    assert choose_exemplar(picks, {}) == (picks[1], "contig_end")


def test_strains_without_n50_rank_after_known_n50_then_by_contig_genes():
    assert quality_key("X", {"Y": 1}, 500) > quality_key("Y", {"Y": 1}, 10)
    assert quality_key("X", {}, 500) < quality_key("W", {}, 100)


def test_no_carrier_gives_no_exemplar():
    assert choose_exemplar([], {}) is None
```


- [ ] **Step 2: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus.py
```

Expected: FAIL: `ImportError: cannot import name 'Placement' from 'island_locus'`


- [ ] **Step 3: Make the change**

Append to the end of `lib/island_locus.py`:

```python


# ---- 2. exemplar -----------------------------------------------------------
@dataclass(frozen=True)
class Placement:
    """Where one strain carries a locus's root variant."""
    strain: str
    contig: str
    lo: int
    hi: int
    left_avail: int
    right_avail: int
    contig_genes: int


def _min_cover_window(copies: list[tuple[int, str]], need: int) -> tuple[int, int] | None:
    """Smallest rank window over sorted (rank, family) copies that holds
    every one of `need` distinct families. Ties: lowest start."""
    best = None
    count: dict[str, int] = collections.Counter()
    have = 0
    left = 0
    for right, (r_rank, r_fam) in enumerate(copies):
        count[r_fam] += 1
        if count[r_fam] == 1:
            have += 1
        while have == need:
            lo, hi = copies[left][0], r_rank
            if best is None or hi - lo < best[1] - best[0]:
                best = (lo, hi)
            l_fam = copies[left][1]
            count[l_fam] -= 1
            if count[l_fam] == 0:
                have -= 1
            left += 1
    return best


def carrier_placements(members: list[str], positions: Positions, strains: list[str],
                       spans: ContigSpans, k: int = DEFAULT_K) -> list[Placement]:
    """Strains that carry the variant `members`, one Placement each.

    A strain carries it when one contig holds a copy of every member family
    within a span of at most len(set(members)) - 1 + k ranks (Ruling R3).
    Per strain the tightest window wins; ties go to the contig name.
    """
    fams = sorted(set(members))
    need = len(fams)
    out = []
    for s in strains:
        by_contig: dict[str, list[tuple[int, str]]] = collections.defaultdict(list)
        missing = False
        for f in fams:
            copies = positions.get((s, f))
            if not copies:
                missing = True
                break
            for contig, rank in copies:
                by_contig[contig].append((rank, f))
        if missing:
            continue
        best = None
        for contig in sorted(by_contig):
            copies = sorted(by_contig[contig])
            win = _min_cover_window(copies, need)
            if win is None or win[1] - win[0] > need - 1 + k:
                continue
            if best is None or win[1] - win[0] < best[2] - best[1]:
                best = (contig, win[0], win[1])
        if best is None:
            continue
        contig, lo, hi = best
        cmin, cmax = spans.get((s, contig), (lo, hi))
        out.append(Placement(s, contig, lo, hi, lo - cmin, cmax - hi, cmax - cmin + 1))
    return out


def quality_key(strain: str, n50: dict[str, int], contig_genes: int) -> tuple:
    """Sort key, best first. Strains with an N50 rank before strains without
    one; with an N50: N50 descending, then name (spec section 2). Without an
    N50 (assembly_quality_vs_content.tsv covers only the ingroup): gene count
    of the locus contig descending, then name (Ruling R4)."""
    if strain in n50:
        return (0, -n50[strain], strain)
    return (1, -contig_genes, strain)


def choose_exemplar(placements: list[Placement], n50: dict[str, int],
                    flank: int = DEFAULT_FLANK,
                    flank_min: int = DEFAULT_FLANK_MIN) -> tuple[Placement, str] | None:
    """(exemplar, tier). Tier "full": >= `flank` genes on both sides.
    "short_flanks": >= `flank_min` on both sides. "contig_end": the strain(s)
    with the most flank genes in total. None when no strain carries it."""
    if not placements:
        return None

    def best(pool):
        return min(pool, key=lambda p: quality_key(p.strain, n50, p.contig_genes))

    tier = [p for p in placements if p.left_avail >= flank and p.right_avail >= flank]
    if tier:
        return best(tier), "full"
    tier = [p for p in placements if p.left_avail >= flank_min and p.right_avail >= flank_min]
    if tier:
        return best(tier), "short_flanks"
    top = max(p.left_avail + p.right_avail for p in placements)
    return best([p for p in placements if p.left_avail + p.right_avail == top]), "contig_end"
```


- [ ] **Step 4: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus.py
```

Expected: 16 passed


- [ ] **Step 5: Commit**

```bash
git add lib/island_locus.py tests/test_island_locus.py
git commit -m "island locus view: carriers and exemplar choice with F_min fallback (spec section 2) (#${ISSUE_A})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 3: Columns from the exemplar's gene order (spec section 3)

**Files:**
- Modify: `lib/island_locus.py` (append section 3)
- Modify: `tests/test_island_locus.py` (append)

**Interfaces:**
- Consumes: the exemplar's contig order `list[(rank, family)]` and its `Placement` block.
- Produces: `Columns(left: tuple, locus: tuple, right: tuple, left_avail: int, right_avail: int)` with `.families -> list[str]`; `locus_columns(order, members, flank, block=None) -> Columns | None` (`block=None` is the feasibility script's rule, used by the regression in Task 8).


- [ ] **Step 1: Write the failing test**

Append to the end of `tests/test_island_locus.py`:

```python


from island_locus import locus_columns  # noqa: E402


# ---- columns ------------------------------------------------------------------

ORDER = [(r, f"g{r}") for r in range(20)]


def test_columns_take_the_block_and_flank_genes_each_side():
    cols = locus_columns(ORDER, ["g9", "g10"], 5, block=(9, 10))
    assert cols.left == ("g4", "g5", "g6", "g7", "g8")
    assert cols.locus == ("g9", "g10")
    assert cols.right == ("g11", "g12", "g13", "g14", "g15")
    assert (cols.left_avail, cols.right_avail) == (9, 9)


def test_columns_clip_at_the_contig_end():
    cols = locus_columns(ORDER, ["g1"], 5, block=(1, 1))
    assert cols.left == ("g0",) and len(cols.right) == 5


def test_block_holds_every_gene_between_its_ends_even_non_members():
    cols = locus_columns(ORDER, ["g9", "g11"], 1, block=(9, 11))
    assert cols.locus == ("g9", "g10", "g11")


def test_without_block_every_member_copy_on_the_contig_sets_the_block():
    cols = locus_columns(ORDER, ["g3", "g15"], 1)
    assert cols.locus[0] == "g3" and cols.locus[-1] == "g15"


def test_no_member_on_the_contig_gives_none():
    assert locus_columns(ORDER, ["zz"], 5) is None
```


- [ ] **Step 2: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus.py
```

Expected: FAIL: `ImportError: cannot import name 'locus_columns'`


- [ ] **Step 3: Make the change**

Append to the end of `lib/island_locus.py`:

```python


# ---- 3. columns ------------------------------------------------------------
@dataclass(frozen=True)
class Columns:
    left: tuple
    locus: tuple
    right: tuple
    left_avail: int
    right_avail: int

    @property
    def families(self) -> list[str]:
        return list(self.left) + list(self.locus) + list(self.right)


def locus_columns(order: list[tuple[int, str]], members, flank: int,
                  block: tuple[int, int] | None = None) -> Columns | None:
    """Columns from the exemplar's contig order [(rank, family), ...].

    `block` (lo, hi ranks) is the exemplar's Placement. Without it the block
    is every copy of a member family on the contig, which is the feasibility
    script's rule (used by the regression check). The locus block is every
    gene between the block ends, in rank order; the flanks are up to `flank`
    genes on each side. None when no member sits on the contig.
    """
    seq = sorted(order)
    mem = set(members)
    if block is None:
        idx = [i for i, (_, f) in enumerate(seq) if f in mem]
    else:
        idx = [i for i, (r, _) in enumerate(seq) if block[0] <= r <= block[1]]
    if not idx:
        return None
    lo, hi = min(idx), max(idx)
    return Columns(
        left=tuple(f for _, f in seq[max(0, lo - flank):lo]),
        locus=tuple(f for _, f in seq[lo:hi + 1]),
        right=tuple(f for _, f in seq[hi + 1:hi + 1 + flank]),
        left_avail=lo,
        right_avail=len(seq) - 1 - hi,
    )
```


- [ ] **Step 4: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus.py
```

Expected: 21 passed


- [ ] **Step 5: Commit**

```bash
git add lib/island_locus.py tests/test_island_locus.py
git commit -m "island locus view: flank + locus + flank columns (spec section 3) (#${ISSUE_A})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 4: Cell states (spec section 4)

**Files:**
- Modify: `lib/island_locus.py` (append section 4)
- Modify: `tests/test_island_locus.py` (append)

**Interfaces:**
- Consumes: `positions`, `spans` (Task 2 shapes), the column family list.
- Produces: `StrainCells(base: list[str], codes: list[str], detail: list, in_place_copies: list[(col, contig, rank)])`; `strain_cells(strain, columns, positions, spans, k=10, genome_only=frozenset(), present_unplaced=frozenset()) -> StrainCells`. `base` holds only `'0'/'1'/'2'` (the feasibility rule); `codes` adds rescue (`'3'`, `'4'`) and contig break (`'5'`); `detail[i]` is `None` or `(contig, rank, neighbour_col | None, delta | None)`.


- [ ] **Step 1: Write the failing test**

Append to the end of `tests/test_island_locus.py`:

```python


from island_locus import (  # noqa: E402
    ABSENT, CONTIG_BREAK, ELSEWHERE, IN_PLACE, RESCUE_ELSEWHERE, RESCUE_IN_PLACE, strain_cells,
)


# ---- cell states ------------------------------------------------------------------

COLS = ["L1", "A", "B", "R1"]


def test_in_place_needs_another_column_within_k_on_the_same_contig():
    pos = {("S1", "L1"): [("c1", 5)], ("S1", "A"): [("c1", 6)], ("S1", "B"): [("c2", 6)],
           ("S1", "R1"): [("c1", 30)]}
    sc = strain_cells("S1", COLS, pos, {}, k=10)
    assert sc.codes == [IN_PLACE, IN_PLACE, ELSEWHERE, ELSEWHERE]
    assert sc.detail[0] == ("c1", 5, 1, 1)
    assert sc.detail[2] == ("c2", 6, None, None)


def test_absent_when_no_copy():
    sc = strain_cells("S1", COLS, {("S1", "L1"): [("c1", 5)], ("S1", "R1"): [("c1", 7)]}, {})
    assert sc.codes == [IN_PLACE, ABSENT, ABSENT, IN_PLACE]


def test_a_genome_only_family_is_drawn_as_rescue():
    pos = {("S1", "L1"): [("c1", 5)], ("S1", "A"): [("c1", 6)], ("S1", "B"): [("c9", 1)]}
    sc = strain_cells("S1", COLS, pos, {}, genome_only=frozenset({"A", "B"}))
    assert sc.codes[:3] == [IN_PLACE, RESCUE_IN_PLACE, RESCUE_ELSEWHERE]
    assert sc.base[:3] == [IN_PLACE, IN_PLACE, ELSEWHERE]


def test_present_without_a_position_is_elsewhere_but_base_absent():
    sc = strain_cells("S1", COLS, {}, {}, present_unplaced=frozenset({"A"}))
    assert sc.codes[1] == ELSEWHERE and sc.base[1] == ABSENT


def test_contig_break_beyond_the_last_anchor_near_a_contig_end():
    pos = {("S1", "L1"): [("c1", 97)], ("S1", "A"): [("c1", 98)]}
    sc = strain_cells("S1", COLS, pos, {("S1", "c1"): (0, 99)}, k=10)
    assert sc.codes == [IN_PLACE, IN_PLACE, CONTIG_BREAK, CONTIG_BREAK]
    assert sc.base[2:] == [ABSENT, ABSENT]


def test_absence_between_anchors_is_never_a_contig_break():
    pos = {("S1", "L1"): [("c1", 97)], ("S1", "R1"): [("c1", 98)]}
    sc = strain_cells("S1", COLS, pos, {("S1", "c1"): (0, 99)})
    assert sc.codes == [IN_PLACE, ABSENT, ABSENT, IN_PLACE]


def test_absence_far_from_a_contig_end_is_absent():
    pos = {("S1", "L1"): [("c1", 50)], ("S1", "A"): [("c1", 51)]}
    sc = strain_cells("S1", COLS, pos, {("S1", "c1"): (0, 99)})
    assert sc.codes[2:] == [ABSENT, ABSENT]


def test_a_tandem_paralog_is_not_its_own_neighbour():
    # Review Focus 1: family P fills two columns (two copies in the exemplar).
    # A strain with one lone copy of P must not be "in place" because the
    # same copy sits in both columns.
    cols = ["L1", "P", "P", "R1"]
    sc = strain_cells("S1", cols, {("S1", "P"): [("c1", 40)]}, {})
    assert sc.codes == [ABSENT, ELSEWHERE, ELSEWHERE, ABSENT]
```


- [ ] **Step 2: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus.py
```

Expected: FAIL: `ImportError: cannot import name 'strain_cells'`


- [ ] **Step 3: Make the change**

Append to the end of `lib/island_locus.py`:

```python


# ---- 4. cell states ----------------------------------------------------------
@dataclass
class StrainCells:
    """One strain's view of one locus."""
    base: list  # ABSENT / IN_PLACE / ELSEWHERE only (the feasibility rule)
    codes: list  # display codes, adds rescue and contig break
    detail: list  # per column: None or (contig, rank, neighbour_col, delta)
    in_place_copies: list  # [(col, contig, rank)] for every in-place copy


def _nearest_other_family(copies: list[tuple[int, int]], idx: int, k: int,
                          columns: list[str]) -> tuple[int, int] | None:
    """(column, signed rank delta) of the nearest copy of a DIFFERENT FAMILY
    within `k` ranks of copies[idx], in a rank-sorted list of (rank, column).
    A family that fills two columns (a tandem paralog in the exemplar) is
    never its own neighbour. Ties: the left one."""
    rank, ci = copies[idx]
    fam = columns[ci]
    left = right = None
    j = idx - 1
    while j >= 0 and rank - copies[j][0] <= k:
        if columns[copies[j][1]] != fam:
            left = (copies[j][1], copies[j][0] - rank)
            break
        j -= 1
    j = idx + 1
    while j < len(copies) and copies[j][0] - rank <= k:
        if columns[copies[j][1]] != fam:
            right = (copies[j][1], copies[j][0] - rank)
            break
        j += 1
    if left is None:
        return right
    if right is None or -left[1] <= right[1]:
        return left
    return right


def strain_cells(strain: str, columns: list[str], positions: Positions,
                 spans: ContigSpans, k: int = DEFAULT_K,
                 genome_only: frozenset = frozenset(),
                 present_unplaced: frozenset = frozenset()) -> StrainCells:
    """Cell states for one strain (spec section 4).

    A copy is in place when a copy of a different column's family sits on
    the same contig within `k` ranks (Ruling R9). `genome_only` names the
    families this strain carries only as a TBLASTN rescue hit (matrix call
    genome_only).
    `present_unplaced` names families the matrix calls present with no
    position row; they are "elsewhere" (Ruling R6).
    """
    n = len(columns)
    by_contig: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    has_copy = [False] * n
    for ci, fam in enumerate(columns):
        for contig, rank in positions.get((strain, fam), ()):
            by_contig[contig].append((rank, ci))
            has_copy[ci] = True
    best_detail: list = [None] * n
    in_place_copies = []
    for contig, copies in by_contig.items():
        copies.sort()
        for idx, (rank, ci) in enumerate(copies):
            nbr = _nearest_other_family(copies, idx, k, columns)
            if nbr is None:
                continue
            in_place_copies.append((ci, contig, rank))
            cur = best_detail[ci]
            cand = (contig, rank, nbr[0], nbr[1])
            if cur is None or (abs(cand[3]), cand[0], cand[1]) < (abs(cur[3]), cur[0], cur[1]):
                best_detail[ci] = cand
    base = []
    codes = []
    detail = []
    for ci, fam in enumerate(columns):
        if best_detail[ci] is not None:
            b = IN_PLACE
            d = best_detail[ci]
        elif has_copy[ci]:
            b = ELSEWHERE
            contig, rank = min(positions[(strain, fam)])
            d = (contig, rank, None, None)
        elif fam in present_unplaced:
            b = ELSEWHERE
            d = None
        else:
            b = ABSENT
            d = None
        base.append(b if (has_copy[ci] or b == ABSENT) else ABSENT)
        rescued = fam in genome_only and b != ABSENT
        if rescued:
            codes.append(RESCUE_IN_PLACE if b == IN_PLACE else RESCUE_ELSEWHERE)
        else:
            codes.append(b)
        detail.append(d)
    _mark_contig_breaks(codes, best_detail, strain, spans, k)
    return StrainCells(base=base, codes=codes, detail=detail, in_place_copies=in_place_copies)


def _mark_contig_breaks(codes: list, anchors: list, strain: str, spans: ContigSpans, k: int) -> None:
    """Turn ABSENT/ELSEWHERE into CONTIG_BREAK where the locus runs off the
    assembly (spec section 4, Ruling R7): the column lies outside the column
    range of the strain's in-place columns, and the nearest in-place copy is
    within `k` ranks of either end of its contig."""
    placed = [ci for ci, a in enumerate(anchors) if a is not None]
    if not placed:
        return
    first, last = placed[0], placed[-1]
    for ci, code in enumerate(codes):
        if code not in (ABSENT, ELSEWHERE) or first <= ci <= last:
            continue
        j = first if ci < first else last
        contig, rank = anchors[j][0], anchors[j][1]
        cmin, cmax = spans.get((strain, contig), (None, None))
        if cmin is None:
            continue
        if rank - cmin <= k or cmax - rank <= k:
            codes[ci] = CONTIG_BREAK
```


- [ ] **Step 4: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus.py
```

Expected: 29 passed


- [ ] **Step 5: Commit**

```bash
git add lib/island_locus.py tests/test_island_locus.py
git commit -m "island locus view: in place / elsewhere / rescue / absent / contig break (spec section 4) (#${ISSUE_A})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 5: Flanks intact, row classes, breakpoint track, row collapse, ranking (spec sections 5-6)

**Files:**
- Modify: `lib/island_locus.py` (append section 5)
- Modify: `tests/test_island_locus.py` (append)

**Interfaces:**
- Consumes: `StrainCells.base` / `.codes` from Task 4.
- Produces: `flank_pair(positions, strain, left, right, n_locus, k=10, in_place=None) -> (contig, lo, hi) | None`; `row_class(base, n_left, n_locus, intact, empty_frac=0.8) -> str`; `breakpoint_track(strain_rows: list[(species, row_class, codes)], n_cols) -> list[{'b', 'indel': {species: n}, 'contig_break': n}]`; `collapse_rows(per_strain: {strain: (row_class, codes)}, details=None, max_detail_rows=300) -> list[{'row_class', 'codes', 'count', 'strains', 'rep'}]`; `encode_detail(detail) -> {'c': [contig], 'd': [cell]}`; `informative_score(counts) -> int`; `rank_key(result, rank_by) -> tuple`.


- [ ] **Step 1: Write the failing test**

Append to the end of `tests/test_island_locus.py`:

```python


from island_locus import (  # noqa: E402
    breakpoint_track, collapse_rows, flank_pair, informative_score, rank_key, row_class,
)


# ---- flanks, row classes, breakpoints ------------------------------------------------

def test_flank_pair_on_one_contig_within_locus_plus_2k():
    pos = {("S1", "L1"): [("c1", 5)], ("S1", "R1"): [("c1", 27)]}
    assert flank_pair(pos, "S1", ["L1"], ["R1"], 2, k=10) == ("c1", 5, 27)
    assert flank_pair(pos, "S1", ["L1"], ["R1"], 1, k=10) is None


def test_flank_pair_on_different_contigs_is_not_intact():
    pos = {("S1", "L1"): [("c1", 5)], ("S1", "R1"): [("c2", 6)]}
    assert flank_pair(pos, "S1", ["L1"], ["R1"], 2) is None


def test_flank_pair_can_require_in_place_copies():
    pos = {("S1", "L1"): [("c1", 5)], ("S1", "R1"): [("c1", 6)]}
    assert flank_pair(pos, "S1", ["L1"], ["R1"], 2, in_place={("c1", 5)}) is None
    assert flank_pair(pos, "S1", ["L1"], ["R1"], 2, in_place={("c1", 5), ("c1", 6)}) == ("c1", 5, 6)


def test_row_classes():
    full = [IN_PLACE] * 4
    assert row_class(full, 1, 2, True) == "full"
    assert row_class([IN_PLACE, IN_PLACE, ABSENT, IN_PLACE], 1, 2, True) == "partial"
    assert row_class([IN_PLACE, ABSENT, ABSENT, IN_PLACE], 1, 2, True) == "empty"
    assert row_class(full, 1, 2, False) == "uninformative"


def test_empty_site_threshold_is_80_percent_of_locus_columns():
    four_of_five = [IN_PLACE] + [ABSENT] * 4 + [IN_PLACE] + [IN_PLACE]
    three_of_five = [IN_PLACE] + [ABSENT] * 3 + [IN_PLACE] * 3
    assert row_class(four_of_five, 1, 5, True, 0.8) == "empty"
    assert row_class(three_of_five, 1, 5, True, 0.8) == "partial"


def test_elsewhere_block_is_partial_not_full():
    assert row_class([IN_PLACE, ELSEWHERE, IN_PLACE], 1, 1, True) == "partial"


def test_breakpoint_track_counts_flank_intact_indels_by_species_and_contig_breaks():
    rows = [("sp1", "empty", "1001"), ("sp1", "full", "1111"), ("sp2", "partial", "1101"),
            ("sp2", "uninformative", "1001"), ("sp1", "uninformative", "1155")]
    track = breakpoint_track(rows, 4)
    assert track == [
        {"b": 1, "indel": {"sp1": 1}, "contig_break": 0},
        {"b": 2, "indel": {"sp2": 1}, "contig_break": 1},
        {"b": 3, "indel": {"sp1": 1, "sp2": 1}, "contig_break": 0},
    ]


def test_collapse_groups_same_class_and_codes_and_orders_by_class():
    rows = collapse_rows({"a": ("empty", "100"), "b": ("full", "111"), "c": ("full", "111"),
                          "d": ("uninformative", "111")},
                         details={"b": [("c1", 1, 1, 1), None, ("c1", 3, None, None)]})
    assert [(r["row_class"], r["count"], r["strains"]) for r in rows] == [
        ("full", 2, ["b", "c"]), ("empty", 1, ["a"]), ("uninformative", 1, ["d"])]
    assert rows[0]["rep"] == {"c": ["c1"], "d": [[0, 1, 1, 1], None, [0, 3]]}


def test_detail_is_kept_only_for_the_largest_rows():
    per = {f"s{i}": ("full", format(i, "03b")) for i in range(4)}
    per["s0b"] = ("full", "000")
    rows = collapse_rows(per, details={s: [None, None, None] for s in per}, max_detail_rows=1)
    assert [r["rep"] is not None for r in rows] == [True, False, False, False]


def test_informative_score_needs_10_empty_and_2_full():
    assert informative_score({"empty": 10, "full": 2}) == 2
    assert informative_score({"empty": 9, "full": 50}) == -1
    assert informative_score({"empty": 30, "full": 1}) == -1


def test_rank_key_orders():
    a = {"informative_score": 5, "n_carriers": 10, "size": 2, "locus_id": "a"}
    b = {"informative_score": -1, "n_carriers": 90, "size": 9, "locus_id": "b"}
    assert sorted([b, a], key=lambda r: rank_key(r, "informative"))[0] is a
    assert sorted([a, b], key=lambda r: rank_key(r, "strains"))[0] is b
    assert sorted([a, b], key=lambda r: rank_key(r, "size"))[0] is b
```


- [ ] **Step 2: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus.py
```

Expected: FAIL: `ImportError: cannot import name 'breakpoint_track'`


- [ ] **Step 3: Make the change**

Append to the end of `lib/island_locus.py`:

```python


# ---- 5. row classes and breakpoints -----------------------------------------
def flank_pair(positions: Positions, strain: str, left: list[str], right: list[str],
               n_locus: int, k: int = DEFAULT_K,
               in_place: set | None = None) -> tuple[str, int, int] | None:
    """(contig, lower rank, upper rank) of the closest left-flank/right-flank
    copy pair on one contig within n_locus + 2k ranks, or None (spec section
    5, "flanks intact"). With `in_place` (a set of (contig, rank) in-place
    copies) both copies must be in place; without it any copies count, which
    is the feasibility script's rule."""
    span = n_locus + 2 * k
    best = None
    for f in left:
        for c1, r1 in positions.get((strain, f), ()):
            if in_place is not None and (c1, r1) not in in_place:
                continue
            for g in right:
                for c2, r2 in positions.get((strain, g), ()):
                    if c2 != c1 or abs(r1 - r2) > span:
                        continue
                    if in_place is not None and (c2, r2) not in in_place:
                        continue
                    cand = (abs(r1 - r2), c1, min(r1, r2), max(r1, r2))
                    if best is None or cand < best:
                        best = cand
    if best is None:
        return None
    return best[1], best[2], best[3]


def row_class(base: list, n_left: int, n_locus: int, intact: bool,
              empty_frac: float = DEFAULT_EMPTY_FRAC) -> str:
    """full / partial / empty / uninformative (spec section 5), from the
    base states (rescue counts as its base state; contig break as absent)."""
    if not intact:
        return "uninformative"
    block = base[n_left:n_left + n_locus]
    if block and sum(b == ABSENT for b in block) / len(block) >= empty_frac:
        return "empty"
    if all(b == IN_PLACE for b in block):
        return "full"
    return "partial"


def breakpoint_track(strain_rows: list[tuple[str, str, str]], n_cols: int) -> list[dict]:
    """Per boundary b (between columns b-1 and b): flank-intact strains whose
    state changes between in place and absent, by species; and strains of any
    class whose state changes between in place and contig break (Ruling R8).
    `strain_rows` = [(species, row_class, codes)]. Only non-zero boundaries."""
    out = []
    for b in range(1, n_cols):
        indel: dict[str, int] = collections.Counter()
        brk = 0
        for species, cls, codes in strain_rows:
            pair = {codes[b - 1], codes[b]}
            placed = bool(pair & IN_PLACE_CODES)
            if placed and ABSENT in pair and cls != "uninformative":
                indel[species] += 1
            if placed and CONTIG_BREAK in pair:
                brk += 1
        if indel or brk:
            out.append({"b": b, "indel": dict(sorted(indel.items())), "contig_break": brk})
    return out


def collapse_rows(per_strain: dict[str, tuple[str, str]],
                  details: dict[str, list] | None = None,
                  max_detail_rows: int = MAX_DETAIL_ROWS) -> list[dict]:
    """Collapse strains with the same (row class, codes) into one row.
    Sorted by ROW_CLASSES order, then codes. `rep` is the per-column detail
    of the row's first strain, kept for the `max_detail_rows` rows with the
    most strains (None beyond that, to bound the payload); see
    encode_detail() for its shape."""
    groups: dict[tuple[str, str], list[str]] = collections.defaultdict(list)
    for strain, key in per_strain.items():
        groups[key].append(strain)
    rows = []
    for (cls, codes), strains in groups.items():
        strains.sort()
        rows.append({"row_class": cls, "codes": codes, "count": len(strains), "strains": strains,
                     "rep": None})
    rows.sort(key=lambda r: (ROW_CLASSES.index(r["row_class"]), r["codes"]))
    if details is not None:
        by_count = sorted(range(len(rows)), key=lambda i: (-rows[i]["count"], i))
        for i in by_count[:max_detail_rows]:
            rep = details.get(rows[i]["strains"][0])
            rows[i]["rep"] = encode_detail(rep) if rep else None
    return rows


def encode_detail(detail: list) -> dict:
    """Compact per-column detail for one strain: {"c": [contig names],
    "d": [cell per column]}. A cell is null (no copy), [contig index, rank]
    (a copy with no neighbour) or [contig index, rank, neighbour column,
    signed rank delta] (in place)."""
    contigs: list[str] = []
    index: dict[str, int] = {}
    cells = []
    for d in detail:
        if d is None:
            cells.append(None)
            continue
        contig, rank, nbr, delta = d
        if contig not in index:
            index[contig] = len(contigs)
            contigs.append(contig)
        cells.append([index[contig], rank] if nbr is None else [index[contig], rank, nbr, delta])
    return {"c": contigs, "d": cells}


def informative_score(counts: dict[str, int]) -> int:
    """Spec section 6: min(empty, full) when >= 10 empty-site and >= 2
    full-locus strains, else -1."""
    empty, full = counts.get("empty", 0), counts.get("full", 0)
    if empty >= INFORMATIVE_MIN_EMPTY and full >= INFORMATIVE_MIN_FULL:
        return min(empty, full)
    return -1


def rank_key(result: dict, rank_by: str) -> tuple:
    """Sort key over computed locus results (dicts with informative_score,
    n_carriers, size, locus_id); best first."""
    if rank_by == "size":
        return (-result["size"], -result["n_carriers"], result["locus_id"])
    if rank_by == "strains":
        return (-result["n_carriers"], -result["size"], result["locus_id"])
    return (-result["informative_score"], -result["n_carriers"], -result["size"], result["locus_id"])
```


- [ ] **Step 4: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus.py
```

Expected: 40 passed


- [ ] **Step 5: Commit**

```bash
git add lib/island_locus.py tests/test_island_locus.py
git commit -m "island locus view: row classes, breakpoints, collapsed rows, ranking (spec sections 5-6) (#${ISSUE_A})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 6: One locus across all strains, and its page entry

**Files:**
- Modify: `lib/island_locus.py` (append section 6)
- Modify: `tests/test_island_locus.py` (append)

**Interfaces:**
- Consumes: everything from Tasks 1-5; `matrix` with `.call(family, strain) -> 'present' | 'genome_only' | 'absent'` (`lib/pangenome_matrix.PresenceMatrix`).
- Produces: `compute_locus(locus, placement, tier, columns, strains, positions, spans, matrix, species_of, k=10, empty_frac=0.8) -> dict` with the payload keys `locus_id size n_variants variant_strains exemplar exemplar_contig tier families n_left n_locus n_right counts counts_by_species n_carriers informative_score rows breakpoints` plus private `_cells`, `_pairs`, `_row_class`; `locus_payload(result, key, bins, family_classes, family_domains, dominant, family_locations=None, exemplar_span=None) -> dict` (drops `_` keys, adds `key family_bins family_classes family_domains dominant_class exemplar_span [family_locations]`).


- [ ] **Step 1: Write the failing test**

Append to the end of `tests/test_island_locus.py`:

```python


from island_locus import Columns, compute_locus  # noqa: E402


# ---- compute_locus end to end --------------------------------------------------------

class FakeMatrix:
    def __init__(self, calls):
        self.calls = calls

    def call(self, fam, strain):
        return self.calls.get((fam, strain), "absent")


def test_compute_locus_classifies_four_synthetic_strains():
    cols = Columns(left=("L1",), locus=("A", "B"), right=("R1",), left_avail=5, right_avail=5)
    pos = {}
    calls = {}

    def put(strain, contig, ranks):
        for fam, r in zip(["L1", "A", "B", "R1"], ranks):
            if r is not None:
                pos[(strain, fam)] = [(contig, r)]
                calls[(fam, strain)] = "present"

    put("FULL", "c1", [10, 11, 12, 13])
    put("EMPTY", "c1", [10, None, None, 11])
    put("PART", "c1", [10, 11, None, 13])
    put("FRAG", "c2", [0, 1, None, None])
    spans = {("FULL", "c1"): (0, 50), ("EMPTY", "c1"): (0, 50), ("PART", "c1"): (0, 50),
             ("FRAG", "c2"): (0, 2)}
    loc = Locus(root={"members": ["A", "B"], "locus_id": "FULL:c1:1-9"},
                variants=[{"n_strains": "2"}])
    res = compute_locus(loc, Placement("FULL", "c1", 11, 12, 11, 38, 51), "full", cols,
                        ["EMPTY", "FRAG", "FULL", "PART"], pos, spans, FakeMatrix(calls),
                        {"FULL": "sp1", "EMPTY": "sp2", "PART": "sp1"})
    assert res["counts"] == {"full": 1, "partial": 1, "empty": 1, "uninformative": 1}
    assert res["_row_class"] == {"EMPTY": "empty", "FRAG": "uninformative", "FULL": "full",
                                 "PART": "partial"}
    frag = [r for r in res["rows"] if r["strains"] == ["FRAG"]][0]
    assert frag["codes"] == "1155"
    assert res["n_carriers"] == 3
    assert res["counts_by_species"]["sp2"]["empty"] == 1
    assert res["breakpoints"] == [
        {"b": 1, "indel": {"sp2": 1}, "contig_break": 0},
        {"b": 2, "indel": {"sp1": 1}, "contig_break": 1},
        {"b": 3, "indel": {"sp1": 1, "sp2": 1}, "contig_break": 0},
    ]
    assert res["informative_score"] == -1
```


- [ ] **Step 2: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus.py
```

Expected: FAIL: `ImportError: cannot import name 'compute_locus'`


- [ ] **Step 3: Make the change**

Append to the end of `lib/island_locus.py`:

```python


# ---- 6. one locus, all strains -------------------------------------------------
def compute_locus(locus: Locus, placement: Placement, tier: str, columns: Columns,
                  strains: list[str], positions: Positions, spans: ContigSpans,
                  matrix, species_of: dict[str, str], k: int = DEFAULT_K,
                  empty_frac: float = DEFAULT_EMPTY_FRAC) -> dict:
    """Every strain's cells, row class and flank pair for one locus, plus the
    collapsed rows, counts and breakpoint track. `matrix` has
    .call(family, strain) -> "present" | "genome_only" | "absent" (a
    lib.pangenome_matrix.PresenceMatrix). Keys starting with "_" are for
    the caller (Part B regions) and never go into the page payload."""
    fams = columns.families
    n_left, n_locus = len(columns.left), len(columns.locus)
    per_strain: dict[str, tuple[str, str]] = {}
    details: dict[str, list] = {}
    cells: dict[str, StrainCells] = {}
    pairs: dict[str, tuple[str, int, int]] = {}
    counts = {c: 0 for c in ROW_CLASSES}
    by_species: dict[str, dict[str, int]] = {}
    track_rows = []
    n_carriers = 0
    for s in strains:
        calls = {f: matrix.call(f, s) for f in set(fams)}
        genome_only = frozenset(f for f, c in calls.items() if c == "genome_only")
        unplaced = frozenset(f for f, c in calls.items()
                             if c != "absent" and not positions.get((s, f)))
        sc = strain_cells(s, fams, positions, spans, k, genome_only, unplaced)
        in_place = {(c, r) for _, c, r in sc.in_place_copies}
        pair = flank_pair(positions, s, list(columns.left), list(columns.right), n_locus, k,
                          in_place=in_place)
        cls = row_class(sc.base, n_left, n_locus, pair is not None, empty_frac)
        codes = "".join(sc.codes)
        per_strain[s] = (cls, codes)
        details[s] = sc.detail
        cells[s] = sc
        if pair is not None:
            pairs[s] = pair
        counts[cls] += 1
        sp = species_of.get(s, "")
        by_species.setdefault(sp, {c: 0 for c in ROW_CLASSES})[cls] += 1
        track_rows.append((sp, cls, codes))
        if any(b == IN_PLACE for b in sc.base[n_left:n_left + n_locus]):
            n_carriers += 1
    return {
        "locus_id": locus.locus_id,
        "size": locus.size,
        "n_variants": len(locus.variants),
        "variant_strains": locus.carrier_proxy(len(strains)),
        "exemplar": placement.strain,
        "exemplar_contig": placement.contig,
        "tier": tier,
        "families": fams,
        "n_left": n_left,
        "n_locus": n_locus,
        "n_right": len(columns.right),
        "counts": counts,
        "counts_by_species": dict(sorted(by_species.items())),
        "n_carriers": n_carriers,
        "informative_score": informative_score(counts),
        "rows": collapse_rows(per_strain, details),
        "breakpoints": breakpoint_track(track_rows, len(fams)),
        "_cells": cells,
        "_pairs": pairs,
        "_row_class": {s: v[0] for s, v in per_strain.items()},
    }


def locus_payload(result: dict, key: str, bins: dict[str, str],
                  family_classes: list[str], family_domains: list[str],
                  dominant: str, family_locations: list | None = None,
                  exemplar_span: tuple[int, int] | None = None) -> dict:
    """The JSON-safe page entry for one computed locus (drops "_" keys)."""
    out = {k: v for k, v in result.items() if not k.startswith("_")}
    out["key"] = key
    out["family_bins"] = [bins.get(f, "") for f in result["families"]]
    out["family_classes"] = family_classes
    out["family_domains"] = family_domains
    out["dominant_class"] = dominant
    out["exemplar_span"] = ({"start": exemplar_span[0], "end": exemplar_span[1]}
                            if exemplar_span else None)
    if family_locations is not None:
        out["family_locations"] = family_locations
    return out
```


- [ ] **Step 4: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus.py
```

Expected: 41 passed


- [ ] **Step 5: Commit**

```bash
git add lib/island_locus.py tests/test_island_locus.py
git commit -m "island locus view: compute one locus over all strains (#${ISSUE_A})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 7: CLI `bin/pangenome_island_loci.py` (three filtered passes, `island_loci.json`)

**Files:**
- Create: `bin/pangenome_island_loci.py`
- Create: `tests/test_pangenome_island_loci.py`

**Interfaces:**
- Consumes: Tasks 1-6; `bin/pangenome_island_synteny.py::load_gene_locations`, `lib/island_synteny.py::column_locations`, `bin/pangenome_domain_enrichment.py::parse_domtblout`.
- Produces: `scan_family_positions(path, families=None, contigs=None, spans=False) -> PositionScan(positions, orders, spans)`; `read_matrix_strains(path) -> list[str]`; `read_n50(path) -> dict[str, int]`; `read_bins(path) -> dict[str, str]`; `build(args) -> dict`; `main(argv=None) -> int`. Output `island_loci.json` = `{project, locus_params, n_loci_total, n_loci_candidates, n_loci_unplaced, n_strains_with_n50, n_strains, loci: [page entries]}`; locus keys `L001`, `L002`, ... in drawn order (Ruling R12).


- [ ] **Step 1: Write the failing test**

Create `tests/test_pangenome_island_loci.py`:

```python
"""CLI tests for bin/pangenome_island_loci.py on a four-strain fixture."""
import json
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).parent.parent / "bin"
sys.path.insert(0, str(BIN))
from pangenome_island_loci import main, read_n50, scan_family_positions  # noqa: E402

LEFT = ["F1", "F2", "F3", "F4", "F5"]
RIGHT = ["G1", "G2", "G3", "G4", "G5"]


def write_fixture(d: Path) -> dict:
    (d / "islands.tsv").write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\tlocus_id\tlocus_contig\t"
        "locus_start\tlocus_end\tpfam_domains\n"
        "2\tS1\t2\tA,B\tS1:c1:100-200\tc1\t100\t200\t-\n"
        "1\tS2\t1\tA\tS2:c1:100-150\tc1\t100\t150\t-\n")
    rows = ["Short\tfamily\tcontig\trank"]
    full = LEFT + ["A", "B"] + RIGHT
    for s in ("S1", "S2"):
        rows += [f"{s}\t{f}\tc1\t{r}" for r, f in enumerate(full)]
    rows += [f"S3\t{f}\tc1\t{r}" for r, f in enumerate(LEFT + RIGHT)]
    rows += ["S4\tF1\tc9\t0"]
    (d / "family_positions.tsv").write_text("\n".join(rows) + "\n")
    fams = LEFT + ["A", "B"] + RIGHT
    calls = {"S1": "present", "S2": "present", "S3": "present", "S4": "absent"}
    lines = ["family\tS1\tS2\tS3\tS4"]
    for f in fams:
        cells = []
        for s in ("S1", "S2", "S3", "S4"):
            if f in ("A", "B") and s == "S3":
                cells.append("absent")
            elif s == "S4":
                cells.append("present" if f == "F1" else "absent")
            else:
                cells.append(calls[s])
        lines.append(f + "\t" + "\t".join(cells))
    (d / "matrix.tsv").write_text("\n".join(lines) + "\n")
    (d / "aq.tsv").write_text("Short\tn_contigs\tn50\nS1\t10\t500\nS2\t10\t900\n")
    (d / "config.csv").write_text(
        "GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup\n"
        "IN,Sp one,S1,S1.pep.fa,S1.dna.fa,S1.gff3,S1,t\nIN,Sp one,S2,S2.pep.fa,S2.dna.fa,S2.gff3,S2,t\n"
        "OUT,Sp two,S3,S3.pep.fa,S3.dna.fa,S3.gff3,S3,t\nOUT,Sp two,S4,S4.pep.fa,S4.dna.fa,S4.gff3,S4,t\n")
    (d / "freq.tsv").write_text("family\tfrequency\tstrain_count\tbin\n" +
                                "".join(f"{f}\t1\t4\tcore\n" for f in LEFT + RIGHT) +
                                "A\t0.5\t2\tshell\nB\t0.5\t2\tshell\n")
    return {"islands": d / "islands.tsv", "fp": d / "family_positions.tsv", "matrix": d / "matrix.tsv"}


def run(d: Path, *extra) -> dict:
    write_fixture(d)
    out = d / "island_loci.json"
    rc = main(["--islands_with_domains", str(d / "islands.tsv"),
               "--presence_matrix", str(d / "matrix.tsv"),
               "--family_positions", str(d / "family_positions.tsv"),
               "--frequency_table", str(d / "freq.tsv"),
               "--assembly_quality", str(d / "aq.tsv"),
               "--config", str(d / "config.csv"),
               "--project", "demo", "--output", str(out), *extra])
    assert rc == 0
    return json.loads(out.read_text())


def test_one_locus_with_both_islands_as_variants(tmp_path):
    data = run(tmp_path)
    assert data["n_loci_total"] == 1
    (locus,) = data["loci"]
    assert locus["key"] == "L001"
    assert locus["n_variants"] == 2
    assert locus["locus_id"] == "S1:c1:100-200"


def test_exemplar_is_the_higher_n50_carrier(tmp_path):
    (locus,) = run(tmp_path)["loci"]
    assert locus["exemplar"] == "S2"
    assert locus["tier"] == "full"
    assert locus["families"] == LEFT + ["A", "B"] + RIGHT
    assert (locus["n_left"], locus["n_locus"], locus["n_right"]) == (5, 2, 5)


def test_row_classes_and_counts(tmp_path):
    (locus,) = run(tmp_path)["loci"]
    assert locus["counts"] == {"full": 2, "partial": 0, "empty": 1, "uninformative": 1}
    by_class = {r["row_class"]: r for r in locus["rows"]}
    assert by_class["empty"]["strains"] == ["S3"]
    assert by_class["empty"]["codes"] == "11111" + "00" + "11111"
    assert by_class["full"]["count"] == 2


def test_bins_species_and_params_reach_the_payload(tmp_path):
    data = run(tmp_path)
    (locus,) = data["loci"]
    assert locus["family_bins"][5] == "shell"
    assert locus["counts_by_species"]["Sp two"]["empty"] == 1
    assert data["locus_params"]["flank"] == 5 and data["locus_params"]["k"] == 10
    assert data["n_strains_with_n50"] == 2


def test_top_loci_zero_draws_nothing_and_still_writes_json(tmp_path):
    data = run(tmp_path, "--top_loci", "0")
    assert data["loci"] == [] and data["n_loci_total"] == 1


def test_no_islands_gives_an_empty_valid_payload(tmp_path):
    write_fixture(tmp_path)
    (tmp_path / "islands.tsv").write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\tlocus_id\tlocus_contig\n")
    out = tmp_path / "o.json"
    assert main(["--islands_with_domains", str(tmp_path / "islands.tsv"),
                 "--presence_matrix", str(tmp_path / "matrix.tsv"),
                 "--family_positions", str(tmp_path / "family_positions.tsv"),
                 "--project", "demo", "--output", str(out)]) == 0
    assert json.loads(out.read_text())["loci"] == []


def test_scan_keeps_spans_orders_and_family_copies(tmp_path):
    write_fixture(tmp_path)
    scan = scan_family_positions(str(tmp_path / "family_positions.tsv"), families={"A"},
                                 contigs={("S1", "c1")}, spans=True)
    assert scan.positions[("S1", "A")] == [("c1", 5)]
    assert scan.orders[("S1", "c1")][0] == (0, "F1")
    assert scan.spans[("S3", "c1")] == (0, 9)


def test_reads_zst_family_positions(tmp_path):
    write_fixture(tmp_path)
    subprocess.run(["zstd", "-q", "-f", str(tmp_path / "family_positions.tsv")], check=True)
    scan = scan_family_positions(str(tmp_path / "family_positions.tsv.zst"), families={"B"})
    assert scan.positions[("S2", "B")] == [("c1", 6)]


def test_missing_n50_file_is_empty(tmp_path):
    assert read_n50(str(tmp_path / "nope.tsv")) == {}


def test_cli_runs_as_a_script(tmp_path):
    write_fixture(tmp_path)
    out = tmp_path / "x.json"
    subprocess.run([sys.executable, str(BIN / "pangenome_island_loci.py"),
                    "--islands_with_domains", str(tmp_path / "islands.tsv"),
                    "--presence_matrix", str(tmp_path / "matrix.tsv"),
                    "--family_positions", str(tmp_path / "family_positions.tsv"),
                    "--project", "demo", "--output", str(out)], check=True)
    assert json.loads(out.read_text())["project"] == "demo"


def test_without_config_all_strains_are_one_unknown_species_group(tmp_path):
    # Review Focus 5: no --config means no species map; the breakpoint track
    # then has one "unknown species" group and nothing crashes.
    write_fixture(tmp_path)
    out = tmp_path / "o.json"
    assert main(["--islands_with_domains", str(tmp_path / "islands.tsv"),
                 "--presence_matrix", str(tmp_path / "matrix.tsv"),
                 "--family_positions", str(tmp_path / "family_positions.tsv"),
                 "--project", "demo", "--output", str(out)]) == 0
    (locus,) = json.loads(out.read_text())["loci"]
    assert list(locus["counts_by_species"]) == [""]
    assert all(set(bp["indel"]) <= {""} for bp in locus["breakpoints"])
```


- [ ] **Step 2: Run**

```bash
pixi run python -m pytest -q tests/test_pangenome_island_loci.py
```

Expected: FAIL: `ModuleNotFoundError: No module named 'pangenome_island_loci'`


- [ ] **Step 3: Make the change**

Create `bin/pangenome_island_loci.py`:

```python
#!/usr/bin/env python3
"""Compute the island locus view (spec
docs/superpowers/specs/2026-09-24-island-locus-view-design.md, sections 1-6)
and write it as island_loci.json for bin/pangenome_island_synteny.py.

family_positions.tsv[.zst] is streamed three times, each pass filtered:
  1. copies of the candidate loci's member families + every contig's
     (min, max) rank, to find carriers and choose each exemplar;
  2. the exemplars' own contigs, to read their gene order (the columns);
  3. copies of every column family, for the per-strain cell states.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import open_maybe_compressed  # noqa: E402
from config_parser import parse_config  # noqa: E402
from island_locus import (  # noqa: E402
    DEFAULT_CONTAINMENT, DEFAULT_EMPTY_FRAC, DEFAULT_FLANK, DEFAULT_FLANK_MIN, DEFAULT_K,
    candidate_loci, carrier_placements, choose_exemplar, compute_locus, group_loci,
    locus_columns, locus_payload, rank_key,
)
from pangenome_matrix import PresenceMatrix  # noqa: E402
from pfam_classes import dominant_class  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from pangenome_domain_enrichment import parse_domtblout  # noqa: E402
from pangenome_island_synteny import load_gene_locations  # noqa: E402
from island_synteny import column_locations  # noqa: E402


@dataclass
class PositionScan:
    positions: dict = field(default_factory=dict)  # {(strain, family): [(contig, rank)]}
    orders: dict = field(default_factory=dict)  # {(strain, contig): [(rank, family)]}
    spans: dict = field(default_factory=dict)  # {(strain, contig): (min_rank, max_rank)}


def scan_family_positions(path: str, families: set[str] | None = None,
                          contigs: set[tuple[str, str]] | None = None,
                          spans: bool = False) -> PositionScan:
    """One streaming pass over family_positions.tsv[.zst] (Short, family,
    contig, rank). Keeps every copy of a family in `families`, the full
    order of each (strain, contig) in `contigs`, and with `spans` the
    (min, max) rank of every (strain, contig)."""
    out = PositionScan()
    with open_maybe_compressed(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_s, i_f, i_c, i_r = (header.index(c) for c in ("Short", "family", "contig", "rank"))
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_s, i_f, i_c, i_r):
                continue
            try:
                rank = int(parts[i_r])
            except ValueError:
                continue
            strain, fam, contig = parts[i_s], parts[i_f], parts[i_c]
            if families is not None and fam in families:
                out.positions.setdefault((strain, fam), []).append((contig, rank))
            if contigs is not None and (strain, contig) in contigs:
                out.orders.setdefault((strain, contig), []).append((rank, fam))
            if spans:
                key = (strain, contig)
                cur = out.spans.get(key)
                if cur is None:
                    out.spans[key] = (rank, rank)
                elif rank < cur[0] or rank > cur[1]:
                    out.spans[key] = (min(cur[0], rank), max(cur[1], rank))
    return out


def read_islands(path: str) -> list[dict]:
    with open_maybe_compressed(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_matrix_strains(path: str) -> list[str]:
    with open_maybe_compressed(path) as fh:
        return fh.readline().rstrip("\n").split("\t")[1:]


def read_n50(path: str | None) -> dict[str, int]:
    """{Short: n50} from assembly_quality_vs_content.tsv. Missing or empty
    file gives {}."""
    if not path or not Path(path).is_file() or Path(path).stat().st_size == 0:
        return {}
    out = {}
    with open_maybe_compressed(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                out[row["Short"]] = int(row["n50"])
            except (KeyError, ValueError):
                continue
    return out


def read_bins(path: str | None) -> dict[str, str]:
    if not path:
        return {}
    with open_maybe_compressed(path) as fh:
        return {r["family"]: r.get("bin", "") for r in csv.DictReader(fh, delimiter="\t")}


def build(args) -> dict:
    islands = read_islands(args.islands_with_domains)
    strains = read_matrix_strains(args.presence_matrix)
    loci = group_loci(islands, args.containment)
    cands = candidate_loci(loci, len(strains), args.rank_by, args.candidates, args.min_strains)
    n50 = read_n50(args.assembly_quality)
    species_of = ({s.short: s.species for s in parse_config(args.config)} if args.config else {})

    members = {f for loc in cands for f in loc.members}
    scan1 = scan_family_positions(args.family_positions, families=members, spans=True)
    chosen = []
    n_unplaced = 0
    for loc in cands:
        places = carrier_placements(loc.members, scan1.positions, strains, scan1.spans, args.k)
        pick = choose_exemplar(places, n50, args.flank, args.flank_min)
        if pick is None:
            n_unplaced += 1
            continue
        chosen.append((loc, pick[0], pick[1]))

    scan2 = scan_family_positions(args.family_positions,
                                  contigs={(p.strain, p.contig) for _, p, _ in chosen})
    with_cols = []
    for loc, place, tier in chosen:
        flank = args.flank_min if tier == "short_flanks" else args.flank
        cols = locus_columns(scan2.orders.get((place.strain, place.contig), []), loc.members,
                             flank, block=(place.lo, place.hi))
        if cols is None:
            n_unplaced += 1
            continue
        with_cols.append((loc, place, tier, cols))

    col_fams = {f for _, _, _, cols in with_cols for f in cols.families}
    scan3 = scan_family_positions(args.family_positions, families=col_fams)
    matrix = PresenceMatrix.from_tsv(args.presence_matrix, families=col_fams)
    results = [compute_locus(loc, place, tier, cols, strains, scan3.positions, scan1.spans,
                             matrix, species_of, args.k, args.empty_frac)
               for loc, place, tier, cols in with_cols]
    results.sort(key=lambda r: rank_key(r, args.rank_by))
    drawn = results[:args.top_loci]

    bins = read_bins(args.frequency_table)
    fam_domains = (parse_domtblout([args.domtblout], max_ievalue=args.domain_evalue)
                   if args.domtblout else {})
    gene_locs = None
    if args.gene_positions and args.cluster_tsv:
        gene_locs = load_gene_locations(args.cluster_tsv, args.gene_positions,
                                        {f for r in drawn for f in r["families"]},
                                        {r["exemplar"] for r in drawn}, args.id_sep)
    out_loci = []
    for i, r in enumerate(drawn):
        dom_sets = [sorted(fam_domains.get(f, set())) for f in r["families"]]
        locs = span = None
        if gene_locs is not None:
            locs = column_locations(r["families"], r["exemplar"], r["exemplar_contig"], -1, -1,
                                    gene_locs)
            block = [x for x in locs[r["n_left"]:r["n_left"] + r["n_locus"]]
                     if x and not x.get("rescued")]
            if block:
                span = (min(x["start"] for x in block), max(x["end"] for x in block))
        out_loci.append(locus_payload(
            r, "L%03d" % (i + 1), bins,
            [dominant_class(d) for d in dom_sets], [",".join(d) for d in dom_sets],
            dominant_class(sorted({d for ds in dom_sets for d in ds})), locs, span))
    return {
        "project": args.project,
        "locus_params": {"flank": args.flank, "flank_min": args.flank_min, "k": args.k,
                         "empty_frac": args.empty_frac, "containment": args.containment,
                         "rank_by": args.rank_by, "top_loci": args.top_loci,
                         "candidates": args.candidates, "min_strains": args.min_strains},
        "n_loci_total": len(loci),
        "n_loci_candidates": len(cands),
        "n_loci_unplaced": n_unplaced,
        "n_strains_with_n50": sum(1 for s in strains if s in n50),
        "n_strains": len(strains),
        "loci": out_loci,
        "_results": drawn,
    }


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--islands_with_domains", required=True)
    ap.add_argument("--presence_matrix", required=True)
    ap.add_argument("--family_positions", required=True)
    ap.add_argument("--project", required=True)
    ap.add_argument("--frequency_table", default=None)
    ap.add_argument("--assembly_quality", default=None,
                    help="assembly_quality_vs_content.tsv (Short, n50, ...) for exemplar ties")
    ap.add_argument("--config", default=None, help="samplesheet CSV, for species")
    ap.add_argument("--domtblout", default=None)
    ap.add_argument("--domain_evalue", type=float, default=1e-3)
    ap.add_argument("--gene_positions", default=None)
    ap.add_argument("--cluster_tsv", default=None)
    ap.add_argument("--id_sep", default="|")
    ap.add_argument("--flank", type=int, default=DEFAULT_FLANK)
    ap.add_argument("--flank_min", type=int, default=DEFAULT_FLANK_MIN)
    ap.add_argument("--k", type=int, default=DEFAULT_K)
    ap.add_argument("--empty_frac", type=float, default=DEFAULT_EMPTY_FRAC)
    ap.add_argument("--containment", type=float, default=DEFAULT_CONTAINMENT)
    ap.add_argument("--rank_by", choices=["informative", "strains", "size"], default="informative")
    ap.add_argument("--top_loci", type=int, default=50)
    ap.add_argument("--candidates", type=int, default=200)
    ap.add_argument("--min_strains", type=int, default=2)
    ap.add_argument("--output", required=True)
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    payload = build(args)
    payload.pop("_results")
    Path(args.output).write_text(json.dumps(payload, separators=(",", ":")))
    print(f"pangenome_island_loci: {payload['n_loci_total']} loci, "
          f"{payload['n_loci_candidates']} candidates, {payload['n_loci_unplaced']} without "
          f"a carrier, {len(payload['loci'])} drawn (--rank_by {args.rank_by}); wrote "
          f"{args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

Make it executable: `chmod +x bin/pangenome_island_loci.py`.


- [ ] **Step 4: Run**

```bash
pixi run python -m pytest -q tests/test_pangenome_island_loci.py
```

Expected: 11 passed


- [ ] **Step 5: Commit**

```bash
git add bin/pangenome_island_loci.py tests/test_pangenome_island_loci.py
git commit -m "island locus view: pangenome_island_loci.py writes island_loci.json (#${ISSUE_A})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 8: Regression against the feasibility prototype (spec validation item 2)

**Files:**
- Create: `tests/test_island_locus_regression.py`

**Interfaces:**
- Consumes: `strain_cells`, `flank_pair`, `row_class`, `locus_columns`, `group_loci` (Tasks 1-5), `scan_family_positions`, `read_matrix_strains` (Task 7).
- Produces: an opt-in test; no new library names.


- [ ] **Step 1: Read this first**

The expected numbers are the spec's Feasibility table, re-measured on 2026-09-26 by running the prototype itself (`island_locus_view_feasibility.py --flank 5 --k 10 --top 50`, `--sort size|strains`, `--empty_frac 0.8|1.0`). This implementation already passed this test during planning (3 passed, 67 s, 238 MB peak RSS).


- [ ] **Step 2: Write the test**

Create `tests/test_island_locus_regression.py`:

```python
"""Regression: lib/island_locus.py reproduces the feasibility prototype.

Spec validation plan item 2. The prototype is
NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/analysis/
island_locus_view_feasibility.py; its numbers on the 529-strain run
`rescue_freqpol_immitis_in_posadasii_out` (NovInvenio 91e3157) are in the
spec's Feasibility table and are pinned below (re-measured 2026-09-26 with
the prototype, --flank 5 --k 10 --top 50).

The prototype's definitions differ from the pipeline's in three stated ways,
so this test drives the library in "prototype mode": islands (not loci) are
the units, the exemplar is the island's example_strain, the locus block is
every member copy on the exemplar contig (locus_columns(block=None)), and
flanks intact accepts any copy pair (flank_pair(in_place=None)). Cell states
use StrainCells.base, which has no rescue or contig-break state.

Opt-in: set NOVINVENIO_LOCUS_REGRESSION_RUN to that run's output/pangenome
directory. About 3 minutes (two sort orders, two passes each).
"""
import collections
import csv
import os
import statistics
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
from island_locus import (  # noqa: E402
    ABSENT, IN_PLACE, flank_pair, group_loci, locus_columns, row_class, strain_cells,
)
from pangenome_island_loci import read_matrix_strains, scan_family_positions  # noqa: E402

RUN = os.environ.get("NOVINVENIO_LOCUS_REGRESSION_RUN")
pytestmark = pytest.mark.skipif(not RUN, reason="set NOVINVENIO_LOCUS_REGRESSION_RUN")

EXPECTED = {
    "size": {
        "both_sides_ge5": 21, "touch_end": 21, "median_left": 6.0, "median_right": 7.0,
        "core_softcore_pct": 66, "in_place": 131158, "elsewhere": 3965, "absent": 304476,
        "median_intact": 61.0, "median_empty_080": 0.0, "median_full": 2.0,
        "informative_080": 17, "informative_100": 0,
    },
    "strains": {
        "both_sides_ge5": 30, "touch_end": 7, "median_left": 15.0, "median_right": 16.0,
        "core_softcore_pct": 76, "in_place": 50604, "elsewhere": 297, "absent": 22101,
        "median_intact": 500.5, "median_empty_080": 52.5, "median_full": 245.5,
        "informative_080": 39, "informative_100": 39,
    },
}


def _islands():
    with open(Path(RUN) / "report_tables" / "islands_with_domains.tsv") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["locus_id"] not in ("-", "")]
    for r in rows:
        r["members"] = [f for f in r["member_families"].split(",") if f]
    return rows


def _measure(sort: str) -> dict:
    run = Path(RUN)
    fp = str(run / "family_positions.tsv.zst")
    qual = [r for r in _islands() if int(r["n_strains"]) >= 2]
    if sort == "size":
        qual.sort(key=lambda r: -int(r["island_size"]))
    else:
        qual.sort(key=lambda r: (-int(r["n_strains"]), -int(r["island_size"])))
    top = qual[:50]
    with open(run / "frequency_table.tsv") as fh:
        bins = {r["family"]: r["bin"] for r in csv.DictReader(fh, delimiter="\t")}
    scan = scan_family_positions(fp, contigs={(r["example_strain"], r["locus_contig"]) for r in top})
    cols = []
    for r in top:
        c = locus_columns(scan.orders.get((r["example_strain"], r["locus_contig"]), []),
                          r["members"], 5)
        if c is not None:
            cols.append(c)
    need = {f for c in cols for f in c.families}
    pos = scan_family_positions(fp, families=need).positions
    strains = read_matrix_strains(str(run / "presence_matrix.rescued.tsv"))
    lefts = [c.left_avail for c in cols]
    rights = [c.right_avail for c in cols]
    flank_bins = collections.Counter(bins.get(f, "?") for c in cols for f in c.left + c.right)
    states = collections.Counter()
    intact_n, empty80, empty100, full = [], [], [], []
    for c in cols:
        fams = c.families
        nl, nb = len(c.left), len(c.locus)
        n_int = e80 = e100 = nf = 0
        for s in strains:
            base = strain_cells(s, fams, pos, {}).base
            states.update(base[nl:nl + nb])
            pair = flank_pair(pos, s, list(c.left), list(c.right), nb) if c.left and c.right else None
            if pair is None:
                continue
            n_int += 1
            if row_class(base, nl, nb, True, 0.8) == "empty":
                e80 += 1
            if row_class(base, nl, nb, True, 1.0) == "empty":
                e100 += 1
            if row_class(base, nl, nb, True, 0.8) == "full":
                nf += 1
        intact_n.append(n_int)
        empty80.append(e80)
        empty100.append(e100)
        full.append(nf)
    tot_bins = sum(flank_bins.values())
    return {
        "both_sides_ge5": sum(1 for a, b in zip(lefts, rights) if a >= 5 and b >= 5),
        "touch_end": sum(1 for a, b in zip(lefts, rights) if a == 0 or b == 0),
        "median_left": statistics.median(lefts), "median_right": statistics.median(rights),
        "core_softcore_pct": round((flank_bins["core"] + flank_bins["soft_core"]) / tot_bins * 100),
        "in_place": states[IN_PLACE], "elsewhere": states["2"], "absent": states[ABSENT],
        "median_intact": statistics.median(intact_n),
        "median_empty_080": statistics.median(empty80), "median_full": statistics.median(full),
        "informative_080": sum(1 for e, f in zip(empty80, full) if e >= 10 and f >= 2),
        "informative_100": sum(1 for e, f in zip(empty100, full) if e >= 10 and f >= 2),
    }


@pytest.mark.parametrize("sort", ["size", "strains"])
def test_matches_feasibility_prototype(sort):
    assert _measure(sort) == EXPECTED[sort]


def test_locus_grouping_matches_m4():
    loci = group_loci(_islands())
    joined = sum(len(loc.variants) - 1 for loc in loci)
    assert (joined, sum(len(loc.variants) for loc in loci)) == (10579, 11880)
```


- [ ] **Step 3: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus_regression.py
```

Expected: 3 skipped (the run directory is not set)


- [ ] **Step 4: Run**

```bash
NOVINVENIO_LOCUS_REGRESSION_RUN=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome pixi run python -m pytest -q tests/test_island_locus_regression.py
```

Expected: 3 passed (about 70 s). A failure prints the differing key: fix the library, never the expected numbers.


- [ ] **Step 5: Commit**

```bash
git add tests/test_island_locus_regression.py
git commit -m "island locus view: regression test against the feasibility prototype (#${ISSUE_A})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 9: Page: locus view fragments in `island_synteny.html` (spec section 7)

**Files:**
- Create: `lib/island_locus_template.py`
- Modify: `lib/island_synteny_template.py`
- Create: `tests/test_island_locus_template.py`

**Interfaces:**
- Consumes: the payload keys `loci`, `locus_meta.locus_params`, `assembly_confound` (Task 10 adds them); the island IIFE's `el css palette classColor classLabel sizeCanvas ellipsize positionTip tipEl appendLocation strainPopupLines speciesCounts haplotypeSpecies SPECIES ROW_LABEL_FONT GUTTER_MIN GUTTER_MAX LABEL_FONT LABEL_MAX_W LABEL_ANGLE GLYPH_H state renderMain`.
- Produces: `LOCUS_VIEW_CSS`, `LOCUS_VIEW_HTML`, `LOCUS_VIEW_JS` (strings); DOM ids `lv-switch lv-btn-loci lv-btn-islands locus-view island-view lv-list lv-search lv-sort lv-count lv-title lv-tier lv-note lv-legend lv-track lv-head lv-grid lv-confound`; pure JS helpers `locusStateStyle locusStateLabel locusClassLabel locusSortedRows locusSidebarOrder locusSpeciesList breakpointBars cellReasonLines locusTitle locusTierText`. `renderLocusMain()` calls `renderClinkerPanel(locus)` only when that function exists (Task 21 defines it).


- [ ] **Step 1: Write the failing test**

Create `tests/test_island_locus_template.py`:

```python
"""Locus view JS in island_synteny.html: pure helpers run under node, plus
structural checks. Same brace-matching extraction as
tests/test_island_synteny.py."""
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from island_locus_template import LOCUS_VIEW_CSS, LOCUS_VIEW_HTML, LOCUS_VIEW_JS  # noqa: E402
from island_synteny_template import ISLAND_SYNTENY_TEMPLATE  # noqa: E402


def _extract_function(src: str, name: str) -> str:
    start = src.index("function " + name + "(")
    depth = 0
    for i in range(src.index("{", start), len(src)):
        if src[i] == "{":
            depth += 1
        elif src[i] == "}":
            depth -= 1
            if depth == 0:
                return src[start:i + 1]
    raise AssertionError(f"unterminated function {name}()")


def run_node(names, cases_js):
    node = shutil.which("node")
    if not node:
        pytest.skip("node not available")
    funcs = "\n".join(_extract_function(ISLAND_SYNTENY_TEMPLATE, n) for n in names)
    proc = subprocess.run([node, "--input-type=commonjs", "-e", funcs + "\n" + cases_js],
                          capture_output=True, text=True, check=False)
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout.strip().splitlines()[-1])


def test_page_embeds_the_three_fragments():
    assert LOCUS_VIEW_CSS in ISLAND_SYNTENY_TEMPLATE
    assert LOCUS_VIEW_HTML in ISLAND_SYNTENY_TEMPLATE
    assert LOCUS_VIEW_JS in ISLAND_SYNTENY_TEMPLATE
    assert ISLAND_SYNTENY_TEMPLATE.index(LOCUS_VIEW_JS) < ISLAND_SYNTENY_TEMPLATE.index("// ---- init")


def test_island_view_is_wrapped_so_it_can_be_hidden():
    page = ISLAND_SYNTENY_TEMPLATE
    assert page.index('<div id="island-view">') < page.index('id="isv-body"')


def test_locus_fragments_use_no_innerhtml_and_no_hex_colours():
    assert "innerHTML" not in LOCUS_VIEW_JS
    for frag in (LOCUS_VIEW_CSS, LOCUS_VIEW_HTML, LOCUS_VIEW_JS):
        assert not re.findall(r"(?<![\w-])#[0-9a-fA-F]{3,8}(?![\w-])", frag)


def test_state_styles_encode_evidence_by_token():
    out = run_node(["locusStateStyle"], """
      console.log(JSON.stringify(["0","1","2","3","4","5"].map(locusStateStyle)));""")
    assert [o["token"] for o in out] == ["--grid", "--series-1", "--series-1", "--series-2",
                                         "--series-2", "--warn"]
    assert [o["hatch"] for o in out] == [False, False, False, True, True, True]
    assert out[2]["alpha"] < 1 and out[1]["alpha"] == 1


def test_rows_group_by_class_then_species_then_count():
    rows = [
        {"row_class": "empty", "codes": "10", "count": 5, "strains": ["e1"]},
        {"row_class": "full", "codes": "11", "count": 1, "strains": ["b1"]},
        {"row_class": "full", "codes": "12", "count": 9, "strains": ["a1"]},
        {"row_class": "uninformative", "codes": "00", "count": 50, "strains": ["u1"]},
        {"row_class": "full", "codes": "13", "count": 3, "strains": ["x1"]},
    ]
    species = {"a1": "sp2", "b1": "sp1", "e1": "sp1", "u1": "sp1"}
    out = run_node(["speciesCounts", "haplotypeSpecies", "locusSortedRows"], f"""
      var rows = {json.dumps(rows)};
      console.log(JSON.stringify(locusSortedRows(rows, {json.dumps(species)}).map(function (r) {{ return r.strains[0]; }})));""")
    assert out == ["b1", "a1", "x1", "e1", "u1"]


def test_sidebar_orders():
    loci = [{"informative_score": 5, "n_carriers": 9, "size": 2, "locus_id": "b"},
            {"informative_score": 40, "n_carriers": 3, "size": 7, "locus_id": "a"},
            {"informative_score": -1, "n_carriers": 99, "size": 1, "locus_id": "c"}]
    out = run_node(["locusSidebarOrder"], f"""
      var L = {json.dumps(loci)};
      console.log(JSON.stringify(["informative","strains","size","name"].map(function (k) {{
        return locusSidebarOrder(L, [0, 1, 2], k); }})));""")
    assert out == [[1, 0, 2], [2, 0, 1], [1, 0, 2], [1, 0, 2]]


def test_breakpoint_bars_are_one_per_species_then_contig_break():
    out = run_node(["breakpointBars"], """
      console.log(JSON.stringify(breakpointBars({b: 2, indel: {sp2: 4}, contig_break: 3}, ["sp1", "sp2"])));""")
    assert out == [{"kind": "indel", "species": "sp1", "n": 0},
                   {"kind": "indel", "species": "sp2", "n": 4},
                   {"kind": "contig_break", "species": "", "n": 3}]


def test_cell_reason_names_the_neighbour_and_rank():
    locus = {"families": ["F1", "A", "B"]}
    row = {"codes": "112", "count": 3, "strains": ["S1", "S2", "S3"],
           "rep": {"c": ["scaffold_9"], "d": [[0, 40, 1, 1], [0, 41, 0, -1], [0, 90]]}}
    out = run_node(["locusStateLabel", "cellReasonLines"], f"""
      var locus = {json.dumps(locus)}, row = {json.dumps(row)};
      console.log(JSON.stringify([cellReasonLines(locus, row, 1, 10), cellReasonLines(locus, row, 2, 10)]));""")
    assert out[0] == ["State: in place", "In S1 (first of 3 strains): scaffold_9, gene rank 41",
                      "Neighbour F1 at rank 40 (-1)"]
    assert out[1][2] == "No other column within 10 genes on this contig"


def test_cell_reason_for_contig_break_and_missing_detail():
    out = run_node(["locusStateLabel", "cellReasonLines"], """
      var locus = {families: ["F1"]};
      var r1 = {codes: "5", count: 1, strains: ["S1"], rep: {c: [], d: [null]}};
      var r2 = {codes: "0", count: 1, strains: ["S1"], rep: null};
      console.log(JSON.stringify([cellReasonLines(locus, r1, 0, 10), cellReasonLines(locus, r2, 0, 10)]));""")
    assert "not evidence" in out[0][1]
    assert out[1][1] == "Position detail not stored for this row"


def test_tier_text():
    out = run_node(["locusTierText"], """
      var p = {flank_min: 3};
      console.log(JSON.stringify([locusTierText("full", p), locusTierText("short_flanks", p), locusTierText("contig_end", p)]));""")
    assert out == ["", "short flanks (3 genes per side)", "exemplar at contig end"]

```


- [ ] **Step 2: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus_template.py
```

Expected: FAIL: `ModuleNotFoundError: No module named 'island_locus_template'`


- [ ] **Step 3: Make the change**

Create `lib/island_locus_template.py`:

```python
"""Locus view fragments for island_synteny.html (spec
docs/superpowers/specs/2026-09-24-island-locus-view-design.md, section 7).

lib/island_synteny_template.py inserts three fragments:
  LOCUS_VIEW_CSS   into the page <style>
  LOCUS_VIEW_HTML  before the island view, inside .wrap
  LOCUS_VIEW_JS    inside the page's one IIFE, after the island-view code,
                   so it reuses el(), css(), palette(), classColor(),
                   classLabel(), sizeCanvas(), ellipsize(), positionTip(),
                   tipEl, appendLocation(), strainPopupLines(),
                   speciesCounts(), haplotypeSpecies() and SPECIES.

The page shows the locus view by default when the payload has `loci`; the
island (presence) view stays one click away (Ruling R1). Colours come from
skin tokens only: --series-1 = in place (protein-search presence),
--series-2 = TBLASTN rescue, --grid = absent, --warn = contig break.
Species are told apart by position and label, never by hue.
"""

LOCUS_VIEW_CSS = r"""
  .lv-switch { display: inline-flex; gap: 0; margin: 0 0 14px; border: 1px solid var(--border); border-radius: 8px; overflow: hidden; }
  .lv-switch button { border: none; border-radius: 0; padding: 6px 14px; background: var(--surface-1); }
  .lv-switch button[aria-pressed="true"] { background: var(--wash); box-shadow: inset 0 -2px 0 var(--series-1); font-weight: 600; }
  .lv-badge { display: inline-block; margin-left: 6px; padding: 1px 7px; border-radius: 999px; font-size: 10.5px; border: 1px solid var(--warn); color: var(--text-primary); }
  .lv-track-wrap { border-bottom: 1px solid var(--border); }
  canvas.lv-track, canvas.lv-head { display: block; }
  .lv-confound { margin: 0 0 12px; font-size: 12px; color: var(--text-secondary); }
"""

LOCUS_VIEW_HTML = r"""
  <p class="lv-confound hidden" id="lv-confound"></p>
  <div class="lv-switch hidden" id="lv-switch" role="group" aria-label="Choose view">
    <button type="button" id="lv-btn-loci" aria-pressed="true">Loci (flank-anchored)</button>
    <button type="button" id="lv-btn-islands" aria-pressed="false">Islands (presence)</button>
  </div>
  <div id="locus-view" class="hidden">
    <div class="isv-explorer">
      <aside class="isv-sidebar">
        <div class="isv-sidebar-controls">
          <input type="search" id="lv-search" placeholder="Search locus, strain or family ID…" aria-label="Search loci">
          <select id="lv-sort" aria-label="Sort loci by">
            <option value="informative">Sort: informative polymorphism</option>
            <option value="strains">Sort: strains with the locus</option>
            <option value="size">Sort: size (families)</option>
            <option value="name">Sort: locus ID</option>
          </select>
          <span class="count" id="lv-count" role="status" aria-live="polite"></span>
        </div>
        <div class="isv-list" id="lv-list" role="listbox" aria-label="Accessory loci"></div>
      </aside>
      <section class="card" id="lv-main" style="padding:14px">
        <div class="isv-main-head">
          <h2 id="lv-title"></h2>
          <span class="lv-badge hidden" id="lv-tier"></span>
        </div>
        <p class="isv-main-note" id="lv-note"></p>
        <div class="isv-legend" id="lv-legend"></div>
        <div class="isv-hscroll" id="lv-hscroll">
          <div class="lv-track-wrap"><canvas class="lv-track" id="lv-track"></canvas></div>
          <canvas class="lv-head" id="lv-head"></canvas>
          <div class="isv-vscroll"><canvas class="isv-grid" id="lv-grid"></canvas></div>
        </div>
      </section>
    </div>
  </div>
"""

LOCUS_VIEW_JS = r"""
  // ==== locus view ==========================================================
  var LOCI = DATA.loci || [];
  var LMETA = DATA.locus_meta || {};
  var LPARAMS = LMETA.locus_params || {};
  var lstate = { search: "", sort: "informative", selected: LOCI.length ? 0 : -1 };
  var lSidebar = [];
  var lRows = [];
  var L_ROW_H = 20, L_CELL_W = 22, L_TRACK_H = 64, L_GUTTER = 190;
  var lTrack = document.getElementById("lv-track");
  var lHead = document.getElementById("lv-head");
  var lGrid = document.getElementById("lv-grid");
  var ltctx = lTrack.getContext("2d");
  var lhctx = lHead.getContext("2d");
  var lgctx = lGrid.getContext("2d");

  // Pure helpers: no DOM, no module state, so tests can run them under node.
  function locusStateStyle(code) {
    if (code === "1") return { token: "--series-1", alpha: 1, hatch: false };
    if (code === "2") return { token: "--series-1", alpha: 0.35, hatch: false };
    if (code === "3") return { token: "--series-2", alpha: 1, hatch: true };
    if (code === "4") return { token: "--series-2", alpha: 0.35, hatch: true };
    if (code === "5") return { token: "--warn", alpha: 0.55, hatch: true };
    return { token: "--grid", alpha: 1, hatch: false };
  }
  function locusStateLabel(code) {
    var labels = { "0": "absent", "1": "in place", "2": "elsewhere",
      "3": "rescue, in place (TBLASTN hit, no annotated gene)",
      "4": "rescue, elsewhere (TBLASTN hit, no annotated gene)", "5": "contig break" };
    return labels[code] || "unknown";
  }
  function locusClassLabel(cls) {
    var labels = { full: "full locus", partial: "partial", empty: "empty site", uninformative: "uninformative" };
    return labels[cls] || cls;
  }
  function locusSortedRows(rows, speciesMap) {
    var order = ["full", "partial", "empty", "uninformative"];
    return rows.slice().sort(function (a, b) {
      var ca = order.indexOf(a.row_class), cb = order.indexOf(b.row_class);
      if (ca !== cb) return ca - cb;
      var sa = haplotypeSpecies(a, speciesMap) || "￿";
      var sb = haplotypeSpecies(b, speciesMap) || "￿";
      if (sa !== sb) return sa < sb ? -1 : 1;
      if (b.count !== a.count) return b.count - a.count;
      return a.codes < b.codes ? -1 : a.codes > b.codes ? 1 : 0;
    });
  }
  function locusSidebarOrder(loci, idxs, key) {
    var cmp;
    if (key === "strains") cmp = function (a, b) { return loci[b].n_carriers - loci[a].n_carriers || a - b; };
    else if (key === "size") cmp = function (a, b) { return loci[b].size - loci[a].size || a - b; };
    else if (key === "name") cmp = function (a, b) { return loci[a].locus_id < loci[b].locus_id ? -1 : loci[a].locus_id > loci[b].locus_id ? 1 : 0; };
    else cmp = function (a, b) { return loci[b].informative_score - loci[a].informative_score || loci[b].n_carriers - loci[a].n_carriers || a - b; };
    return idxs.slice().sort(cmp);
  }
  function locusSpeciesList(locus) {
    return Object.keys(locus.counts_by_species || {}).sort();
  }
  function breakpointBars(bp, speciesList) {
    var bars = speciesList.map(function (sp) {
      return { kind: "indel", species: sp, n: (bp.indel && bp.indel[sp]) || 0 };
    });
    bars.push({ kind: "contig_break", species: "", n: bp.contig_break || 0 });
    return bars;
  }
  function cellReasonLines(locus, row, ci, k) {
    var code = row.codes.charAt(ci);
    var lines = ["State: " + locusStateLabel(code)];
    var who = row.strains[0] + (row.count > 1 ? " (first of " + row.count + " strains)" : "");
    if (!row.rep) { lines.push("Position detail not stored for this row"); return lines; }
    var d = row.rep.d[ci];
    if (!d) {
      lines.push(code === "5"
        ? "No copy at the locus in " + who + "; the nearest placed column is within " + k + " genes of a contig end, so this absence is not evidence"
        : (code === "2" ? "Present in " + who + " but no position recorded" : "No copy in " + who));
      return lines;
    }
    lines.push("In " + who + ": " + row.rep.c[d[0]] + ", gene rank " + d[1]);
    if (d.length > 2) {
      lines.push("Neighbour " + locus.families[d[2]] + " at rank " + (d[1] + d[3]) +
        " (" + (d[3] > 0 ? "+" : "") + d[3] + ")");
    } else {
      lines.push("No other column within " + k + " genes on this contig");
    }
    return lines;
  }
  function locusTitle(locus) {
    var span = locus.exemplar_span;
    return locus.exemplar + ":" + locus.exemplar_contig +
      (span ? ":" + span.start + "-" + span.end : "");
  }
  function locusTierText(tier, params) {
    if (tier === "short_flanks") return "short flanks (" + params.flank_min + " genes per side)";
    if (tier === "contig_end") return "exemplar at contig end";
    return "";
  }

  function lcolX(i) { return L_GUTTER + i * L_CELL_W; }
  function lTotalWidth(locus) { return lcolX(locus.families.length) + 8; }
  function isFlank(locus, i) { return i < locus.n_left || i >= locus.n_left + locus.n_locus; }

  function locusHaystack(locus) {
    return (locus.locus_id + " " + locus.exemplar + " " + locus.families.join(" ")).toLowerCase();
  }
  function applyLocusFilter() {
    var q = lstate.search.trim().toLowerCase();
    var idxs = [];
    for (var i = 0; i < LOCI.length; i++) {
      if (q && locusHaystack(LOCI[i]).indexOf(q) === -1) continue;
      idxs.push(i);
    }
    lSidebar = locusSidebarOrder(LOCI, idxs, lstate.sort);
  }
  function renderLocusSidebar() {
    var list = document.getElementById("lv-list");
    list.textContent = "";
    if (!lSidebar.length) {
      list.appendChild(el("p", "isv-empty-list", "No loci match the current search."));
    }
    lSidebar.forEach(function (idx) {
      var locus = LOCI[idx];
      var btn = el("button", "isv-item");
      btn.type = "button";
      btn.setAttribute("role", "option");
      btn.setAttribute("aria-selected", idx === lstate.selected ? "true" : "false");
      if (idx === lstate.selected) btn.classList.add("sel");
      btn.appendChild(el("div", "isv-item-id", locus.locus_id));
      var c = locus.counts;
      btn.appendChild(el("div", "isv-item-stats",
        locus.size + " families · " + locus.n_variants + " variants · " +
        "empty " + c.empty + " · full " + c.full + " · partial " + c.partial +
        " · uninformative " + c.uninformative));
      var chip = el("span", "isv-chip", classLabel(locus.dominant_class));
      chip.style.borderColor = classColor(locus.dominant_class);
      chip.style.color = classColor(locus.dominant_class);
      btn.appendChild(chip);
      var tier = locusTierText(locus.tier, LPARAMS);
      if (tier) btn.appendChild(el("span", "lv-badge", tier));
      btn.addEventListener("click", function () { selectLocus(idx); });
      list.appendChild(btn);
    });
    document.getElementById("lv-count").textContent =
      lSidebar.length.toLocaleString() + " of " + LOCI.length.toLocaleString() + " shown";
  }
  function selectLocus(idx) {
    lstate.selected = idx;
    renderLocusSidebar();
    renderLocusMain();
  }

  function lRowLabel(row) {
    return (row.count === 1 ? row.strains[0] : row.strains[0] + " +" + (row.count - 1));
  }
  function locusGutter(rows) {
    lgctx.font = ROW_LABEL_FONT;
    var maxW = 0;
    rows.forEach(function (r) { maxW = Math.max(maxW, lgctx.measureText(lRowLabel(r)).width); });
    return Math.min(GUTTER_MAX, Math.max(GUTTER_MIN, maxW + 60 + 12));
  }

  function hatch(ctx, x, y, w, h, color) {
    ctx.save();
    ctx.strokeStyle = color;
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (var off = -h; off < w; off += 5) {
      ctx.moveTo(x + off, y + h);
      ctx.lineTo(x + off + h, y);
    }
    ctx.stroke();
    ctx.restore();
  }

  function drawLocusTrack(locus, totalW) {
    var P = palette();
    sizeCanvas(lTrack, ltctx, totalW, L_TRACK_H);
    ltctx.clearRect(0, 0, totalW, L_TRACK_H);
    ltctx.fillStyle = P.surface;
    ltctx.fillRect(0, 0, totalW, L_TRACK_H);
    var species = locusSpeciesList(locus);
    var maxN = 1;
    (locus.breakpoints || []).forEach(function (bp) {
      breakpointBars(bp, species).forEach(function (b) { maxN = Math.max(maxN, b.n); });
    });
    var base = L_TRACK_H - 14;
    ltctx.font = ROW_LABEL_FONT;
    ltctx.fillStyle = P.secondary;
    ltctx.textBaseline = "middle";
    ltctx.fillText("breakpoints (max " + maxN + ")", 8, 10);
    ltctx.fillText("species: " + (species.map(function (s) { return s || "unknown"; }).join(" | ")), 8, L_TRACK_H - 6);
    var nBars = species.length + 1;
    var barW = Math.max(2, Math.floor((L_CELL_W - 4) / nBars));
    (locus.breakpoints || []).forEach(function (bp) {
      var bars = breakpointBars(bp, species);
      var x0 = lcolX(bp.b) - (barW * nBars) / 2;
      bars.forEach(function (b, i) {
        if (!b.n) return;
        var h = Math.max(1, Math.round((base - 16) * b.n / maxN));
        ltctx.fillStyle = b.kind === "contig_break" ? css("--warn") : P.secondary;
        ltctx.fillRect(x0 + i * barW, base - h, barW - 1, h);
      });
    });
    ltctx.fillStyle = P.axis;
    ltctx.fillRect(0, base, totalW, 1);
  }

  function drawLocusHead(locus, totalW) {
    var P = palette();
    lhctx.font = LABEL_FONT;
    var maxW = 0;
    locus.families.forEach(function (f) { maxW = Math.max(maxW, Math.min(LABEL_MAX_W, lhctx.measureText(f).width)); });
    var H = Math.ceil(maxW * Math.sin(LABEL_ANGLE)) + GLYPH_H + 32;
    var W = totalW + Math.ceil(maxW * Math.cos(LABEL_ANGLE)) + 12;
    var glyphY = H - GLYPH_H - 10;
    sizeCanvas(lHead, lhctx, W, H);
    lhctx.clearRect(0, 0, W, H);
    lhctx.fillStyle = P.surface;
    lhctx.fillRect(0, 0, W, H);
    locus.families.forEach(function (fam, i) {
      var x = lcolX(i);
      lhctx.fillStyle = classColor((locus.family_classes && locus.family_classes[i]) || "unannotated");
      lhctx.fillRect(x + 3, glyphY, L_CELL_W - 6, GLYPH_H);
      if (isFlank(locus, i)) {
        lhctx.fillStyle = P.axis;
        lhctx.fillRect(x + 1, H - 6, L_CELL_W - 2, 3);
      }
      lhctx.save();
      lhctx.translate(x + L_CELL_W / 2, glyphY - 6);
      lhctx.rotate(-LABEL_ANGLE);
      lhctx.fillStyle = isFlank(locus, i) ? P.secondary : P.primary;
      lhctx.font = LABEL_FONT;
      lhctx.textAlign = "left";
      lhctx.textBaseline = "middle";
      lhctx.fillText(ellipsize(lhctx, fam, LABEL_MAX_W), 0, 0);
      lhctx.restore();
    });
    lhctx.font = ROW_LABEL_FONT;
    lhctx.fillStyle = P.secondary;
    lhctx.textBaseline = "middle";
    lhctx.fillText("flank = anchor (bar under column)", 8, H - 16);
  }

  function drawLocusGrid(locus, rows, totalW) {
    var P = palette();
    var h = Math.max(L_ROW_H, rows.length * L_ROW_H);
    sizeCanvas(lGrid, lgctx, totalW, h);
    lgctx.clearRect(0, 0, totalW, h);
    lgctx.fillStyle = P.surface;
    lgctx.fillRect(0, 0, totalW, h);
    var prevClass = null;
    rows.forEach(function (row, ri) {
      var y = ri * L_ROW_H;
      if (row.row_class !== prevClass) {
        if (prevClass !== null) { lgctx.fillStyle = P.axis; lgctx.fillRect(0, y, totalW, 1); }
        prevClass = row.row_class;
      }
      lgctx.textBaseline = "middle";
      lgctx.textAlign = "left";
      lgctx.font = "600 11px ui-monospace, SFMono-Regular, Menlo, monospace";
      lgctx.fillStyle = P.primary;
      lgctx.fillText("×" + row.count, 8, y + L_ROW_H / 2);
      lgctx.font = ROW_LABEL_FONT;
      lgctx.fillStyle = P.secondary;
      lgctx.fillText(ellipsize(lgctx, lRowLabel(row), L_GUTTER - 60), 50, y + L_ROW_H / 2);
      for (var ci = 0; ci < row.codes.length; ci++) {
        var st = locusStateStyle(row.codes.charAt(ci));
        var x = lcolX(ci) + 1, cw = L_CELL_W - 2, ch = L_ROW_H - 2;
        lgctx.fillStyle = P.grid;
        lgctx.fillRect(x, y + 1, cw, ch);
        lgctx.globalAlpha = st.alpha;
        lgctx.fillStyle = css(st.token);
        lgctx.fillRect(x, y + 1, cw, ch);
        lgctx.globalAlpha = 1;
        if (st.hatch) hatch(lgctx, x, y + 1, cw, ch, P.surface);
      }
    });
    lgctx.fillStyle = P.axis;
    lgctx.fillRect(lcolX(locus.n_left), 0, 1, h);
    lgctx.fillRect(lcolX(locus.n_left + locus.n_locus), 0, 1, h);
  }

  function renderLocusLegend() {
    var legend = document.getElementById("lv-legend");
    legend.textContent = "";
    ["1", "2", "3", "0", "5"].forEach(function (code) {
      var st = locusStateStyle(code);
      var item = el("span", "isv-legend-item");
      var sw = el("span", "isv-swatch");
      sw.style.background = css(st.token);
      sw.style.opacity = String(st.alpha);
      if (st.hatch) sw.style.backgroundImage = "repeating-linear-gradient(45deg, transparent 0 3px, " + css("--surface-1") + " 3px 4px)";
      item.appendChild(sw);
      item.appendChild(el("span", null, locusStateLabel(code).split(" (")[0]));
      legend.appendChild(item);
    });
    var bp = el("span", "isv-legend-item",
      "Bars above the columns: flank-intact strains changing between in place and absent, one bar per species (left to right as named); the last bar counts contig breaks");
    legend.appendChild(bp);
  }

  function renderLocusMain() {
    var locus = LOCI[lstate.selected];
    if (!locus) return;
    document.getElementById("lv-title").textContent = locusTitle(locus);
    var tierEl = document.getElementById("lv-tier");
    var tierText = locusTierText(locus.tier, LPARAMS);
    tierEl.textContent = tierText;
    tierEl.classList.toggle("hidden", !tierText);
    var c = locus.counts;
    document.getElementById("lv-note").textContent =
      "Exemplar " + locus.exemplar + ": carries the locus's largest variant with at least " +
      LPARAMS.flank + " genes on both sides (else " + LPARAMS.flank_min + "), ties by N50 then name. " +
      locus.n_left + " + " + locus.n_locus + " + " + locus.n_right +
      " columns (flank + locus + flank) in the exemplar's gene order. A cell is in place when " +
      "the strain has another column's gene within k = " + LPARAMS.k + " genes on the same contig. " +
      "Rows: strains with identical states, grouped by row class, then species. " +
      "Full " + c.full + ", partial " + c.partial + ", empty site " + c.empty +
      ", uninformative " + c.uninformative + " strains.";
    renderLocusLegend();
    lRows = locusSortedRows(locus.rows, SPECIES);
    L_GUTTER = locusGutter(lRows);
    var w = lTotalWidth(locus);
    drawLocusTrack(locus, w);
    drawLocusHead(locus, w);
    drawLocusGrid(locus, lRows, w);
    if (typeof renderClinkerPanel === "function") renderClinkerPanel(locus);
  }

  function lColAt(canvas, clientX, locus) {
    var x = clientX - canvas.getBoundingClientRect().left;
    for (var i = 0; i < locus.families.length; i++) {
      if (x >= lcolX(i) && x < lcolX(i) + L_CELL_W) return i;
    }
    return -1;
  }
  function appendLocusColumn(locus, ci) {
    var role = ci < locus.n_left ? "left flank (anchor)"
      : (ci >= locus.n_left + locus.n_locus ? "right flank (anchor)" : "locus gene");
    tipEl.appendChild(el("div", "tip-id", locus.families[ci]));
    tipEl.appendChild(el("div", null, "Column " + (ci + 1) + " of " + locus.families.length + ", " + role));
    tipEl.appendChild(el("div", null, "Frequency bin: " + ((locus.family_bins && locus.family_bins[ci]) || "unknown")));
    appendLocation({ family_locations: locus.family_locations, example_strain: locus.exemplar }, ci);
    tipEl.appendChild(el("div", null, classLabel((locus.family_classes && locus.family_classes[ci]) || "unannotated")));
    var doms = (locus.family_domains && locus.family_domains[ci]) || "";
    tipEl.appendChild(el("div", null, doms ? "Domains: " + doms : "No annotated Pfam domain"));
  }
  lHead.addEventListener("mousemove", function (e) {
    var locus = LOCI[lstate.selected];
    if (!locus) return;
    var ci = lColAt(lHead, e.clientX, locus);
    if (ci < 0) { tipEl.style.display = "none"; return; }
    tipEl.textContent = "";
    appendLocusColumn(locus, ci);
    positionTip(e);
  });
  lHead.addEventListener("mouseleave", function () { tipEl.style.display = "none"; });
  lGrid.addEventListener("mousemove", function (e) {
    var locus = LOCI[lstate.selected];
    var row = lRows[Math.floor((e.clientY - lGrid.getBoundingClientRect().top) / L_ROW_H)];
    if (!locus || !row) { tipEl.style.display = "none"; return; }
    tipEl.textContent = "";
    tipEl.appendChild(el("div", "tip-id", row.count + (row.count === 1 ? " strain" : " strains") +
      " · " + locusClassLabel(row.row_class)));
    strainPopupLines(row).slice(1).forEach(function (line) { tipEl.appendChild(el("div", null, line)); });
    if (Object.keys(SPECIES).length) {
      var counts = speciesCounts(row, SPECIES);
      Object.keys(counts).sort().forEach(function (sp) {
        tipEl.appendChild(el("div", null, (sp || "Unknown species") + ": " + counts[sp]));
      });
    }
    var ci = lColAt(lGrid, e.clientX, locus);
    if (ci >= 0) {
      tipEl.appendChild(el("div", "tip-id", locus.families[ci]));
      cellReasonLines(locus, row, ci, LPARAMS.k).forEach(function (line) {
        tipEl.appendChild(el("div", null, line));
      });
    }
    positionTip(e);
  });
  lGrid.addEventListener("mouseleave", function () { tipEl.style.display = "none"; });

  function setView(view) {
    var loci = view === "loci";
    document.getElementById("locus-view").classList.toggle("hidden", !loci);
    document.getElementById("island-view").classList.toggle("hidden", loci);
    document.getElementById("lv-btn-loci").setAttribute("aria-pressed", loci ? "true" : "false");
    document.getElementById("lv-btn-islands").setAttribute("aria-pressed", loci ? "false" : "true");
    if (loci) renderLocusMain();
    else if (state.selected >= 0) renderMain();
  }
  document.getElementById("lv-btn-loci").addEventListener("click", function () { setView("loci"); });
  document.getElementById("lv-btn-islands").addEventListener("click", function () { setView("islands"); });
  document.getElementById("lv-search").addEventListener("input", function (e) {
    lstate.search = e.target.value;
    applyLocusFilter();
    renderLocusSidebar();
  });
  document.getElementById("lv-sort").addEventListener("change", function (e) {
    lstate.sort = e.target.value;
    applyLocusFilter();
    renderLocusSidebar();
  });
  var islandSkinChange = window.onSkinChange;
  window.onSkinChange = function () {
    if (islandSkinChange) islandSkinChange();
    if (LOCI.length) renderLocusMain();
  };
  if (DATA.assembly_confound) {
    var conf = document.getElementById("lv-confound");
    conf.textContent = "Accessory content is confounded with assembly quality on this run " +
      "(assembly_quality_confound triggered), so uninformative rows are expected.";
    conf.classList.remove("hidden");
  }
  if (LOCI.length) {
    document.getElementById("lv-switch").classList.remove("hidden");
    applyLocusFilter();
    renderLocusSidebar();
    setView("loci");
  }
"""
```

In `lib/island_synteny_template.py`, replace this text (it occurs once):

```python
from skins import SKINS
```

with:

```python
from island_locus_template import LOCUS_VIEW_CSS, LOCUS_VIEW_HTML, LOCUS_VIEW_JS
from skins import SKINS
```

In `lib/island_synteny_template.py`, replace this text (it occurs once):

```python
  .isv-empty-state h3 { margin: 0 0 8px; color: var(--text-primary); font-size: 15px; }
```

with:

```python
  .isv-empty-state h3 { margin: 0 0 8px; color: var(--text-primary); font-size: 15px; }
""" + LOCUS_VIEW_CSS + r"""
```

In `lib/island_synteny_template.py`, replace this text (it occurs once):

```python
  </header>

  <section class="card">
    <h2 class="card-title">Accessory island synteny</h2>
```

with:

```python
  </header>
""" + LOCUS_VIEW_HTML + r"""
  <div id="island-view">
  <section class="card">
    <h2 class="card-title">Accessory island synteny</h2>
```

In `lib/island_synteny_template.py`, replace this text (it occurs once):

```python
      <p id="isv-empty-text"></p>
    </div>
  </div>
```

with:

```python
      <p id="isv-empty-text"></p>
    </div>
  </div>
  </div>
```

In `lib/island_synteny_template.py`, replace this text (it occurs once):

```python
""" + SKIN_PICKER_JS + r"""
```

with:

```python
""" + LOCUS_VIEW_JS + r"""

""" + SKIN_PICKER_JS + r"""
```


- [ ] **Step 4: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus_template.py tests/test_report_templates.py tests/test_island_synteny.py tests/test_island_synteny_species_sort.py
```

Expected: all pass (10 new tests; the existing template tests still pass, including `test_page_javascript_parses[island_synteny]`, which runs `node --check`)


- [ ] **Step 5: Commit**

```bash
git add lib/island_locus_template.py lib/island_synteny_template.py tests/test_island_locus_template.py
git commit -m "island locus view: locus view in island_synteny.html, default when loci exist (Ruling R1) (#${ISSUE_A})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 10: `bin/pangenome_island_synteny.py`: `--loci_json` and `--diagnostics_tsv`

**Files:**
- Modify: `bin/pangenome_island_synteny.py`
- Modify: `tests/test_pangenome_island_synteny.py` (append)

**Interfaces:**
- Consumes: `island_loci.json` (Task 7), `diagnostics/diagnostics.tsv` (DIAGNOSTICS).
- Produces: `assembly_confound_triggered(path) -> bool`; `add_locus_view(payload, loci_json, diagnostics_tsv) -> None` which sets `payload['assembly_confound']` always and `payload['loci']`, `payload['locus_meta']` (everything in the JSON except `loci`) when the JSON is non-empty.


- [ ] **Step 1: Write the failing test**

Append to the end of `tests/test_pangenome_island_synteny.py`:

```python


# ---- locus view (docs/superpowers/plans/2026-09-26-island-locus-view-clinker.md) ----

LOCI_JSON = {
    "project": "demo",
    "locus_params": {"flank": 5, "flank_min": 3, "k": 10},
    "n_loci_total": 1,
    "loci": [{"key": "L001", "locus_id": "S1:c1:100-400", "families": ["famA", "famB"],
              "rows": [], "counts": {"full": 1, "partial": 0, "empty": 0, "uninformative": 1}}],
}


def test_cli_without_loci_json_has_no_loci_key(tmp_path):
    payload = payload_of(run_cli(tmp_path))
    assert "loci" not in payload
    assert payload["assembly_confound"] is False


def test_cli_with_loci_json_embeds_loci_and_meta(tmp_path):
    lj = tmp_path / "island_loci.json"
    lj.write_text(json.dumps(LOCI_JSON))
    payload = payload_of(run_cli(tmp_path, "--loci_json", str(lj)))
    assert payload["loci"][0]["key"] == "L001"
    assert payload["locus_meta"]["locus_params"]["k"] == 10
    assert "loci" not in payload["locus_meta"]


def test_cli_empty_loci_json_is_ignored(tmp_path):
    lj = tmp_path / "island_loci.json"
    lj.write_text("")
    assert "loci" not in payload_of(run_cli(tmp_path, "--loci_json", str(lj)))


def test_cli_reads_the_assembly_confound_diagnostic(tmp_path):
    diag = tmp_path / "diagnostics.tsv"
    diag.write_text("diagnostic_id\tstatus\twould_fail_strict\tdetail\tprune_options\n"
                    "assembly_quality_confound\ttriggered\tTrue\tmax rho\t\n")
    assert payload_of(run_cli(tmp_path, "--diagnostics_tsv", str(diag)))["assembly_confound"] is True
    diag.write_text("diagnostic_id\tstatus\twould_fail_strict\tdetail\tprune_options\n"
                    "assembly_quality_confound\tok\tFalse\tmax rho\t\n")
    assert payload_of(run_cli(tmp_path, "--diagnostics_tsv", str(diag)))["assembly_confound"] is False
```


- [ ] **Step 2: Run**

```bash
pixi run python -m pytest -q tests/test_pangenome_island_synteny.py
```

Expected: FAIL: 4 new tests fail (`KeyError: 'assembly_confound'` and `error: unrecognized arguments: --loci_json`)


- [ ] **Step 3: Make the change**

In `bin/pangenome_island_synteny.py`, replace this text (it occurs once):

```python
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
```

with:

```python
    ap.add_argument("--loci_json", default=None,
                    help="island_loci.json from bin/pangenome_island_loci.py; adds the "
                    "locus view (default view of the page). An empty file is ignored.")
    ap.add_argument("--diagnostics_tsv", default=None,
                    help="diagnostics.tsv from bin/pangenome_diagnostics.py; when "
                    "assembly_quality_confound is 'triggered' the page says so.")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
```

In `bin/pangenome_island_synteny.py`, replace this text (it occurs once):

```python
    # Escape `</` so a Pfam description
```

with:

```python
    add_locus_view(payload, args.loci_json, args.diagnostics_tsv)

    # Escape `</` so a Pfam description
```

In `bin/pangenome_island_synteny.py`, replace this text (it occurs once):

```python
def main() -> int:
```

with:

```python
def assembly_confound_triggered(path: str | None) -> bool:
    """True when diagnostics.tsv has assembly_quality_confound = triggered."""
    if not path or not Path(path).is_file() or Path(path).stat().st_size == 0:
        return False
    with open_maybe_compressed(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("diagnostic_id") == "assembly_quality_confound":
                return row.get("status") == "triggered"
    return False


def add_locus_view(payload: dict, loci_json: str | None, diagnostics_tsv: str | None) -> None:
    """Add `loci`, `locus_meta` and `assembly_confound` to the page payload.
    Without --loci_json (or with an empty file) the payload has no `loci`
    and the page shows the island view only, as before."""
    payload["assembly_confound"] = assembly_confound_triggered(diagnostics_tsv)
    if not loci_json or not Path(loci_json).is_file() or Path(loci_json).stat().st_size == 0:
        return
    data = json.loads(Path(loci_json).read_text())
    payload["loci"] = data.pop("loci", [])
    payload["locus_meta"] = data


def main() -> int:
```


- [ ] **Step 4: Run**

```bash
pixi run python -m pytest -q tests/test_pangenome_island_synteny.py
```

Expected: 33 passed


- [ ] **Step 5: Commit**

```bash
git add bin/pangenome_island_synteny.py tests/test_pangenome_island_synteny.py
git commit -m "island locus view: embed island_loci.json and the confound flag in the page (#${ISSUE_A})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 11: Browser behaviour under jsdom

**Files:**
- Modify: `tests/test_report_js_behaviour.py`
- Modify: `tests/js/drive_reports.mjs`

**Interfaces:**
- Consumes: the page from Tasks 9-10.
- Produces: fixture page `island_synteny_loci.html` and 15 `locus view:` checks.


- [ ] **Step 1: Read this first**

This suite skips unless jsdom resolves. Install it outside the repo first:

```bash
npm install --prefix "$SCRATCH/novinvenio-jsdom" jsdom --no-audit --no-fund --silent
export NOVINVENIO_JSDOM="$SCRATCH/novinvenio-jsdom/node_modules/jsdom"
```

On a login node without `$SCRATCH`, use the session scratchpad directory instead.


- [ ] **Step 2: Add the fixture page and the checks**

In `tests/test_report_js_behaviour.py`, replace this text (it occurs once):

```python
def _fasta() -> str:
```

with:

```python
# Locus view (docs/superpowers/plans/2026-09-26-island-locus-view-clinker.md):
# a hand-written island_loci.json with two loci. L001 ranks first by
# informative polymorphism, L002 first by locus ID ("A:..." < "B:...").
def _locus(key, locus_id, score, exemplar):
    return {
        "key": key, "locus_id": locus_id, "exemplar": exemplar, "exemplar_contig": "c1",
        "exemplar_span": {"start": 10, "end": 90}, "tier": "full", "size": 2,
        "n_variants": 1, "variant_strains": 3, "families": ["F1", "famA", "famB", "G1"],
        "family_bins": ["core", "shell", "shell", "core"],
        "family_classes": ["unannotated"] * 4, "family_domains": [""] * 4,
        "dominant_class": "unannotated", "n_left": 1, "n_locus": 2, "n_right": 1,
        "counts": {"full": 2, "partial": 0, "empty": 1, "uninformative": 0},
        "counts_by_species": {"": {"full": 2, "partial": 0, "empty": 1, "uninformative": 0}},
        "n_carriers": 2, "informative_score": score,
        "rows": [
            {"row_class": "full", "codes": "1111", "count": 2, "strains": ["S1", "S2"],
             "rep": {"c": ["c1"], "d": [[0, 1, 1, 1], [0, 2, 0, -1], [0, 3, 1, -1], [0, 4, 2, -1]]}},
            {"row_class": "empty", "codes": "1001", "count": 1, "strains": ["S3"], "rep": None},
        ],
        "breakpoints": [{"b": 1, "indel": {"": 1}, "contig_break": 0},
                        {"b": 3, "indel": {"": 1}, "contig_break": 0}],
        "clinker_strains": [{"strain": exemplar, "reason": "exemplar", "row_class": "full",
                             "species": "", "contig": "c1", "rank_lo": 0, "rank_hi": 5}],
    }


ISV_LOCI = {
    "project": "demo",
    "locus_params": {"flank": 5, "flank_min": 3, "k": 10, "empty_frac": 0.8},
    "n_loci_total": 2,
    "loci": [_locus("L001", "B:c1:1-9", 7, "S1"), _locus("L002", "A:c1:1-9", 1, "S2")],
}


def _fasta() -> str:
```

In `tests/test_report_js_behaviour.py`, replace this text (it occurs once):

```python
    (d / 'isv_empty_islands.tsv').write_text(ISV_EMPTY_ISLANDS)
```

with:

```python
    (d / 'isv_loci.json').write_text(json.dumps(ISV_LOCI))
    run('pangenome_island_synteny.py',
        '--islands_with_domains', str(d / 'isv_islands.tsv'),
        '--presence_matrix', str(d / 'isv_matrix.tsv'),
        '--family_positions', str(d / 'isv_positions.tsv'),
        '--loci_json', str(d / 'isv_loci.json'),
        '--project', 'demo', '--output', str(d / 'island_synteny_loci.html'))

    (d / 'isv_empty_islands.tsv').write_text(ISV_EMPTY_ISLANDS)
```

In `tests/js/drive_reports.mjs`, replace this text (it occurs once):

```javascript
console.log(failures === 0 ? 'ALL PASSED' : failures + ' FAILED');
```

with:

```javascript
// ------------------------------------------------ island synteny: locus view
{
  const dom = boot(path.join(FX, 'island_synteny_loci.html'));
  const w = dom.window, d = w.document;
  const errors = [];
  w.addEventListener('error', (e) => errors.push(String(e.error)));
  await sleep(60);
  check('locus view: loads without error', errors.length === 0, errors.join('; '));
  check('locus view: is the default view',
        !d.getElementById('locus-view').classList.contains('hidden') &&
        d.getElementById('island-view').classList.contains('hidden'));
  check('locus view: view switch is shown',
        !d.getElementById('lv-switch').classList.contains('hidden'));
  const items = () => [...d.querySelectorAll('#lv-list .isv-item')];
  check('locus view: sidebar lists both loci', items().length === 2, items().length);
  check('locus view: default order is informative polymorphism (L001 first)',
        items()[0].textContent.includes('B:c1:1-9'), items()[0].textContent);
  check('locus view: title is the exemplar locus',
        d.getElementById('lv-title').textContent === 'S1:c1:10-90',
        d.getElementById('lv-title').textContent);
  items()[1].dispatchEvent(ev(w, 'click'));
  check('locus view: selecting the second locus updates the title',
        d.getElementById('lv-title').textContent === 'S2:c1:10-90',
        d.getElementById('lv-title').textContent);
  const sort = d.getElementById('lv-sort');
  sort.value = 'name';
  sort.dispatchEvent(ev(w, 'change'));
  check('locus view: sort by locus ID puts A:... first',
        items()[0].textContent.includes('A:c1:1-9'), items()[0].textContent);
  const search = d.getElementById('lv-search');
  search.value = 'B:c1';
  search.dispatchEvent(ev(w, 'input'));
  check('locus view: search filters the list', items().length === 1, items().length);
  check('locus view: legend has the five cell states',
        d.querySelectorAll('#lv-legend .isv-swatch').length === 5,
        d.querySelectorAll('#lv-legend .isv-swatch').length);
  const drawn = w.__fillTextCalls.join('|');
  check('locus view: grid rows were drawn', drawn.includes('S1 +1') && drawn.includes('S3'), drawn.slice(0, 200));
  d.getElementById('lv-btn-islands').dispatchEvent(ev(w, 'click'));
  check('locus view: islands button shows the island view',
        d.getElementById('locus-view').classList.contains('hidden') &&
        !d.getElementById('island-view').classList.contains('hidden'));
  check('locus view: islands button is pressed',
        d.getElementById('lv-btn-islands').getAttribute('aria-pressed') === 'true');
}

// ------------------------------- island synteny: no loci keeps the old page
{
  const dom = boot(path.join(FX, 'island_synteny.html'));
  const d = dom.window.document;
  await sleep(60);
  check('locus view: absent without loci (switch hidden)',
        d.getElementById('lv-switch').classList.contains('hidden'));
  check('locus view: absent without loci (island view shown)',
        !d.getElementById('island-view').classList.contains('hidden'));
}

console.log(failures === 0 ? 'ALL PASSED' : failures + ' FAILED');
```


- [ ] **Step 3: Run**

```bash
pixi run python -m pytest -q tests/test_report_js_behaviour.py
```

Expected: 2 passed. To see the new checks: rerun the driver by hand on the fixture dir the test writes (`--basetemp="$SCRATCH/jsfx"`, then `node tests/js/drive_reports.mjs "$SCRATCH"/jsfx/reports0 "$NOVINVENIO_JSDOM" | grep 'locus view'`): 15 `PASS locus view:` lines, then `ALL PASSED`.


- [ ] **Step 4: Read this first**

This driver passed during planning (all 15 locus-view checks). A check fails if the default view, the sidebar order, selection, sort, search, the five-state legend or the view switch is broken.


- [ ] **Step 5: Commit**

```bash
git add tests/test_report_js_behaviour.py tests/js/drive_reports.mjs
git commit -m "island locus view: jsdom checks for the locus view (#${ISSUE_A})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 12: Nextflow: `ISLAND_LOCI`, params, help, docs

**Files:**
- Create: `modules/pangenome/island_loci.nf`
- Modify: `modules/pangenome/island_synteny.nf`
- Modify: `workflows/pangenome_profile.nf`
- Modify: `nextflow.config` (params block only)
- Modify: `pangenome.nf` (help text, param validation)
- Modify: `docs/pangenome-assumptions.md`
- Modify: `CHANGES.md`
- Modify: `.living/decisions.md`

**Interfaces:**
- Consumes: `REPORT_TABLES.out.islands_with_domains`, `rescued_matrix`, `FAMILY_POSITIONS.out.positions`, `FREQUENCY_BINS.out.table`, `ASSEMBLY_QUALITY_QC.out.table`, `samplesheet`, `FAMILY_PFAM_SCAN.out.domtblout`, `GENE_POSITIONS.out.positions`, `CLUSTER_TIER1.out.cluster_tsv`, `DIAGNOSTICS.out.tsv`.
- Produces: process `ISLAND_LOCI` with `emit: loci`; `ISLAND_SYNTENY` gains inputs `path(loci_json)`, `path(diagnostics_tsv)`; params `pangenome_locus_flank=5`, `pangenome_locus_flank_min=3`, `pangenome_locus_k=10`, `pangenome_locus_empty_frac=0.8`, `pangenome_locus_containment=0.5`, `pangenome_locus_rank='informative'`, `pangenome_top_loci=50`, `pangenome_locus_candidates=200`.


- [ ] **Step 1: Read this first**

There is no pangenome integration test in this repo. The checks are `nextflow lint` (Nextflow 26 strict parser) and `--help`, which compiles every included module. Do not touch `lib/Helpers.groovy` or the `process {}` block of `nextflow.config`: PR branch `fix-194-189-project-name-presence-matrix-mem` changes those. This task only adds lines at the end of the `params {}` block.


- [ ] **Step 2: Make the change**

Create `modules/pangenome/island_loci.nf`:

```groovy
// ISLAND_LOCI -- island locus view (spec
// docs/superpowers/specs/2026-09-24-island-locus-view-design.md): groups
// islands into loci, picks each locus's exemplar, computes every strain's
// cell states, row class and the breakpoint track, and writes
// island_loci.json for ISLAND_SYNTENY.
//
// Streams family_positions three times, filtered each time (see
// bin/pangenome_island_loci.py). Measured on the 529-strain Coccidioides run
// (2026-09-26, 200 candidate loci, 50 drawn): 1 min 41 s wall, 1.36 GB peak
// RSS, 1.94 MB island_loci.json -- inside low_cpu's 4 GB.
//
// --assembly_quality is ASSEMBLY_QUALITY_QC's table; it covers the ingroup
// only, so outgroup strains fall back to locus-contig gene count for ties
// (lib/island_locus.py::quality_key, plan Ruling R4).
process ISLAND_LOCI {
    label 'low_cpu'
    tag "island_loci"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(islands_with_domains)
    path(presence_matrix)
    path(family_positions)
    path(frequency_table)
    path(assembly_quality)
    path(samplesheet)
    path(domtblout)
    path(gene_positions)
    path(cluster_tsv)

    output:
    path("island_loci.json"), emit: loci

    script:
    """
    pangenome_island_loci.py \
        --islands_with_domains ${islands_with_domains} \
        --presence_matrix ${presence_matrix} \
        --family_positions ${family_positions} \
        --frequency_table ${frequency_table} \
        --assembly_quality ${assembly_quality} \
        --config ${samplesheet} \
        --domtblout ${domtblout} \
        --domain_evalue ${params.pangenome_pfam_domain_evalue} \
        --gene_positions ${gene_positions} \
        --cluster_tsv ${cluster_tsv} \
        --id_sep '${params.pangenome_id_sep}' \
        --flank ${params.pangenome_locus_flank} \
        --flank_min ${params.pangenome_locus_flank_min} \
        --k ${params.pangenome_locus_k} \
        --empty_frac ${params.pangenome_locus_empty_frac} \
        --containment ${params.pangenome_locus_containment} \
        --rank_by ${params.pangenome_locus_rank} \
        --top_loci ${params.pangenome_top_loci} \
        --candidates ${params.pangenome_locus_candidates} \
        --min_strains ${params.pangenome_top_islands_min_strains} \
        --project '${Helpers.projectName(params)}' \
        --output island_loci.json
    """
}
```

In `modules/pangenome/island_synteny.nf`, replace this text (it occurs once):

```groovy
    path(rescue_positions)

    output:
```

with:

```groovy
    path(rescue_positions)
    path(loci_json)
    path(diagnostics_tsv)

    output:
```

In `modules/pangenome/island_synteny.nf`, replace this text (it occurs once):

```groovy
        --id_sep '${params.pangenome_id_sep}' \
        --output island_synteny.html
```

with:

```groovy
        --id_sep '${params.pangenome_id_sep}' \
        --loci_json ${loci_json} \
        --diagnostics_tsv ${diagnostics_tsv} \
        --output island_synteny.html
```

In `workflows/pangenome_profile.nf`, replace this text (it occurs once):

```groovy
include { ISLAND_SYNTENY }                                                 from '../modules/pangenome/island_synteny'
```

with:

```groovy
include { ISLAND_SYNTENY }                                                 from '../modules/pangenome/island_synteny'
include { ISLAND_LOCI }                                                    from '../modules/pangenome/island_loci'
```

In `workflows/pangenome_profile.nf`, replace this text (it occurs once):

```groovy
    if (params.pangenome_island_pfam_hmm) {
        ISLAND_SYNTENY(
            REPORT_TABLES.out.islands_with_domains,
            rescued_matrix,
            FAMILY_POSITIONS.out.positions,
            FAMILY_PFAM_SCAN.out.domtblout,
            samplesheet,
            DIAGNOSTICS.out.banner_html,
            GENE_POSITIONS.out.positions,
            CLUSTER_TIER1.out.cluster_tsv,
            rescue_positions,
        )
    }
```

with:

```groovy
    if (params.pangenome_island_pfam_hmm) {
        // Island locus view (docs/superpowers/specs/2026-09-24-island-locus-view-design.md).
        ISLAND_LOCI(
            REPORT_TABLES.out.islands_with_domains,
            rescued_matrix,
            FAMILY_POSITIONS.out.positions,
            FREQUENCY_BINS.out.table,
            ASSEMBLY_QUALITY_QC.out.table,
            samplesheet,
            FAMILY_PFAM_SCAN.out.domtblout,
            GENE_POSITIONS.out.positions,
            CLUSTER_TIER1.out.cluster_tsv,
        )
        ISLAND_SYNTENY(
            REPORT_TABLES.out.islands_with_domains,
            rescued_matrix,
            FAMILY_POSITIONS.out.positions,
            FAMILY_PFAM_SCAN.out.domtblout,
            samplesheet,
            DIAGNOSTICS.out.banner_html,
            GENE_POSITIONS.out.positions,
            CLUSTER_TIER1.out.cluster_tsv,
            rescue_positions,
            ISLAND_LOCI.out.loci,
            DIAGNOSTICS.out.tsv,
        )
    }
```

In `nextflow.config`, replace this text (it occurs once):

```groovy
    pangenome_viz_top_islands = 50
}
```

with:

```groovy
    pangenome_viz_top_islands = 50

    // Island locus view (docs/superpowers/specs/2026-09-24-island-locus-view-design.md,
    // ISLAND_LOCI). Defaults are the spec's, taken from the 2026-09-24
    // feasibility run and not swept (spec open question 1).
    pangenome_locus_flank       = 5              // F: flank genes per side
    pangenome_locus_flank_min   = 3              // F_min: fallback when no carrier has F per side
    pangenome_locus_k           = 10             // k: gene window of the "in place" test
    pangenome_locus_empty_frac  = 0.8            // empty site: >= this fraction of locus columns absent
    pangenome_locus_containment = 0.5            // an island joins a larger island's locus at >= this
    pangenome_locus_rank        = 'informative'  // informative | strains | size
    pangenome_top_loci          = 50             // loci drawn in island_synteny.html
    pangenome_locus_candidates  = 200            // loci whose states are computed before ranking
}
```

In `pangenome.nf`, replace this text (it occurs once):

```groovy
      --pangenome_viz_top_islands      Islands drawn in island_synteny.html
                                       (default: 50), selected by size after
                                       the min-strains filter above.
```

with:

```groovy
      --pangenome_viz_top_islands      Islands drawn in island_synteny.html
                                       (default: 50), selected by size after
                                       the min-strains filter above.
      --pangenome_top_loci             Loci drawn in island_synteny.html's locus
                                       view (default: 50).
      --pangenome_locus_rank           Locus ranking: informative (default; >= 10
                                       empty-site and >= 2 full-locus strains,
                                       by the smaller count), strains, or size.
      --pangenome_locus_flank          Flank genes per side (default: 5); falls back
                                       to --pangenome_locus_flank_min (default: 3).
      --pangenome_locus_k              Gene window of the "in place" test (default: 10).
      --pangenome_locus_empty_frac     Empty site: this fraction of locus columns
                                       absent (default: 0.8).
      --pangenome_locus_containment    An island joins a larger island's locus when
                                       this fraction of its families is in it (default: 0.5).
      --pangenome_locus_candidates     Loci scored before ranking (default: 200).
```

In `pangenome.nf`, replace this text (it occurs once):

```groovy
    if (params.pangenome_cluster_backend !in ['mmseqs', 'diamond'])
```

with:

```groovy
    if (params.pangenome_locus_rank !in ['informative', 'strains', 'size'])
        error "ERROR: --pangenome_locus_rank must be informative, strains or size (got: ${params.pangenome_locus_rank})"
    if (params.pangenome_locus_flank_min > params.pangenome_locus_flank)
        error "ERROR: --pangenome_locus_flank_min (${params.pangenome_locus_flank_min}) must not exceed --pangenome_locus_flank (${params.pangenome_locus_flank})"
    if (params.pangenome_cluster_backend !in ['mmseqs', 'diamond'])
```


- [ ] **Step 3: Run**

```bash
pixi run nextflow lint modules/pangenome/island_loci.nf modules/pangenome/island_synteny.nf workflows/pangenome_profile.nf pangenome.nf
```

Expected: `Nextflow linting complete!` with no errors (warnings that already exist on main are fine)


- [ ] **Step 4: Run**

```bash
pixi run nextflow run pangenome.nf --help | grep -E 'pangenome_top_loci|pangenome_locus_rank'
```

Expected: both help lines print; no `Script compilation failed`


- [ ] **Step 5: Run**

```bash
pixi run nextflow run pangenome.nf --pangenome_samplesheet tests/data/test.csv --pangenome_data_dir tests/data --pangenome_locus_rank bogus --outdir "$SCRATCH/lv_badparam" 2>&1 | tail -3
```

Expected: `ERROR: --pangenome_locus_rank must be informative, strains or size (got: bogus)` before any task runs (outdir is outside the repo, per CLAUDE.md)


- [ ] **Step 6: Make the change**

In `docs/pangenome-assumptions.md`, replace this text (it occurs once):

```markdown
| `pangenome_viz_top_islands` | `50` | `theory+practical` | Payload-size constraint, measured directly: at top-50, the embedded `island_synteny.html` payload is 0.55MB (0.28MB after a proposed `strain_hap` index-array optimization not yet implemented); top-50 references 298 of 54,421 distinct families. Browsers parse a 5-10MB inline JSON block well under a second from `file://`; the current encoding does not cross 10MB below ~1,500 strains at top-500, or ~4,000 strains at top-500 with the `strain_hap` optimization. | 529-strain Coccidioides, real payload rebuild from `lib/island_synteny.build_payload` | A much larger strain count (thousands) would need this revisited alongside the `strain_hap` encoding change; both are quantified in the storage-and-compute note. | Rebuild the payload at the target `top_islands` value and measure its size directly (`lib/island_synteny.build_payload`) before assuming 50 is still small enough. |
```

with:

```markdown
| `pangenome_viz_top_islands` | `50` | `theory+practical` | Payload-size constraint, measured directly: at top-50, the embedded `island_synteny.html` payload is 0.55MB (0.28MB after a proposed `strain_hap` index-array optimization not yet implemented); top-50 references 298 of 54,421 distinct families. Browsers parse a 5-10MB inline JSON block well under a second from `file://`; the current encoding does not cross 10MB below ~1,500 strains at top-500, or ~4,000 strains at top-500 with the `strain_hap` optimization. | 529-strain Coccidioides, real payload rebuild from `lib/island_synteny.build_payload` | A much larger strain count (thousands) would need this revisited alongside the `strain_hap` encoding change; both are quantified in the storage-and-compute note. | Rebuild the payload at the target `top_islands` value and measure its size directly (`lib/island_synteny.build_payload`) before assuming 50 is still small enough. |
| `pangenome_locus_flank` | `5` | `unvalidated-assumption` | Spec default F from the 2026-09-24 feasibility run: with F = 5 the exemplar has 5 genes on both sides for 30/50 loci (top 50 by strain count); not swept. | 529-strain Coccidioides `rescue_freqpol_immitis_in_posadasii_out` | Short contigs everywhere (most loci fall to the F_min tier) | Count loci per tier in `island_loci.json` (`tier`) on the new study |
| `pangenome_locus_flank_min` | `3` | `unvalidated-assumption` | Reviewer decision 2026-09-24 (spec section 2): keeps a flank anchor near contig ends. 1 of 50 drawn loci used it on the Coccidioides run (2026-09-26). | same | — | as above |
| `pangenome_locus_k` | `10` | `unvalidated-assumption` | Same window as `pangenome_pair_class_k`; not swept (spec open question 1). | same | Gene-dense islands where unrelated genes sit within 10 genes | Recompute states at k = 5 and 20 and compare row-class counts |
| `pangenome_locus_empty_frac` | `0.8` | `unvalidated-assumption` | Feasibility: with the top 50 by size, 17/50 loci qualify at 80% and 0/50 at 100%; top 50 by strain count 39/50 at both. | same | — | Compare qualifying loci at 0.8 and 1.0 |
| `pangenome_locus_containment` | `0.5` | `unvalidated-assumption` | The M4 criterion: 10,579 of 11,880 located islands join a larger island at 50%. Not swept. | same | — | Count loci at 0.3/0.5/0.7 |
| `pangenome_locus_rank` / `pangenome_top_loci` | `informative` / `50` | `theory+practical` | Spec section 6. Caution: on the Coccidioides run 0 of 6 sequence-checked empty-site calls on the top 3 loci were real deletions (gene-model splits; plan Task 14). | same | Gene-model differences dominate the informative ranking | Spot check (plan Task 14) on the new study |
| `pangenome_locus_candidates` | `200` | `theory+practical` | Compute bound: states for 200 candidates took 1 min 41 s and 1.36 GB (2026-09-26); 1301 loci on that run. | same | Informative loci outside the top 200 by carrier proxy | Rerun with 400 and compare the drawn set |
```

In `CHANGES.md`, replace this text (it occurs once):

```markdown
## Unreleased

```

with:

```markdown
## Unreleased

### New: island locus view (default view of island_synteny.html)

- **`lib/island_locus.py`**, **`bin/pangenome_island_loci.py`**, **`modules/pangenome/island_loci.nf`**,
  **`lib/island_locus_template.py`** (spec `docs/superpowers/specs/2026-09-24-island-locus-view-design.md`,
  plan `docs/superpowers/plans/2026-09-26-island-locus-view-clinker.md`) -- islands are grouped into
  loci; each locus is drawn on one exemplar's genes (5 flank genes each side, 3 near contig ends) and
  every strain gets a state per column: in place, elsewhere, rescue (hatched), absent or contig break.
  Rows are grouped into full locus / partial / empty site / uninformative, and a track above the
  columns counts shared breakpoints per species. The page opens on this view; the island presence
  view is one click away.
- New params: `--pangenome_top_loci` (50), `--pangenome_locus_rank` (`informative`), `--pangenome_locus_flank` (5),
  `--pangenome_locus_flank_min` (3), `--pangenome_locus_k` (10), `--pangenome_locus_empty_frac` (0.8),
  `--pangenome_locus_containment` (0.5), `--pangenome_locus_candidates` (200).
- Caution: on the Coccidioides run the top loci's "empty site" calls were gene-model splits, not
  deletions, in all 6 sequence-checked cases (plan Task 14).

```

Append to the end of `.living/decisions.md`:

```markdown


## 2026-09-26 — Island locus view: page, ranking, carriers and N50 gaps

**Context**: Spec `docs/superpowers/specs/2026-09-24-island-locus-view-design.md` left the page choice (open question 3) open and assumed N50 for every strain and a list of carrier strains, which the inputs do not provide.
**Decision**: Locus view is the default view inside the existing `island_synteny.html`; the island presence view stays one click away (plan Ruling R1). Carriers are strains with every root member family on one contig within `len - 1 + k` ranks (R3). Strains without an N50 (ASSEMBLY_QUALITY_QC covers the ingroup only: 169 of 529) rank after strains with one, by locus-contig gene count (R4). States are computed for 200 candidate loci before ranking (R5).
**Alternatives**: a separate locus page; compute N50 for all strains in ASSEMBLY_QUALITY_QC (a change to a module outside this spec); rank every locus (1301 loci, memory grows with column families).
**Rationale**: keeps the current page and links working; uses only existing inputs; measured cost 1 min 41 s, 1.36 GB.
```


- [ ] **Step 7: Commit**

```bash
git add modules/pangenome/island_loci.nf modules/pangenome/island_synteny.nf workflows/pangenome_profile.nf nextflow.config pangenome.nf docs/pangenome-assumptions.md CHANGES.md .living/decisions.md
git commit -m "island locus view: ISLAND_LOCI process, params and docs (#${ISSUE_A})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 13: Real-data run and page render (controller task)

**Files:**
- No repo files. Outputs go to `$SCRATCH/lv_real/` (never under the repo).

**Interfaces:**
- Consumes: the real run `rescue_freqpol_immitis_in_posadasii_out` and Tasks 7 and 10.
- Produces: measured time, peak RSS and sizes to paste into the PR description.


- [ ] **Step 1: Read this first**

Run on a compute node (`srun -p short -c 2 --mem 8G --pty bash`), from the worktree root.


- [ ] **Step 2: Run**

```bash
export PATH="$PWD/.pixi/envs/default/bin:$PATH"
ST=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome
P=$ST/results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome
O="${SCRATCH:?}/lv_real"; mkdir -p "$O"
/usr/bin/time -v python bin/pangenome_island_loci.py \
  --islands_with_domains $P/report_tables/islands_with_domains.tsv \
  --presence_matrix $P/presence_matrix.rescued.tsv --family_positions $P/family_positions.tsv.zst \
  --frequency_table $P/frequency_table.tsv --assembly_quality $P/assembly_quality_vs_content.tsv \
  --config $ST/config_immitis_in_posadasii_out.csv --domtblout $P/pfam.domtblout \
  --gene_positions $P/gene_positions.tsv.zst --cluster_tsv $P/cluster/tier1_cluster.tsv \
  --project cocci --output "$O/island_loci.json" 2> "$O/loci.time"
grep -E 'pangenome_island_loci|Elapsed|Maximum' "$O/loci.time"
python bin/pangenome_island_synteny.py \
  --islands_with_domains $P/report_tables/islands_with_domains.tsv \
  --presence_matrix $P/presence_matrix.rescued.tsv --family_positions $P/family_positions.tsv.zst \
  --domtblout $P/pfam.domtblout --config $ST/config_immitis_in_posadasii_out.csv \
  --diagnostics_banner $P/diagnostics/diagnostics_banner.html \
  --diagnostics_tsv $P/diagnostics/diagnostics.tsv --loci_json "$O/island_loci.json" \
  --project cocci --output "$O/island_synteny.html"
ls -la "$O"
```


- [ ] **Step 3: Read this first**

Expected, from the planning run of this exact code (2026-09-26): `1301 loci, 200 candidates, 0 without a carrier, 50 drawn (--rank_by informative)`; about 1 min 41 s wall and 1.36 GB peak RSS; `island_loci.json` about 1.94 MB; `island_synteny.html` about 2.7 MB. Tiers: 48 `full`, 1 `short_flanks`, 1 `contig_end`. First locus `L001` = `1M0:scaffold_217:2489-4804`, exemplar `UTAH_20380X16`, counts full 221, partial 12, empty 292, uninformative 4. Record the actual values; a difference above 10% in time or memory, or any difference in the counts, needs an explanation before the PR.


- [ ] **Step 4: Read this first**

Render check in a real browser engine (headless Chromium from the Playwright cache on the cluster):


- [ ] **Step 5: Run**

```bash
H=$(ls ~/.cache/ms-playwright/chromium_headless_shell-*/chrome-headless-shell-linux64/chrome-headless-shell | tail -1)
"$H" --headless --no-sandbox --disable-gpu --virtual-time-budget=5000 --window-size=1600,1400 \
  --screenshot="$O/page.png" "$O/island_synteny.html"
"$H" --headless --no-sandbox --disable-gpu --virtual-time-budget=5000 --dump-dom "$O/island_synteny.html" \
  | grep -oE 'id="lv-count"[^>]*>[^<]*|id="lv-title">[^<]*|id="island-view" class="[^"]*"'
```


- [ ] **Step 6: Read this first**

Expected: `50 of 50 shown`, the L001 exemplar title (`UTAH_20380X16:scaffold_5:1081614-1083929` in planning) and `island-view` with class `hidden`. Open `page.png`: the breakpoint bars, the rotated column labels with anchor bars under the flank columns, and the grid with five state colours must be visible. Then open the page in a desktop browser from `file://` and hover a column and a cell; the popups must show the family, bin, domains, exemplar location, and the cell's state, contig, rank and neighbour.


### Task 14: Sequence spot checks of empty-site and partial calls (spec validation item 3, controller task)

**Files:**
- Create (outside the repo): `$SCRATCH/lv_real/spotcheck_locus_calls.py` (throwaway; never commit)

**Interfaces:**
- Consumes: `$SCRATCH/lv_real/island_loci.json` from Task 13; the study `data_dir`, `gene_positions.tsv.zst`, tier-1 cluster TSV.
- Produces: a per-column agree/disagree table for 5 empty-site and 5 partial strains.


- [ ] **Step 1: Read this first**

This test checks the neighbour rule against the DNA, which the rule itself cannot do. **Planning result (2026-09-26): 0 of 6 empty-site calls on L001-L003 were confirmed.** In every case the empty-site strain's region DNA aligned across the reference strain's locus genes at 97.0-99.9% identity (blastn). At L001 the reference `UTAH_20380X16` has two locus genes (277 aa and 446 aa, families `CA25|C4385CC_002993-T1` and `UCSF_5AZ|C4980C8_004519-T1`); the empty-site strain `UTAH_20380X10` has one 705 aa gene model over the same DNA, in a third family. So on these loci an "empty site" is a gene-model split or merge, not a deletion. Repeat the check on the loci the final code draws and report the result to the user at the Checkpoint.


- [ ] **Step 2: Make the change**

Create `$SCRATCH/lv_real/spotcheck_locus_calls.py`:

```python
#!/usr/bin/env python3
"""Sequence spot check of locus-view calls (spec validation plan item 3).

Throwaway controller tool: keep it in $SCRATCH, never commit it.

For one locus of island_loci.json, a reference strain (a full-locus strain)
and a query strain: take each strain's region on the contig holding most of
its flank-family genes (bp span of those genes, from gene_positions + the
tier-1 cluster TSV), extract the DNA from the study data_dir, align query to
reference with blastn, and for every locus column report whether the
reference gene is covered (>= 80% of its length) by a query HSP. A column
agrees when "covered" matches the query's call (in place/rescue in place =
covered; absent/contig break = not covered).

Usage (NovInvenio pixi env: biopython + blastn):
  python spotcheck_locus_calls.py --loci_json island_loci.json --locus L001 \
      --ref UTAH_20380X16 --query UTAH_20380X10 --config config.csv \
      --data_dir data_dir --gene_positions gene_positions.tsv.zst \
      --cluster_tsv cluster/tier1_cluster.tsv
"""
import argparse
import collections
import json
import subprocess
import sys
import tempfile
from pathlib import Path

NOVINVENIO = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/locus-view-a")
sys.path.insert(0, str(NOVINVENIO / "lib"))
sys.path.insert(0, str(NOVINVENIO / "bin"))
from Bio import SeqIO  # noqa: E402
from compressed_io import open_maybe_compressed  # noqa: E402
from config_parser import parse_config  # noqa: E402
from pangenome_island_synteny import load_gene_locations  # noqa: E402


def region(locs, strain, flank_fams, pad=0):
    by_contig = collections.defaultdict(list)
    for f in flank_fams:
        for _pid, contig, start, end in locs.get((strain, f), []):
            by_contig[contig].append((start, end))
    if not by_contig:
        return None
    contig = max(by_contig, key=lambda c: len(by_contig[c]))
    spans = by_contig[contig]
    return contig, max(1, min(s for s, _ in spans) - pad), max(e for _, e in spans) + pad


def fetch(fasta, contig, start, end):
    with open_maybe_compressed(str(fasta)) as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            if rec.id == contig:
                return str(rec.seq[start - 1:end])
    raise SystemExit(f"{contig} not in {fasta}")


def main() -> int:
    ap = argparse.ArgumentParser()
    for a in ("--loci_json", "--locus", "--ref", "--query", "--config", "--data_dir",
              "--gene_positions", "--cluster_tsv"):
        ap.add_argument(a, required=True)
    ap.add_argument("--id_sep", default="|")
    a = ap.parse_args()
    loc = next(x for x in json.loads(Path(a.loci_json).read_text())["loci"] if x["key"] == a.locus)
    fams = loc["families"]
    nl, nb = loc["n_left"], loc["n_locus"]
    flank = fams[:nl] + fams[nl + nb:]
    locus_fams = fams[nl:nl + nb]
    codes = next(r["codes"] for r in loc["rows"] if a.query in r["strains"])
    locs = load_gene_locations(a.cluster_tsv, a.gene_positions, set(fams), {a.ref, a.query}, a.id_sep)
    samples = {s.short: s for s in parse_config(a.config)}
    seqs = {}
    regs = {}
    for s in (a.ref, a.query):
        reg = region(locs, s, flank)
        if reg is None:
            raise SystemExit(f"{s}: no flank genes placed")
        regs[s] = reg
        dna = Path(a.data_dir) / "dna" / samples[s].dna
        seqs[s] = fetch(dna, *reg)
    r_contig, r_start, _ = regs[a.ref]
    with tempfile.TemporaryDirectory() as t:
        for s in (a.ref, a.query):
            Path(t, f"{s}.fa").write_text(f">{s}\n{seqs[s]}\n")
        out = subprocess.run(["blastn", "-query", f"{t}/{a.query}.fa", "-subject", f"{t}/{a.ref}.fa",
                              "-outfmt", "6 sstart send pident length", "-evalue", "1e-20"],
                             capture_output=True, text=True, check=True).stdout
    hsps = [tuple(float(v) for v in line.split("\t")) for line in out.splitlines() if line]
    print(f"{a.locus} ref {a.ref} {regs[a.ref]}  query {a.query} {regs[a.query]}  codes {codes}")
    agree = 0
    for i, f in enumerate(locus_fams):
        copies = [c for c in locs.get((a.ref, f), []) if c[1] == r_contig]
        if not copies:
            print(f"  {f}: reference has no copy on {r_contig}")
            continue
        _pid, _c, gs, ge = min(copies, key=lambda c: c[2])
        gs, ge = gs - r_start + 1, ge - r_start + 1
        cov = 0
        for s0, s1, _pid2, _len in hsps:
            lo, hi = min(s0, s1), max(s0, s1)
            cov = max(cov, max(0, min(hi, ge) - max(lo, gs) + 1))
        covered = cov >= 0.8 * (ge - gs + 1)
        call = codes[nl + i]
        expect_cov = call in "13"
        ok = covered == expect_cov
        agree += ok
        print(f"  {f}: call {call}, ref gene {gs:.0f}-{ge:.0f}, covered {covered} -> "
              f"{'agrees' if ok else 'DISAGREES'}")
    print(f"{agree}/{len(locus_fams)} locus columns agree")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```


- [ ] **Step 3: Read this first**

Pick the strains from `island_loci.json`: for each of loci `L001`-`L005`, one full-locus strain as `--ref`, one empty-site strain and one partial strain as `--query` (the first strain of the largest `empty` and `partial` rows). Run each pair:


- [ ] **Step 4: Run**

```bash
cd "${SCRATCH:?}/lv_real"
export PATH="/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/locus-view-a/.pixi/envs/default/bin:$PATH"
ST=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome
P=$ST/results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome
python - <<'EOF' > pairs.txt
import json
d = json.load(open("island_loci.json"))
for loc in d["loci"][:5]:
    rows = sorted(loc["rows"], key=lambda r: -r["count"])
    best = {c: next((r["strains"][0] for r in rows if r["row_class"] == c), None)
            for c in ("full", "empty", "partial")}
    for cls in ("empty", "partial"):
        if best["full"] and best[cls]:
            print(loc["key"], best["full"], best[cls], cls)
EOF
while read key ref qry cls; do
  python spotcheck_locus_calls.py --loci_json island_loci.json --locus $key --ref $ref --query $qry \
    --config $ST/config_immitis_in_posadasii_out.csv --data_dir $ST/data_dir \
    --gene_positions $P/gene_positions.tsv.zst --cluster_tsv $P/cluster/tier1_cluster.tsv | tail -4
done < pairs.txt | tee spotcheck.txt
```


- [ ] **Step 5: Read this first**

Expected output per pair: one line per locus column (`call 0, ... covered True -> DISAGREES` or `-> agrees`) and `N/M locus columns agree`. Paste `spotcheck.txt` into the PR description. Do not change code in response; the result goes to the user (Checkpoint).


## Checkpoint after Task 14 (stop here)

Open the Part A PR (after the user approves the push), then stop and report to the user:

1. Task 13's measured time, memory, sizes and top-locus counts next to the planning numbers.
2. Task 14's spot-check table.
3. The planning result: 0 of 6 empty-site calls on L001-L003 were confirmed by the DNA; the locus DNA was present in every "empty" strain. On these loci the informative ranking finds gene-model differences, not insertions or deletions.

Ask the user whether to start Part B as written, or to change the empty-site rule first. A DNA-level presence test is not in the spec, so it needs a spec change before any code. Do not start Part B without an answer.

---

## Part B: GenBank slices, clinker and the synteny panel (spec section 8)

Setup: after Part A is merged, `git -C /bigdata/stajichlab/jstajich/projects/NovInvenio fetch -q origin && git -C /bigdata/stajichlab/jstajich/projects/NovInvenio worktree add /bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/locus-view-b -b locus-view-b origin/main`. If the user wants Part B before Part A merges, branch from `locus-view-a` instead. Ask the user to approve one issue for Part B, then `export ISSUE_B=<number>`.

Part B ends in working software: each drawn locus has `pangenome/clinker/<key>.html`, the page shows it in the "Synteny (clinker)" panel, and NII publishes `clinker/` with the page.

### Task 15: Clinker strains and regions (spec section 8)

**Files:**
- Modify: `lib/island_locus.py` (keep `_placement` in `compute_locus`; append section 8)
- Create: `tests/test_island_locus_clinker.py`

**Interfaces:**
- Consumes: `compute_locus` result keys `_cells`, `_pairs`, `_row_class`, and the new `_placement` (the exemplar's `Placement`).
- Produces: `DEFAULT_CLINKER_MAX_STRAINS = 12`; `strain_region(result, strain, spans, flank, k=10) -> {'contig', 'rank_lo', 'rank_hi', 'anchors': [(rank, family)]} | None`; `select_clinker_strains(result, n50, species_of, spans, flank, k=10, max_strains=12) -> list[{'strain', 'reason', 'row_class', 'species', 'contig', 'rank_lo', 'rank_hi', 'anchors'}]` with `reason` in `exemplar | best_full | best_partial | best_empty | fill`.


- [ ] **Step 1: Write the failing test**

Create `tests/test_island_locus_clinker.py`:

```python
"""Clinker strain selection and per-strain regions (spec section 8)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from island_locus import (  # noqa: E402
    Columns, Locus, Placement, compute_locus, select_clinker_strains, strain_region,
)

FAMS = ["L1", "L2", "A", "B", "R1", "R2"]


class FakeMatrix:
    def __init__(self, calls):
        self.calls = calls

    def call(self, fam, strain):
        return self.calls.get((fam, strain), "absent")


def build(strains):
    """strains: {name: (contig, [rank or None per column], species)}."""
    pos, calls, spans, species = {}, {}, {}, {}
    for s, (contig, ranks, sp) in strains.items():
        species[s] = sp
        for fam, r in zip(FAMS, ranks):
            if r is not None:
                pos[(s, fam)] = [(contig, r)]
                calls[(fam, s)] = "present"
        spans[(s, contig)] = (0, 100)
    cols = Columns(left=("L1", "L2"), locus=("A", "B"), right=("R1", "R2"),
                   left_avail=20, right_avail=20)
    loc = Locus(root={"members": ["A", "B"], "locus_id": "X"}, variants=[{"n_strains": "3"}])
    ex = sorted(strains)[0]
    res = compute_locus(loc, Placement(ex, strains[ex][0], 22, 23, 20, 77, 101), "full", cols,
                        sorted(strains), pos, spans, FakeMatrix(calls), species)
    return res, spans, species


FULL = [20, 21, 22, 23, 24, 25]
EMPTY = [20, 21, None, None, 22, 23]
PART = [20, 21, 22, None, 24, 25]


def test_region_spans_in_place_copies_plus_flank_genes():
    res, spans, _ = build({"A1": ("c1", FULL, "sp")})
    assert strain_region(res, "A1", spans, flank=5) == {
        "contig": "c1", "rank_lo": 15, "rank_hi": 30,
        "anchors": [(20, "L1"), (21, "L2"), (22, "A"), (23, "B"), (24, "R1"), (25, "R2")]}


def test_empty_site_region_is_the_two_flank_blocks():
    res, spans, _ = build({"A1": ("c1", FULL, "sp"), "E1": ("c1", EMPTY, "sp")})
    reg = strain_region(res, "E1", spans, flank=5)
    assert (reg["rank_lo"], reg["rank_hi"]) == (15, 28)
    assert [f for _, f in reg["anchors"]] == ["L1", "L2", "R1", "R2"]


def test_region_is_clipped_at_the_contig_end():
    res, spans, _ = build({"A1": ("c1", [1, 2, 3, 4, 5, 6], "sp")})
    assert strain_region(res, "A1", spans, flank=5)["rank_lo"] == 0


def test_unanchored_strain_has_no_region():
    res, spans, _ = build({"A1": ("c1", FULL, "sp"), "U1": ("c2", [None, None, 5, None, None, None], "sp")})
    assert strain_region(res, "U1", spans, flank=5) is None


def test_selection_order_exemplar_then_class_species_then_fill():
    strains = {"A1": ("c1", FULL, "sp1"), "B1": ("c1", FULL, "sp2"), "C1": ("c1", PART, "sp1"),
               "D1": ("c1", EMPTY, "sp2"), "E1": ("c1", FULL, "sp1"), "F1": ("c1", PART, "sp1"),
               "U1": ("c2", [None, None, 5, None, None, None], "sp1")}
    res, spans, species = build(strains)
    picks = select_clinker_strains(res, {"E1": 900, "A1": 100, "C1": 50, "F1": 60}, species,
                                   spans, flank=5)
    assert [(p["strain"], p["reason"]) for p in picks] == [
        ("A1", "exemplar"), ("E1", "best_full"), ("B1", "best_full"), ("F1", "best_partial"),
        ("D1", "best_empty"), ("C1", "fill")]
    assert "U1" not in [p["strain"] for p in picks]


def test_selection_is_capped():
    strains = {f"S{i:02d}": ("c1", FULL, "sp") for i in range(20)}
    res, spans, species = build(strains)
    assert len(select_clinker_strains(res, {}, species, spans, flank=5, max_strains=12)) == 12


def test_each_pick_carries_its_region_and_class():
    res, spans, species = build({"A1": ("c1", FULL, "sp")})
    (pick,) = select_clinker_strains(res, {}, species, spans, flank=5)
    assert pick["row_class"] == "full" and pick["contig"] == "c1" and pick["rank_hi"] == 30
```


- [ ] **Step 2: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus_clinker.py
```

Expected: FAIL: `ImportError: cannot import name 'select_clinker_strains'`


- [ ] **Step 3: Make the change**

In `lib/island_locus.py`, replace this text (it occurs once):

```python
        "_pairs": pairs,
```

with:

```python
        "_pairs": pairs,
        "_placement": placement,
```

Append to the end of `lib/island_locus.py`:

```python


# ---- 8. clinker panel: strains and regions (spec section 8) --------------------
DEFAULT_CLINKER_MAX_STRAINS = 12
CLINKER_CLASSES = ("full", "partial", "empty")


def strain_region(result: dict, strain: str, spans: ContigSpans, flank: int,
                  k: int = DEFAULT_K) -> dict | None:
    """The rank range to draw for one strain (spec section 8, Ruling R11).

    On the contig of the strain's flank pair: from the lowest to the highest
    in-place column copy within n_locus + 2k ranks of the pair, then `flank`
    genes further on each side, clipped at the contig ends. The exemplar
    without a flank pair uses its Placement block instead. None for any
    other strain without a flank pair (it cannot be anchored)."""
    cells = result["_cells"][strain]
    pair = result["_pairs"].get(strain)
    fams = result["families"]
    if pair is not None:
        contig, p_lo, p_hi = pair
        w = result["n_locus"] + 2 * k
        anchors = sorted((r, fams[ci]) for ci, c, r in cells.in_place_copies
                         if c == contig and p_lo - w <= r <= p_hi + w)
        lo, hi = anchors[0][0], anchors[-1][0]
    elif strain == result["exemplar"]:
        place = result["_placement"]
        contig, lo, hi = place.contig, place.lo, place.hi
        anchors = sorted((r, fams[ci]) for ci, c, r in cells.in_place_copies
                         if c == contig and lo <= r <= hi)
    else:
        return None
    cmin, cmax = spans.get((strain, contig), (lo, hi))
    return {"contig": contig, "rank_lo": max(cmin, lo - flank), "rank_hi": min(cmax, hi + flank),
            "anchors": anchors}


def select_clinker_strains(result: dict, n50: dict[str, int], species_of: dict[str, str],
                           spans: ContigSpans, flank: int, k: int = DEFAULT_K,
                           max_strains: int = DEFAULT_CLINKER_MAX_STRAINS) -> list[dict]:
    """Up to `max_strains` strains for the clinker figure, in spec order:
    the exemplar; per row class (full, partial, empty) and species the
    flank-intact strain with the best quality_key(); then more full or
    partial strains by quality_key() until the cap. Uninformative strains
    are never chosen (the exemplar is kept even when its flanks are not
    intact). Each pick carries its reason and its strain_region()."""
    pairs = result["_pairs"]
    row_cls = result["_row_class"]

    def qkey(s):
        contig = pairs[s][0] if s in pairs else result["exemplar_contig"]
        cmin, cmax = spans.get((s, contig), (0, -1))
        return quality_key(s, n50, cmax - cmin + 1)

    picks: list[dict] = []
    chosen: set[str] = set()

    def add(strain, reason):
        if len(picks) >= max_strains or strain in chosen:
            return
        region = strain_region(result, strain, spans, flank, k)
        if region is None:
            return
        chosen.add(strain)
        picks.append({"strain": strain, "reason": reason, "row_class": row_cls[strain],
                      "species": species_of.get(strain, ""), **region})

    add(result["exemplar"], "exemplar")
    species = sorted({species_of.get(s, "") for s in pairs})
    for cls in CLINKER_CLASSES:
        for sp in species:
            pool = [s for s in pairs if row_cls[s] == cls and species_of.get(s, "") == sp
                    and s not in chosen]
            if pool:
                add(min(pool, key=qkey), "best_" + cls)
    fill = sorted((s for s in pairs if row_cls[s] in ("full", "partial") and s not in chosen),
                  key=qkey)
    for s in fill:
        if len(picks) >= max_strains:
            break
        add(s, "fill")
    return picks
```


- [ ] **Step 4: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus_clinker.py tests/test_island_locus.py
```

Expected: 48 passed


- [ ] **Step 5: Commit**

```bash
git add lib/island_locus.py tests/test_island_locus_clinker.py
git commit -m "clinker panel: choose strains and regions per locus (spec section 8) (#${ISSUE_B})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 16: `pangenome_island_loci.py --regions_out` and `--clinker_max_strains`

**Files:**
- Modify: `bin/pangenome_island_loci.py`
- Modify: `tests/test_pangenome_island_loci.py` (append)

**Interfaces:**
- Consumes: `select_clinker_strains` (Task 15).
- Produces: `island_regions.tsv` with columns `locus_key strain reason row_class species contig rank_lo rank_hi anchors` (`anchors` = `rank:family;...`); each page entry gains `clinker_strains` (the picks without `anchors`); `write_regions(path, regions)`, `REGION_COLUMNS`.


- [ ] **Step 1: Write the failing test**

Append to the end of `tests/test_pangenome_island_loci.py`:

```python


def test_regions_out_lists_clinker_strains_with_anchors(tmp_path):
    import csv
    run(tmp_path, "--regions_out", str(tmp_path / "regions.tsv"))
    rows = list(csv.DictReader(open(tmp_path / "regions.tsv"), delimiter="\t"))
    assert [(r["strain"], r["reason"]) for r in rows] == [
        ("S2", "exemplar"), ("S1", "best_full"), ("S3", "best_empty")]
    s3 = rows[2]
    assert (s3["contig"], s3["rank_lo"], s3["rank_hi"]) == ("c1", "0", "9")
    assert s3["anchors"].split(";")[0] == "0:F1"


def test_clinker_max_strains_zero_selects_none(tmp_path):
    data = run(tmp_path, "--clinker_max_strains", "0", "--regions_out", str(tmp_path / "r.tsv"))
    assert data["loci"][0]["clinker_strains"] == []
    assert (tmp_path / "r.tsv").read_text().count("\n") == 1
```


- [ ] **Step 2: Run**

```bash
pixi run python -m pytest -q tests/test_pangenome_island_loci.py
```

Expected: FAIL: `error: unrecognized arguments: --regions_out`


- [ ] **Step 3: Make the change**

In `bin/pangenome_island_loci.py`, replace this text (it occurs once):

```python
    candidate_loci, carrier_placements, choose_exemplar, compute_locus, group_loci,
    locus_columns, locus_payload, rank_key,
)
```

with:

```python
    candidate_loci, carrier_placements, choose_exemplar, compute_locus, group_loci,
    locus_columns, locus_payload, rank_key, select_clinker_strains,
)
```

In `bin/pangenome_island_loci.py`, replace this text (it occurs once):

```python
    out_loci = []
    for i, r in enumerate(drawn):
        dom_sets
```

with:

```python
    out_loci = []
    regions = []
    for i, r in enumerate(drawn):
        key = "L%03d" % (i + 1)
        picks = []
        if args.clinker_max_strains > 0:
            flank = args.flank_min if r["tier"] == "short_flanks" else args.flank
            picks = select_clinker_strains(r, n50, species_of, scan1.spans, flank, args.k,
                                           args.clinker_max_strains)
        for p in picks:
            regions.append({"locus_key": key, **p})
        dom_sets
```

In `bin/pangenome_island_loci.py`, replace this text (it occurs once):

```python
        out_loci.append(locus_payload(
            r, "L%03d" % (i + 1), bins,
            [dominant_class(d) for d in dom_sets], [",".join(d) for d in dom_sets],
            dominant_class(sorted({d for ds in dom_sets for d in ds})), locs, span))
```

with:

```python
        entry = locus_payload(
            r, key, bins,
            [dominant_class(d) for d in dom_sets], [",".join(d) for d in dom_sets],
            dominant_class(sorted({d for ds in dom_sets for d in ds})), locs, span)
        entry["clinker_strains"] = [{k: v for k, v in p.items() if k != "anchors"}
                                    for p in picks]
        out_loci.append(entry)
```

In `bin/pangenome_island_loci.py`, replace this text (it occurs once):

```python
        "loci": out_loci,
        "_results": drawn,
    }
```

with:

```python
        "loci": out_loci,
        "_results": drawn,
        "_regions": regions,
    }


REGION_COLUMNS = ["locus_key", "strain", "reason", "row_class", "species", "contig",
                  "rank_lo", "rank_hi", "anchors"]


def write_regions(path: str, regions: list[dict]) -> None:
    """island_regions.tsv: one row per (locus, clinker strain). `anchors` is
    'rank:family;...' for the strain's in-place column copies in the region,
    which bin/pangenome_island_gbk_slice.py checks against its own rank
    rebuild."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(REGION_COLUMNS)
        for r in regions:
            w.writerow([r["locus_key"], r["strain"], r["reason"], r["row_class"], r["species"],
                        r["contig"], r["rank_lo"], r["rank_hi"],
                        ";".join(f"{rank}:{fam}" for rank, fam in r["anchors"])])
```

In `bin/pangenome_island_loci.py`, replace this text (it occurs once):

```python
    ap.add_argument("--min_strains", type=int, default=2)
    ap.add_argument("--output", required=True)
```

with:

```python
    ap.add_argument("--min_strains", type=int, default=2)
    ap.add_argument("--clinker_max_strains", type=int, default=12,
                    help="strains per locus for the clinker panel; 0 selects none")
    ap.add_argument("--regions_out", default=None,
                    help="write island_regions.tsv (clinker strains and rank ranges)")
    ap.add_argument("--output", required=True)
```

In `bin/pangenome_island_loci.py`, replace this text (it occurs once):

```python
    payload = build(args)
    payload.pop("_results")
```

with:

```python
    payload = build(args)
    payload.pop("_results")
    regions = payload.pop("_regions")
    if args.regions_out:
        write_regions(args.regions_out, regions)
```


- [ ] **Step 4: Run**

```bash
pixi run python -m pytest -q tests/test_pangenome_island_loci.py
```

Expected: 13 passed


- [ ] **Step 5: Commit**

```bash
git add bin/pangenome_island_loci.py tests/test_pangenome_island_loci.py
git commit -m "clinker panel: write island_regions.tsv (#${ISSUE_B})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 17: GenBank slicing library `lib/genbank_slice.py`

**Files:**
- Create: `lib/genbank_slice.py`
- Create: `tests/test_genbank_slice.py`

**Interfaces:**
- Consumes: gene_positions rows `(strain, protein_id, contig, start, end)` and rescue_positions rows `(strain, family, contig, start)` in file order; GFF3 lines.
- Produces: `rank_entries(gene_rows, rescue_rows, wanted: set[(strain, contig)]) -> {(strain, contig): [(rank, protein_id | None, family | None, start, end)]}` (same ranks as `bin/pangenome_build_family_positions.py`); `cds_ids(attrs) -> list[str]`; `cds_exons(lines, contigs, ids) -> {pid: (strand, [(start, end)])}`; `build_record(strain, contig, contig_seq, entries, exons, proteins, family_of, label_of) -> (SeqRecord, summary)` with summary keys `bp_start bp_end n_genes n_rescue n_missing labels`; `safe_name(name) -> str`.


- [ ] **Step 1: Write the failing test**

Create `tests/test_genbank_slice.py`:

```python
"""lib/genbank_slice.py: rank rebuild, GFF3 CDS parsing, GenBank records."""
import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from Bio import SeqIO  # noqa: E402
from genbank_slice import build_record, cds_exons, cds_ids, rank_entries, safe_name  # noqa: E402
from pangenome_build_family_positions import build_family_positions  # noqa: E402
from pangenome_build_gene_positions import _resolve_ids  # noqa: E402

GENES = [
    ("S1", "p3", "c2", 50, 90), ("S1", "p1", "c1", 10, 40), ("S1", "p2", "c1", 60, 99),
    ("S2", "q1", "c1", 5, 30), ("S1", "p4", "c1", 200, 260), ("S1", "p5", "c10", 1, 9),
]
RESCUES = [("S1", "famR", "c1", 60), ("S1", "famZ", "c2", 10)]
MEMBER = {"S1|p1": "famA", "S1|p2": "famB", "S1|p3": "famC", "S1|p4": "famD", "S1|p5": "famE",
          "S2|q1": "famA"}


def test_rank_entries_match_build_family_positions():
    fp = build_family_positions(GENES, MEMBER, RESCUES)["S1"]
    rank_of = {(c, r): fam for fam, copies in fp.items() for c, r in copies}
    idx = rank_entries(iter(GENES), iter(RESCUES), {("S1", "c1"), ("S1", "c2"), ("S1", "c10")})
    for (strain, contig), entries in idx.items():
        for rank, pid, fam, _s, _e in entries:
            expected = rank_of[(contig, rank)]
            assert (MEMBER[f"S1|{pid}"] if pid else fam) == expected, (contig, rank)


def test_rank_entries_order_ties_gene_before_rescue():
    idx = rank_entries(iter(GENES), iter(RESCUES), {("S1", "c1")})
    assert [(e[1], e[2]) for e in idx[("S1", "c1")]] == [
        ("p1", None), ("p2", None), (None, "famR"), ("p4", None)]


def test_rank_offset_counts_contigs_that_sort_first():
    idx = rank_entries(iter(GENES), iter(RESCUES), {("S1", "c2")})
    # c1 (4 entries) and c10 (1 entry) sort before c2.
    assert [e[0] for e in idx[("S1", "c2")]] == [5, 6]


def test_cds_ids_matches_gene_positions_rule():
    for attrs in ["ID=x;Parent=m1", "protein_id=P1;Parent=m1", "Parent=a,b", "ID=x"]:
        assert cds_ids(attrs) == _resolve_ids(attrs)


GFF = """##gff-version 3
c1\tsrc\tgene\t10\t40\t.\t+\t.\tID=g1
c1\tsrc\tmRNA\t10\t40\t.\t+\t.\tID=p1;Parent=g1
c1\tsrc\tCDS\t10\t20\t.\t+\t0\tID=cds1;Parent=p1
c1\tsrc\tCDS\t30\t40\t.\t+\t0\tID=cds1;Parent=p1
c1\tsrc\tCDS\t60\t99\t.\t-\t0\tID=cds2;Parent=p2
c2\tsrc\tCDS\t50\t90\t.\t+\t0\tID=cds3;Parent=p3
"""


def test_cds_exons_keeps_wanted_ids_on_wanted_contigs():
    ex = cds_exons(io.StringIO(GFF), {"c1"}, {"p1", "p2", "p3"})
    assert ex == {"p1": (1, [(10, 20), (30, 40)]), "p2": (-1, [(60, 99)])}


def test_build_record_writes_gene_and_cds_features(tmp_path):
    seq = "A" * 300
    entries = rank_entries(iter(GENES), iter(RESCUES), {("S1", "c1")})[("S1", "c1")][:3]
    exons = cds_exons(io.StringIO(GFF), {"c1"}, {"p1", "p2"})
    rec, summ = build_record("S1", "c1", seq, entries, exons, {"p1": "MK*", "p2": "MR"},
                             {"p1": "famA", "p2": "famB"}, lambda pid: pid)
    assert (summ["bp_start"], summ["bp_end"]) == (10, 99)
    assert (summ["n_genes"], summ["n_rescue"], summ["n_missing"]) == (2, 1, 0)
    assert summ["labels"] == {"p1": "famA", "p2": "famB"}
    out = tmp_path / "S1.gbk"
    SeqIO.write(rec, str(out), "genbank")
    back = SeqIO.read(str(out), "genbank")
    cds = [f for f in back.features if f.type == "CDS"]
    assert [f.qualifiers["locus_tag"][0] for f in cds] == ["p1", "p2"]
    assert cds[0].qualifiers["translation"][0] == "MK"
    assert cds[0].qualifiers["note"][0] == "family=famA"
    assert len(cds[0].location.parts) == 2 and cds[1].location.strand == -1
    assert len([f for f in back.features if f.type == "gene"]) == 2


def test_gene_without_a_cds_row_is_counted_missing():
    entries = [(0, "pX", None, 10, 40)]
    _rec, summ = build_record("S1", "c1", "A" * 50, entries, {}, {}, {}, lambda p: p)
    assert (summ["n_genes"], summ["n_missing"]) == (0, 1)


def test_safe_name():
    assert safe_name("B0858-Guatemala") == "B0858-Guatemala"
    assert safe_name("a/b c") == "a_b_c"
```


- [ ] **Step 2: Run**

```bash
pixi run python -m pytest -q tests/test_genbank_slice.py
```

Expected: FAIL: `ModuleNotFoundError: No module named 'genbank_slice'`


- [ ] **Step 3: Make the change**

Create `lib/genbank_slice.py`:

```python
"""GenBank slices of per-strain locus regions, for clinker (spec section 8).

Pure functions, no I/O at import time; bin/pangenome_island_gbk_slice.py
does the reading and writing.

A region arrives as a family_positions rank range. rank_entries() rebuilds
which gene (or TBLASTN rescue hit) sits at each rank with exactly the rule
bin/pangenome_build_family_positions.py uses: per strain, every
gene_positions row in file order, then every rescue_positions row in file
order, stable-sorted by (contig, start) and enumerated. A rank is therefore
(number of the strain's entries on contigs that sort before this contig) +
(index within this contig).
"""
from __future__ import annotations

import collections
import re

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

_PROTEIN_ID_RE = re.compile(r"(?:^|;)protein_id=([^;\n]+)")
_PARENT_RE = re.compile(r"(?:^|;)Parent=([^;\n]+)")
_UNSAFE_RE = re.compile(r"[^A-Za-z0-9._-]")


def safe_name(name: str) -> str:
    """A file-name-safe form of a strain name."""
    return _UNSAFE_RE.sub("_", name) or "_"


def rank_entries(gene_rows, rescue_rows, wanted: set[tuple[str, str]]) -> dict:
    """{(strain, contig): [(rank, protein_id | None, family | None, start, end), ...]}
    for every (strain, contig) in `wanted`.

    `gene_rows` yields (strain, protein_id, contig, start, end) in
    gene_positions file order; `rescue_rows` yields (strain, family, contig,
    start) in rescue_positions file order. Entries of strains not in
    `wanted` are skipped; entries on other contigs of a wanted strain are
    only counted (they shift the rank offset)."""
    strains = {s for s, _ in wanted}
    counts: dict[tuple[str, str], int] = collections.Counter()
    kept: dict[tuple[str, str], list] = collections.defaultdict(list)
    for strain, pid, contig, start, end in gene_rows:
        if strain not in strains:
            continue
        counts[(strain, contig)] += 1
        if (strain, contig) in wanted:
            kept[(strain, contig)].append((start, pid, None, end))
    for strain, fam, contig, start in rescue_rows:
        if strain not in strains:
            continue
        counts[(strain, contig)] += 1
        if (strain, contig) in wanted:
            kept[(strain, contig)].append((start, None, fam, start))
    by_strain: dict[str, list[str]] = collections.defaultdict(list)
    for strain, contig in counts:
        by_strain[strain].append(contig)
    out = {}
    for strain, contig in wanted:
        offset = sum(counts[(strain, c)] for c in by_strain[strain] if c < contig)
        entries = sorted(kept.get((strain, contig), []), key=lambda e: e[0])
        out[(strain, contig)] = [(offset + i, pid, fam, start, end)
                                 for i, (start, pid, fam, end) in enumerate(entries)]
    return out


def cds_ids(attrs: str) -> list[str]:
    """IDs a GFF3 CDS row belongs to: protein_id= if present, else the
    comma-split Parent= (same rule as
    bin/pangenome_build_gene_positions.py::_resolve_ids)."""
    m = _PROTEIN_ID_RE.search(attrs)
    if m:
        return [m.group(1)]
    m = _PARENT_RE.search(attrs)
    if m:
        return [p for p in m.group(1).split(",") if p]
    return []


def cds_exons(gff3_lines, contigs: set[str], ids: set[str]) -> dict[str, tuple[int, list]]:
    """{protein_id: (strand, [(start, end), ...])} for CDS rows on `contigs`
    whose resolved ID is in `ids`. Coordinates are 1-based inclusive."""
    out: dict[str, list] = {}
    for line in gff3_lines:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS" or f[0] not in contigs:
            continue
        for pid in cds_ids(f[8]):
            if pid not in ids:
                continue
            entry = out.setdefault(pid, [-1 if f[6] == "-" else 1, []])
            entry[1].append((int(f[3]), int(f[4])))
    return {pid: (strand, sorted(ex)) for pid, (strand, ex) in out.items()}


def build_record(strain: str, contig: str, contig_seq: str, entries: list,
                 exons: dict, proteins: dict[str, str], family_of: dict[str, str],
                 label_of) -> tuple[SeqRecord, dict]:
    """One GenBank record for a region.

    `entries` are rank_entries() rows inside the region. Genes become a
    `gene` + `CDS` feature pair (the gene feature stops clinker's
    "Could not find parent gene" warning, spike problem 5); rescue hits have
    no gene model and are not drawn. `label_of(protein_id)` gives the
    /locus_tag. Returns (record, summary) with summary keys bp_start,
    bp_end, n_genes, n_rescue, n_missing, labels ({label: family})."""
    starts = [e[3] for e in entries]
    ends = [e[4] for e in entries]
    bp_start, bp_end = min(starts), max(ends)
    sub = contig_seq[bp_start - 1:bp_end]
    rec = SeqRecord(Seq(sub), id=safe_name(strain)[:16], name=safe_name(strain)[:16],
                    description=f"{strain} {contig}:{bp_start}-{bp_end}")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "linear"
    labels: dict[str, str] = {}
    n_genes = n_rescue = n_missing = 0
    for _rank, pid, fam, _s, _e in entries:
        if pid is None:
            n_rescue += 1
            continue
        if pid not in exons:
            n_missing += 1
            continue
        strand, exs = exons[pid]
        parts = []
        for s, e in exs:
            ns, ne = max(0, s - bp_start), min(len(sub), e - bp_start + 1)
            if ne > ns:
                parts.append(FeatureLocation(ns, ne, strand=strand))
        if not parts:
            n_missing += 1
            continue
        if strand == -1:
            parts = parts[::-1]
        loc = parts[0] if len(parts) == 1 else CompoundLocation(parts)
        label = label_of(pid)
        family = family_of.get(pid, "")
        quals = {"locus_tag": [label], "note": [f"family={family}"]}
        rec.features.append(SeqFeature(FeatureLocation(loc.start, loc.end, strand=strand),
                                       type="gene", qualifiers={"locus_tag": [label]}))
        cds_q = dict(quals)
        if proteins.get(pid):
            cds_q["translation"] = [proteins[pid].rstrip("*")]
        rec.features.append(SeqFeature(loc, type="CDS", qualifiers=cds_q))
        labels[label] = family
        n_genes += 1
    return rec, {"bp_start": bp_start, "bp_end": bp_end, "n_genes": n_genes,
                 "n_rescue": n_rescue, "n_missing": n_missing, "labels": labels}
```


- [ ] **Step 4: Run**

```bash
pixi run python -m pytest -q tests/test_genbank_slice.py
```

Expected: 8 passed


- [ ] **Step 5: Commit**

```bash
git add lib/genbank_slice.py tests/test_genbank_slice.py
git commit -m "clinker panel: GenBank slices with the family_positions rank rule (#${ISSUE_B})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 18: CLI `bin/pangenome_island_gbk_slice.py`

**Files:**
- Create: `bin/pangenome_island_gbk_slice.py`
- Create: `tests/test_pangenome_island_gbk_slice.py`

**Interfaces:**
- Consumes: `island_regions.tsv` (Task 16), samplesheet, data_dir, gff3_dir, `gene_positions.tsv[.zst]`, `rescue_positions.tsv`, tier-1 cluster TSV; `lib/genbank_slice.py`.
- Produces: `<out_dir>/<key>/<safe strain>.gbk`, `<out_dir>/<key>/groups.csv` (`locus_tag,family`, no header), `<out_dir>/island_slices.tsv` (`locus_key strain contig bp_start bp_end n_genes n_rescue n_missing`); exit 1 when an anchor family does not match the rank rebuild.


- [ ] **Step 1: Write the failing test**

Create `tests/test_pangenome_island_gbk_slice.py`:

```python
"""CLI tests for bin/pangenome_island_gbk_slice.py on a one-strain fixture."""
import csv
import sys
from pathlib import Path

BIN = Path(__file__).parent.parent / "bin"
sys.path.insert(0, str(BIN))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from Bio import SeqIO  # noqa: E402
from pangenome_island_gbk_slice import main  # noqa: E402

GFF = """##gff-version 3
c1\tsrc\tCDS\t10\t40\t.\t+\t0\tParent=p1
c1\tsrc\tCDS\t60\t99\t.\t-\t0\tParent=p2
c1\tsrc\tCDS\t120\t150\t.\t+\t0\tParent=p3
"""
REGION_HEADER = "locus_key\tstrain\treason\trow_class\tspecies\tcontig\trank_lo\trank_hi\tanchors\n"


def fixture(d: Path, anchors="0:famA;1:famB") -> list[str]:
    for sub in ("dna", "pep", "gff3"):
        (d / "data" / sub).mkdir(parents=True, exist_ok=True)
    (d / "data" / "dna" / "S1.dna.fa").write_text(">c1\n" + "ACGT" * 50 + "\n>c2\nAAAA\n")
    (d / "data" / "pep" / "S1.pep.fa").write_text(">p1\nMKV*\n>p2\nMRR\n>p3\nMQQ\n")
    (d / "data" / "gff3" / "S1.gff3").write_text(GFF)
    (d / "config.csv").write_text("GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup\n"
                                  "IN,Sp,S1,S1.pep.fa,S1.dna.fa,S1.gff3,S1,t\n")
    (d / "gene_positions.tsv").write_text("Short\tprotein_id\tcontig\tstart\tend\n"
                                          "S1\tp1\tc1\t10\t40\nS1\tp2\tc1\t60\t99\n"
                                          "S1\tp3\tc1\t120\t150\n")
    (d / "rescue.tsv").write_text("Short\tfamily\tcontig\tstart\nS1\tfamR\tc1\t110\n")
    (d / "cluster.tsv").write_text("famA\tS1|p1\nfamB\tS1|p2\nfamC\tS1|p3\n")
    (d / "regions.tsv").write_text(REGION_HEADER +
                                   f"L001\tS1\texemplar\tfull\tSp\tc1\t0\t3\t{anchors}\n")
    return ["--regions", str(d / "regions.tsv"), "--config", str(d / "config.csv"),
            "--data_dir", str(d / "data"), "--gff3_dir", str(d / "data" / "gff3"),
            "--gene_positions", str(d / "gene_positions.tsv"),
            "--rescue_positions", str(d / "rescue.tsv"),
            "--cluster_tsv", str(d / "cluster.tsv"), "--out_dir", str(d / "gbk")]


def test_writes_one_genbank_per_region_with_families(tmp_path):
    assert main(fixture(tmp_path)) == 0
    rec = SeqIO.read(str(tmp_path / "gbk" / "L001" / "S1.gbk"), "genbank")
    cds = [f for f in rec.features if f.type == "CDS"]
    assert [f.qualifiers["note"][0] for f in cds] == ["family=famA", "family=famB", "family=famC"]
    assert cds[0].qualifiers["translation"][0] == "MKV"
    assert len(rec.seq) == 150 - 10 + 1


def test_groups_csv_maps_locus_tags_to_families(tmp_path):
    main(fixture(tmp_path))
    rows = list(csv.reader(open(tmp_path / "gbk" / "L001" / "groups.csv")))
    assert rows == [["p1", "famA"], ["p2", "famB"], ["p3", "famC"]]


def test_slices_table_reports_bp_and_counts(tmp_path):
    main(fixture(tmp_path))
    (row,) = list(csv.DictReader(open(tmp_path / "gbk" / "island_slices.tsv"), delimiter="\t"))
    assert (row["bp_start"], row["bp_end"], row["n_genes"], row["n_rescue"]) == ("10", "150", "3", "1")


def test_anchor_mismatch_fails(tmp_path):
    assert main(fixture(tmp_path, anchors="0:famB")) == 1


def test_strain_without_files_is_skipped(tmp_path):
    args = fixture(tmp_path)
    (tmp_path / "data" / "dna" / "S1.dna.fa").unlink()
    assert main(args) == 0
    lines = (tmp_path / "gbk" / "island_slices.tsv").read_text().splitlines()
    assert len(lines) == 1


def test_no_regions_writes_a_header_only_table(tmp_path):
    args = fixture(tmp_path)
    (tmp_path / "regions.tsv").write_text(REGION_HEADER)
    assert main(args) == 0
    assert (tmp_path / "gbk" / "island_slices.tsv").read_text().startswith("locus_key\t")


def test_protein_id_clash_across_strains_prefixes_labels(tmp_path):
    args = fixture(tmp_path)
    d = tmp_path
    for sub, ext in {"dna": "dna.fa", "pep": "pep.fa", "gff3": "gff3"}.items():
        text = (d / "data" / sub / f"S1.{ext}").read_text()
        (d / "data" / sub / f"S2.{ext}").write_text(text)
    with open(d / "config.csv", "a") as fh:
        fh.write("IN,Sp,S2,S2.pep.fa,S2.dna.fa,S2.gff3,S2,t\n")
    with open(d / "gene_positions.tsv", "a") as fh:
        fh.write("S2\tp1\tc1\t10\t40\nS2\tp2\tc1\t60\t99\nS2\tp3\tc1\t120\t150\n")
    with open(d / "cluster.tsv", "a") as fh:
        fh.write("famA\tS2|p1\nfamB\tS2|p2\nfamC\tS2|p3\n")
    with open(d / "regions.tsv", "a") as fh:
        fh.write("L001\tS2\tfill\tfull\tSp\tc1\t0\t2\t0:famA\n")
    assert main(args) == 0
    labels = [r[0] for r in csv.reader(open(d / "gbk" / "L001" / "groups.csv"))]
    assert "S1|p1" in labels and "S2|p1" in labels


def test_gzipped_genome_and_proteins_are_read(tmp_path):
    # Review Focus 4: study FASTAs may be .gz; lib/compressed_io handles them.
    import gzip
    args = fixture(tmp_path)
    for sub, name in (("dna", "S1.dna.fa"), ("pep", "S1.pep.fa")):
        src = tmp_path / "data" / sub / name
        with gzip.open(str(src) + ".gz", "wt") as fh:
            fh.write(src.read_text())
        src.unlink()
    cfg = tmp_path / "config.csv"
    cfg.write_text(cfg.read_text().replace("S1.pep.fa,S1.dna.fa", "S1.pep.fa.gz,S1.dna.fa.gz"))
    assert main(args) == 0
    rec = SeqIO.read(str(tmp_path / "gbk" / "L001" / "S1.gbk"), "genbank")
    assert len([f for f in rec.features if f.type == "CDS"]) == 3
```


- [ ] **Step 2: Run**

```bash
pixi run python -m pytest -q tests/test_pangenome_island_gbk_slice.py
```

Expected: FAIL: `ModuleNotFoundError: No module named 'pangenome_island_gbk_slice'`


- [ ] **Step 3: Make the change**

Create `bin/pangenome_island_gbk_slice.py`:

```python
#!/usr/bin/env python3
"""ISLAND_GBK_SLICE: one GenBank file per (locus, strain) region listed in
island_regions.tsv (bin/pangenome_island_loci.py --regions_out), for clinker.

Writes, under --out_dir:
  <locus_key>/<strain>.gbk   gene + CDS features, /locus_tag, /translation,
                             /note="family=<tier-1 family>"
  <locus_key>/groups.csv     locus_tag,family (clinker -gf; no header)
  island_slices.tsv          locus_key, strain, contig, bp_start, bp_end,
                             n_genes, n_rescue, n_missing
Exits non-zero if a region's anchor families do not match the families at
those ranks after the rank rebuild (lib/genbank_slice.rank_entries), since
that means the regions and the slices describe different genes.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from Bio import SeqIO  # noqa: E402
from compressed_io import open_maybe_compressed  # noqa: E402
from config_parser import parse_config  # noqa: E402
from genbank_slice import build_record, cds_exons, rank_entries, safe_name  # noqa: E402

DNA_SUBDIRS = ["dna", "genome", "scaffolds"]
PEP_SUBDIRS = ["pep", "proteins"]


def resolve_under(data_dir: Path, name: str, subdirs: list[str]) -> Path | None:
    """Flat layout first, then each subdir -- pangenome.nf's resolve_fa()."""
    for cand in [data_dir / name] + [data_dir / sub / name for sub in subdirs]:
        if cand.exists():
            return cand
    return None


def read_regions(path: str) -> list[dict]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["rank_lo"], r["rank_hi"] = int(r["rank_lo"]), int(r["rank_hi"])
        r["anchor_list"] = [(int(a.split(":", 1)[0]), a.split(":", 1)[1])
                            for a in r["anchors"].split(";") if a]
    return rows


def iter_gene_rows(path: str):
    with open_maybe_compressed(path) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            try:
                yield p[0], p[1], p[2], int(p[3]), int(p[4])
            except ValueError:
                continue


def iter_rescue_rows(path: str | None):
    if not path or not Path(path).is_file() or Path(path).stat().st_size == 0:
        return
    with open_maybe_compressed(path) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            try:
                yield p[0], p[1], p[2], int(p[3])
            except ValueError:
                continue


def read_families(cluster_tsv: str, members: set[str]) -> dict[str, str]:
    """{member_id: family} for the wanted `Short<sep>protein` member IDs."""
    out = {}
    with open_maybe_compressed(cluster_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1] in members:
                out[parts[1]] = parts[0]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--data_dir", required=True)
    ap.add_argument("--gff3_dir", required=True)
    ap.add_argument("--gene_positions", required=True)
    ap.add_argument("--rescue_positions", default=None)
    ap.add_argument("--cluster_tsv", required=True)
    ap.add_argument("--id_sep", default="|")
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args(argv)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    regions = read_regions(args.regions)
    slices_path = out_dir / "island_slices.tsv"
    header = ["locus_key", "strain", "contig", "bp_start", "bp_end", "n_genes", "n_rescue",
              "n_missing"]
    if not regions:
        slices_path.write_text("\t".join(header) + "\n")
        print("pangenome_island_gbk_slice: no regions", file=sys.stderr)
        return 0

    samples = {s.short: s for s in parse_config(args.config)}
    wanted = {(r["strain"], r["contig"]) for r in regions}
    index = rank_entries(iter_gene_rows(args.gene_positions),
                         iter_rescue_rows(args.rescue_positions), wanted)
    region_entries = []
    for r in regions:
        entries = [e for e in index.get((r["strain"], r["contig"]), [])
                   if r["rank_lo"] <= e[0] <= r["rank_hi"]]
        region_entries.append(entries)
    members = {f"{r['strain']}{args.id_sep}{e[1]}"
               for r, ents in zip(regions, region_entries) for e in ents if e[1]}
    fam_of_member = read_families(args.cluster_tsv, members)

    errors = []
    for r, ents in zip(regions, region_entries):
        at = {e[0]: (fam_of_member.get(f"{r['strain']}{args.id_sep}{e[1]}") if e[1] else e[2])
              for e in ents}
        for rank, fam in r["anchor_list"]:
            if at.get(rank) != fam:
                errors.append(f"{r['locus_key']} {r['strain']} rank {rank}: region says "
                              f"{fam}, rank rebuild gives {at.get(rank)}")
    if errors:
        print("ERROR: rank rebuild does not match island_regions.tsv:\n  " +
              "\n  ".join(errors[:20]), file=sys.stderr)
        return 1

    by_strain: dict[str, list[int]] = {}
    for i, r in enumerate(regions):
        by_strain.setdefault(r["strain"], []).append(i)
    rows_out = []
    group_rows: dict[str, dict[str, str]] = {}
    pid_owner: dict[str, dict[str, set[str]]] = {}
    for i, r in enumerate(regions):
        for e in region_entries[i]:
            if e[1]:
                pid_owner.setdefault(r["locus_key"], {}).setdefault(e[1], set()).add(r["strain"])
    for strain, idxs in sorted(by_strain.items()):
        s = samples.get(strain)
        if s is None:
            print(f"WARNING: {strain} not in --config; skipped", file=sys.stderr)
            continue
        dna = resolve_under(Path(args.data_dir), s.dna, DNA_SUBDIRS)
        pep = resolve_under(Path(args.data_dir), s.protein, PEP_SUBDIRS)
        gff = Path(args.gff3_dir) / s.gff3 if s.gff3 else None
        if dna is None or pep is None or gff is None or not gff.exists():
            print(f"WARNING: {strain}: DNA, protein or GFF3 file not found; skipped",
                  file=sys.stderr)
            continue
        contigs = {regions[i]["contig"] for i in idxs}
        pids = {e[1] for i in idxs for e in region_entries[i] if e[1]}
        with open_maybe_compressed(str(dna)) as fh:
            seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(fh, "fasta") if rec.id in contigs}
        with open_maybe_compressed(str(pep)) as fh:
            prots = {rec.id: str(rec.seq) for rec in SeqIO.parse(fh, "fasta") if rec.id in pids}
        with open_maybe_compressed(str(gff)) as fh:
            exons = cds_exons(fh, contigs, pids)
        for i in idxs:
            r = regions[i]
            ents = region_entries[i]
            if not ents or r["contig"] not in seqs:
                print(f"WARNING: {r['locus_key']} {strain}: no genes or contig "
                      f"{r['contig']} not in {dna}; skipped", file=sys.stderr)
                continue
            clash = any(len(owners) > 1 for owners in pid_owner.get(r["locus_key"], {}).values())
            family_of = {e[1]: fam_of_member.get(f"{strain}{args.id_sep}{e[1]}", "")
                         for e in ents if e[1]}
            rec, summ = build_record(
                strain, r["contig"], seqs[r["contig"]], ents, exons, prots, family_of,
                (lambda pid, st=strain: f"{st}{args.id_sep}{pid}") if clash else (lambda pid: pid))
            locus_dir = out_dir / r["locus_key"]
            locus_dir.mkdir(exist_ok=True)
            SeqIO.write(rec, str(locus_dir / f"{safe_name(strain)}.gbk"), "genbank")
            group_rows.setdefault(r["locus_key"], {}).update(summ["labels"])
            rows_out.append([r["locus_key"], strain, r["contig"], summ["bp_start"],
                             summ["bp_end"], summ["n_genes"], summ["n_rescue"], summ["n_missing"]])
    for key, labels in group_rows.items():
        with open(out_dir / key / "groups.csv", "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            for label, fam in sorted(labels.items()):
                w.writerow([label, fam])
    with open(slices_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(sorted(rows_out))
    print(f"pangenome_island_gbk_slice: {len(rows_out)} slices for {len(group_rows)} loci",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

Make it executable: `chmod +x bin/pangenome_island_gbk_slice.py`.


- [ ] **Step 4: Run**

```bash
pixi run python -m pytest -q tests/test_pangenome_island_gbk_slice.py
```

Expected: 8 passed


- [ ] **Step 5: Commit**

```bash
git add bin/pangenome_island_gbk_slice.py tests/test_pangenome_island_gbk_slice.py
git commit -m "clinker panel: ISLAND_GBK_SLICE script (#${ISSUE_B})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 19: clinker 0.0.32 in the environment, slimming, and the per-locus runner

**Files:**
- Modify: `pixi.toml`, `pixi.lock`
- Modify: `Dockerfile`, `novinvenio.def`
- Create: `lib/clinker_html.py`, `bin/pangenome_island_clinker.py`
- Create: `tests/test_clinker_html.py`, `tests/test_pangenome_island_clinker.py`

**Interfaces:**
- Consumes: `<key>/` dirs from Task 18.
- Produces: `slim_clinker_html(html: str) -> str` (raises `ValueError` on a non-clinker page); `bin/pangenome_island_clinker.py --locus_dirs D [D ...] --out_dir . --cpus N [--clinker EXE] [--keep_sequences]` writing `<out_dir>/<key>.html`, skipping a failed locus with a warning (Ruling R15); `clinker` on PATH in the pixi env.


- [ ] **Step 1: Make the change**

In `pixi.toml`, replace this text (it occurs once):

```toml
leidenalg = ">=0.10,<1"
```

with:

```toml
leidenalg = ">=0.10,<1"

# gamcil/clinker (gene-cluster synteny figures, ISLAND_CLINKER) is on PyPI only.
# Bioconda's "clinker" is an unrelated RNA-seq tool: never add it under
# [dependencies].
[pypi-dependencies]
clinker = "==0.0.32"
```

In `Dockerfile`, replace this text (it occurs once):

```bash
        mash \
        zstd

# OpenMPI needs a writable /tmp.
```

with:

```bash
        mash \
        zstd

# gamcil/clinker 0.0.32 for ISLAND_CLINKER is on PyPI only (bioconda's
# "clinker" is an unrelated RNA-seq tool).
RUN pip install --no-cache-dir clinker==0.0.32

# OpenMPI needs a writable /tmp.
```

In `novinvenio.def`, replace this text (it occurs once):

```bash
        mash \
        zstd

%environment
```

with:

```bash
        mash \
        zstd
    pip install --no-cache-dir clinker==0.0.32

%environment
```


- [ ] **Step 2: Run**

```bash
pixi install && pixi run clinker --version
```

Expected: `clinker v0.0.32` (checked during planning with a minimal pixi workspace: the PyPI package resolves under `[pypi-dependencies]`). `pixi.lock` changes; commit it.


- [ ] **Step 3: Write the failing test**

Create `tests/test_clinker_html.py`:

```python
"""lib/clinker_html.py: strip sequences from clinker's embedded data."""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from clinker_html import slim_clinker_html  # noqa: E402


def gene(uid):
    return {"uid": uid, "label": uid, "names": {"locus_tag": uid}, "start": 0, "end": 9,
            "strand": 1, "sequence": "ATG" * 50, "translation": "M" * 50}


DATA = {"clusters": [{"uid": "c1", "name": "S1", "loci": [{"uid": "l1", "genes": [gene("g1")]}]}],
        "links": [{"uid": "k1", "query": gene("g1"), "target": gene("g2"), "identity": 0.9}],
        "groups": [{"uid": "gr1", "label": "famA", "genes": ["g1", "g2"]}]}
PAGE = ("<html><head></head><body><div id=\"plot\"></div><script>const data=" +
        json.dumps(DATA) + ";function serialise(svg) { return 1; }\nplot(data)</script></body></html>")


def embedded(html):
    start = html.index("const data=") + len("const data=")
    return json.JSONDecoder().raw_decode(html, start)[0]


def test_sequences_are_removed_from_clusters_and_links():
    data = embedded(slim_clinker_html(PAGE))
    genes = [data["clusters"][0]["loci"][0]["genes"][0], data["links"][0]["query"],
             data["links"][0]["target"]]
    for g in genes:
        assert "sequence" not in g and "translation" not in g
        assert g["label"] and g["end"] == 9


def test_everything_else_is_kept():
    slim = slim_clinker_html(PAGE)
    assert slim.startswith("<html><head></head><body><div id=\"plot\"></div><script>const data=")
    assert slim.endswith(";function serialise(svg) { return 1; }\nplot(data)</script></body></html>")
    data = embedded(slim)
    assert data["groups"] == DATA["groups"]
    assert data["links"][0]["identity"] == 0.9


def test_a_closing_script_tag_in_a_label_is_escaped():
    d = json.loads(json.dumps(DATA))
    d["groups"][0]["label"] = "</script><b>x"
    page = PAGE.replace(json.dumps(DATA), json.dumps(d))
    blob = slim_clinker_html(page).split("const data=", 1)[1].split(";function", 1)[0]
    assert "</script>" not in blob


def test_non_clinker_page_raises():
    with pytest.raises(ValueError):
        slim_clinker_html("<html></html>")
```

Create `tests/test_pangenome_island_clinker.py`:

```python
"""bin/pangenome_island_clinker.py with a stand-in clinker executable."""
import stat
import sys
from pathlib import Path

BIN = Path(__file__).parent.parent / "bin"
sys.path.insert(0, str(BIN))
from pangenome_island_clinker import main  # noqa: E402

# Writes a minimal clinker-shaped page to the -p argument; fails for L002.
FAKE = """#!/usr/bin/env python3
import json, sys
args = sys.argv[1:]
out = args[args.index("-p") + 1]
if "L002" in args[0]:
    sys.exit(3)
gene = {"uid": "g", "label": "g", "sequence": "ATG" * 30, "translation": "M" * 30}
data = {"clusters": [{"loci": [{"genes": [gene]}]}], "links": [], "groups": []}
open(out, "w").write("<script>const data=" + json.dumps(data) + ";plot(data)</script>")
"""


def setup(tmp_path):
    fake = tmp_path / "fake_clinker"
    fake.write_text(FAKE)
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    dirs = []
    for key in ("L001", "L002", "L003"):
        d = tmp_path / "gbk" / key
        d.mkdir(parents=True)
        (d / "S1.gbk").write_text("LOCUS\n")
        (d / "groups.csv").write_text("g,famA\n")
        dirs.append(str(d))
    (tmp_path / "gbk" / "L003" / "S1.gbk").unlink()
    return fake, dirs


def test_each_locus_gets_a_slim_page_and_failures_are_skipped(tmp_path, capsys):
    fake, dirs = setup(tmp_path)
    out = tmp_path / "out"
    assert main(["--locus_dirs", *dirs, "--out_dir", str(out), "--clinker", str(fake)]) == 0
    assert sorted(p.name for p in out.iterdir()) == ["L001.html"]
    assert "ATGATG" not in (out / "L001.html").read_text()
    err = capsys.readouterr().err
    assert "L002: clinker exited 3" in err and "L003: no .gbk files" in err


def test_keep_sequences_writes_the_page_unchanged(tmp_path):
    fake, dirs = setup(tmp_path)
    out = tmp_path / "out"
    main(["--locus_dirs", dirs[0], "--out_dir", str(out), "--clinker", str(fake), "--keep_sequences"])
    assert "ATGATG" in (out / "L001.html").read_text()


def test_missing_clinker_executable_skips_every_locus(tmp_path, capsys):
    _fake, dirs = setup(tmp_path)
    out = tmp_path / "out"
    assert main(["--locus_dirs", dirs[0], "--out_dir", str(out),
                 "--clinker", str(tmp_path / "no_such_clinker")]) == 0
    assert list(out.iterdir()) == []
    assert "cannot run" in capsys.readouterr().err
```


- [ ] **Step 4: Run**

```bash
pixi run python -m pytest -q tests/test_clinker_html.py tests/test_pangenome_island_clinker.py
```

Expected: FAIL: `ModuleNotFoundError: No module named 'clinker_html'`


- [ ] **Step 5: Make the change**

Create `lib/clinker_html.py`:

```python
"""Slim clinker 0.0.32 HTML pages (spec section 8).

clinker embeds its data as `const data={...};` in one inline <script> and
copies both genes' DNA and protein sequence into every link. clinker's own
page scripts never read the `sequence` or `translation` fields (checked in
clinker 0.0.32's plot/clinker.js and clustermap.min.js), so removing them
only removes bytes: 4.13 MB -> 0.81 MB on the 2026-09-26 spike page and
5.35 MB -> 1.01 MB on a real 12-strain locus, with the same drawn figure
(tests/test_clinker_render.py). No I/O at import time.
"""
from __future__ import annotations

import json

MARKER = "const data="
STRIP = ("sequence", "translation")


def _strip_gene(gene: dict) -> None:
    for key in STRIP:
        gene.pop(key, None)


def slim_clinker_html(html: str) -> str:
    """The page with sequence/translation removed from every gene in
    data.clusters[].loci[].genes[] and data.links[].query/target. `</` in
    the re-serialised data is escaped so a label cannot close the script.
    Raises ValueError when the page has no `const data=` object."""
    start = html.find(MARKER)
    if start < 0:
        raise ValueError("no 'const data=' block: not a clinker HTML page")
    begin = start + len(MARKER)
    data, end = json.JSONDecoder().raw_decode(html, begin)
    for cluster in data.get("clusters", []):
        for locus in cluster.get("loci", []):
            for gene in locus.get("genes", []):
                _strip_gene(gene)
    for link in data.get("links", []):
        for side in ("query", "target"):
            if isinstance(link.get(side), dict):
                _strip_gene(link[side])
    blob = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    return html[:begin] + blob + html[end:]
```

Create `bin/pangenome_island_clinker.py`:

```python
#!/usr/bin/env python3
"""ISLAND_CLINKER: one clinker figure per locus directory (spec section 8).

For each --locus_dirs entry (<key>/ holding <strain>.gbk files and
groups.csv from bin/pangenome_island_gbk_slice.py) runs

    clinker <key>/*.gbk -gf <key>/groups.csv -j <cpus> -p <tmp>.html

and writes <out_dir>/<key>.html, slimmed by lib/clinker_html.py unless
--keep_sequences. A locus whose clinker run fails (or a missing clinker
executable) is skipped with a warning and no page, so one bad locus never
costs the others (plan Ruling R15); the synteny page then says "No synteny
figure for this locus". Exit status is 0 unless an argument is wrong.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from clinker_html import slim_clinker_html  # noqa: E402


def run_locus(locus_dir: Path, out_dir: Path, clinker: str, cpus: int, keep: bool) -> bool:
    key = locus_dir.name
    gbks = sorted(str(p) for p in locus_dir.glob("*.gbk"))
    if not gbks:
        print(f"WARNING: {key}: no .gbk files; skipped", file=sys.stderr)
        return False
    raw = out_dir / f"{key}.raw.html"
    cmd = [clinker, *gbks, "-gf", str(locus_dir / "groups.csv"), "-j", str(cpus), "-p", str(raw)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except OSError as e:
        print(f"WARNING: {key}: cannot run {clinker}: {e}; skipped", file=sys.stderr)
        return False
    if proc.returncode != 0 or not raw.is_file():
        print(f"WARNING: {key}: clinker exited {proc.returncode}; skipped\n{proc.stderr[-2000:]}",
              file=sys.stderr)
        raw.unlink(missing_ok=True)
        return False
    html = raw.read_text()
    (out_dir / f"{key}.html").write_text(html if keep else slim_clinker_html(html))
    raw.unlink()
    return True


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--locus_dirs", nargs="+", required=True)
    ap.add_argument("--out_dir", default=".")
    ap.add_argument("--cpus", type=int, default=1)
    ap.add_argument("--clinker", default="clinker", help="clinker executable")
    ap.add_argument("--keep_sequences", action="store_true",
                    help="write clinker's page unchanged (--pangenome_clinker_slim false)")
    args = ap.parse_args(argv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    done = sum(run_locus(Path(d), out_dir, args.clinker, args.cpus, args.keep_sequences)
               for d in args.locus_dirs)
    print(f"pangenome_island_clinker: {done} of {len(args.locus_dirs)} loci drawn", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

Make it executable: `chmod +x bin/pangenome_island_clinker.py`.


- [ ] **Step 6: Run**

```bash
pixi run python -m pytest -q tests/test_clinker_html.py tests/test_pangenome_island_clinker.py
```

Expected: 7 passed


- [ ] **Step 7: Read this first**

The container image is not rebuilt by this task. The published image `ghcr.io/stajichlab/novinvenio:0.5.0` has no clinker, so under `-profile docker|singularity` every locus is skipped with a warning until a rebuilt image is pushed (a manual step with registry credentials; tell the user). Local and SLURM runs use the pixi env and are not affected.


- [ ] **Step 8: Commit**

```bash
git add pixi.toml pixi.lock Dockerfile novinvenio.def lib/clinker_html.py bin/pangenome_island_clinker.py tests/test_clinker_html.py tests/test_pangenome_island_clinker.py
git commit -m "clinker panel: clinker 0.0.32 from PyPI, slimming and per-locus runner (#${ISSUE_B})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 20: Check that the slimmed clinker page renders the same figure

**Files:**
- Create: `tests/test_clinker_render.py`

**Interfaces:**
- Consumes: `clinker` (Task 19), `slim_clinker_html`, a headless Chromium.
- Produces: an opt-in test that compares drawn gene, cluster and link-path counts.


- [ ] **Step 1: Read this first**

Planning evidence (2026-09-26): headless Chromium `--dump-dom` on the spike page gave 145 genes, 11 clusters, 668 link paths for both the 4.13 MB page and the 0.81 MB slimmed page; on real locus L001 (12 strains) 237 genes, 12 clusters, 1178 paths for both the 5.35 MB page and the 1.01 MB slimmed page; clinker 0.0.32's own scripts never read `sequence` or `translation`. The slimmed L001 page also rendered inside the island page's iframe from `file://`. This test keeps that check in the suite.


- [ ] **Step 2: Write the test**

Create `tests/test_clinker_render.py`:

```python
"""Slimmed clinker pages render the same figure as the full page.

Spec section 8: the 0.81 MB slimmed page was not yet checked to render and
behave the same. This runs clinker 0.0.32 on three small GenBank files,
slims the page with lib/clinker_html.py, loads both pages in a
headless Chromium (`--dump-dom`, which returns the DOM after scripts ran)
and compares the drawn gene, cluster and link-path counts.

Skips unless `clinker` is on PATH and a headless Chromium is found:
NOVINVENIO_HEADLESS_CHROME, else the newest
~/.cache/ms-playwright/chromium_headless_shell-*/chrome-headless-shell-linux64/chrome-headless-shell.
"""
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).parent.parent / "bin"
sys.path.insert(0, str(BIN))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from Bio import SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402
from Bio.SeqFeature import FeatureLocation, SeqFeature  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402
from clinker_html import slim_clinker_html  # noqa: E402


def find_chrome() -> str | None:
    env = os.environ.get("NOVINVENIO_HEADLESS_CHROME")
    if env and Path(env).is_file():
        return env
    hits = sorted(Path.home().glob(".cache/ms-playwright/chromium_headless_shell-*/"
                                   "chrome-headless-shell-linux64/chrome-headless-shell"))
    return str(hits[-1]) if hits else None


CHROME = find_chrome()
pytestmark = pytest.mark.skipif(not (shutil.which("clinker") and CHROME),
                                reason="needs clinker on PATH and a headless Chromium")


def write_gbks(d: Path) -> None:
    rng = random.Random(7)
    prots = ["M" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(120)) for _ in range(4)]
    for strain, keep in (("S1", [0, 1, 2, 3]), ("S2", [0, 1, 3]), ("S3", [0, 3])):
        feats, pos = [], 0
        for i in keep:
            loc = FeatureLocation(pos, pos + 3 * len(prots[i]) + 3, strand=1)
            tag = f"{strain}_g{i}"
            feats.append(SeqFeature(loc, type="gene", qualifiers={"locus_tag": [tag]}))
            feats.append(SeqFeature(loc, type="CDS", qualifiers={
                "locus_tag": [tag], "translation": [prots[i]], "note": [f"family=f{i}"]}))
            pos += 3 * len(prots[i]) + 200
        rec = SeqRecord(Seq("A" * pos), id=strain, name=strain, features=feats)
        rec.annotations["molecule_type"] = "DNA"
        SeqIO.write(rec, str(d / f"{strain}.gbk"), "genbank")
    with open(d / "groups.csv", "w") as fh:
        for strain, keep in (("S1", [0, 1, 2, 3]), ("S2", [0, 1, 3]), ("S3", [0, 3])):
            for i in keep:
                fh.write(f"{strain}_g{i},f{i}\n")


def dom_counts(page: Path) -> dict:
    dom = subprocess.run([CHROME, "--headless", "--no-sandbox", "--disable-gpu",
                          "--virtual-time-budget=10000", "--window-size=1600,1000",
                          "--dump-dom", str(page)],
                         capture_output=True, text=True, timeout=120).stdout
    return {k: dom.count(k) for k in ('class="gene"', 'class="cluster"', "<path")}


def test_slimmed_page_draws_the_same_figure(tmp_path):
    write_gbks(tmp_path)
    full = tmp_path / "full.html"
    subprocess.run(["clinker", *sorted(str(p) for p in tmp_path.glob("*.gbk")),
                    "-gf", str(tmp_path / "groups.csv"), "-j", "1", "-p", str(full)],
                   check=True, capture_output=True)
    slim = tmp_path / "slim.html"
    slim.write_text(slim_clinker_html(full.read_text()))
    assert slim.stat().st_size < full.stat().st_size
    before, after = dom_counts(full), dom_counts(slim)
    assert before['class="gene"'] == 9 and before['class="cluster"'] == 3
    assert after == before
```


- [ ] **Step 3: Run**

```bash
pixi run python -m pytest -q tests/test_clinker_render.py
```

Expected: 1 passed where `clinker` is on PATH (the pixi env) and a headless Chromium exists (`~/.cache/ms-playwright/...` or `NOVINVENIO_HEADLESS_CHROME`); otherwise 1 skipped


- [ ] **Step 4: Read this first**

Manual check, once, in a desktop browser: open a slimmed page from Task 24's output (`$SCRATCH/lv_b/clinker/L001.html`), then its unslimmed twin (`--pangenome_clinker_slim false` output or `bin/pangenome_island_clinker.py --keep_sequences`). Hover a gene (tooltip shows its name and coordinates), hover a link (identity), drag a cluster name to reorder, and click a legend circle to recolour. All four must behave the same in both pages. If any differs, set `pangenome_clinker_slim = false` in `nextflow.config` and tell the user.


- [ ] **Step 5: Commit**

```bash
git add tests/test_clinker_render.py
git commit -m "clinker panel: render check for slimmed pages (#${ISSUE_B})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 21: Page: "Synteny (clinker)" panel and CLI flags

**Files:**
- Modify: `lib/island_locus_template.py`
- Modify: `bin/pangenome_island_synteny.py`
- Modify: `tests/test_island_locus_template.py`, `tests/test_pangenome_island_synteny.py`, `tests/test_report_js_behaviour.py`, `tests/js/drive_reports.mjs`

**Interfaces:**
- Consumes: `clinker_strains` per locus (Task 16), `island_slices.tsv` (Task 18), the list of locus keys that got a page (Task 22).
- Produces: CLI flags `--slices_tsv`, `--clinker_keys`, `--clinker_enabled true|false`; `add_clinker(payload, slices_tsv, keys, enabled)`; payload `locus_meta.clinker = {enabled, keys}` and `bp_start bp_end n_genes` on each pick; JS `clinkerPanelState clinkerReason clinkerRegionText renderClinkerPanel`; DOM ids `lv-clinker lv-clinker-note lv-clinker-list lv-clinker-frame`; iframe `src="clinker/<key>.html"` only for keys matching `^L[0-9]+$`.


- [ ] **Step 1: Write the failing test**

Append to the end of `tests/test_island_locus_template.py`:

```python


# ---- clinker panel (Part B) ----

def test_clinker_panel_states():
    out = run_node(["clinkerPanelState"], """
      var on = {enabled: true, keys: ["L001"]};
      console.log(JSON.stringify([
        clinkerPanelState({key: "L001"}, on), clinkerPanelState({key: "L002"}, on),
        clinkerPanelState({key: "L001"}, {enabled: false, keys: ["L001"]}),
        clinkerPanelState({key: "../x"}, {enabled: true, keys: ["../x"]})]));""")
    assert out == [{"mode": "ok", "src": "clinker/L001.html"}, {"mode": "missing", "src": ""},
                   {"mode": "off", "src": ""}, {"mode": "missing", "src": ""}]


def test_clinker_reason_and_region_text():
    out = run_node(["locusClassLabel", "clinkerReason", "clinkerRegionText"], """
      var a = {reason: "exemplar", row_class: "full", species: "", contig: "c1", rank_lo: 3, rank_hi: 9};
      var b = {reason: "best_empty", row_class: "empty", species: "Sp two", contig: "c1",
               bp_start: 100, bp_end: 900, n_genes: 12};
      var c = {reason: "fill", row_class: "partial", species: "Sp one", contig: "c2", rank_lo: 1, rank_hi: 2};
      console.log(JSON.stringify([clinkerReason(a), clinkerReason(b), clinkerReason(c),
                                  clinkerRegionText(a), clinkerRegionText(b)]));""")
    assert out == ["locus exemplar", "best-assembled empty site strain of Sp two",
                   "more partial strains, best assembly first", "c1 gene ranks 3-9",
                   "c1:100-900, 12 genes"]


def test_clinker_panel_markup_is_below_the_grid():
    page = ISLAND_SYNTENY_TEMPLATE
    assert page.index('id="lv-grid"') < page.index('id="lv-clinker"')
    assert "Synteny (clinker)" in page
    assert "No synteny figure for this locus." in page
```

Append to the end of `tests/test_pangenome_island_synteny.py`:

```python


def test_cli_clinker_flags_fill_meta_and_region_bp(tmp_path):
    data = json.loads(json.dumps(LOCI_JSON))
    data["loci"][0]["clinker_strains"] = [
        {"strain": "S1", "reason": "exemplar", "row_class": "full", "species": "",
         "contig": "c1", "rank_lo": 0, "rank_hi": 9}]
    lj = tmp_path / "island_loci.json"
    lj.write_text(json.dumps(data))
    slices = tmp_path / "island_slices.tsv"
    slices.write_text("locus_key\tstrain\tcontig\tbp_start\tbp_end\tn_genes\tn_rescue\tn_missing\n"
                      "L001\tS1\tc1\t100\t9000\t10\t0\t0\n")
    payload = payload_of(run_cli(tmp_path, "--loci_json", str(lj), "--slices_tsv", str(slices),
                                 "--clinker_keys", "L001", "--clinker_enabled", "true"))
    assert payload["locus_meta"]["clinker"] == {"enabled": True, "keys": ["L001"]}
    pick = payload["loci"][0]["clinker_strains"][0]
    assert (pick["bp_start"], pick["bp_end"], pick["n_genes"]) == (100, 9000, 10)


def test_cli_clinker_disabled_by_default(tmp_path):
    lj = tmp_path / "island_loci.json"
    lj.write_text(json.dumps(LOCI_JSON))
    payload = payload_of(run_cli(tmp_path, "--loci_json", str(lj)))
    assert payload["locus_meta"]["clinker"] == {"enabled": False, "keys": []}
```

In `tests/test_report_js_behaviour.py`, replace this text (it occurs once):

```python
        '--project', 'demo', '--output', str(d / 'island_synteny_loci.html'))
```

with:

```python
        '--project', 'demo', '--output', str(d / 'island_synteny_loci.html'))
    run('pangenome_island_synteny.py',
        '--islands_with_domains', str(d / 'isv_islands.tsv'),
        '--presence_matrix', str(d / 'isv_matrix.tsv'),
        '--family_positions', str(d / 'isv_positions.tsv'),
        '--loci_json', str(d / 'isv_loci.json'),
        '--clinker_enabled', 'true', '--clinker_keys', 'L001',
        '--project', 'demo', '--output', str(d / 'island_synteny_clinker.html'))
```

In `tests/js/drive_reports.mjs`, replace this text (it occurs once):

```javascript
// ------------------------------- island synteny: no loci keeps the old page
```

with:

```javascript
// ------------------------------------- island synteny: clinker panel
{
  const off = boot(path.join(FX, 'island_synteny_loci.html')).window.document;
  await sleep(60);
  check('clinker panel: says the step was not run when disabled',
        /not run/.test(off.getElementById('lv-clinker-note').textContent),
        off.getElementById('lv-clinker-note').textContent);
  check('clinker panel: no iframe when disabled', !off.querySelector('#lv-clinker-frame iframe'));

  const dom = boot(path.join(FX, 'island_synteny_clinker.html'));
  const w = dom.window, d = w.document;
  const errors = [];
  w.addEventListener('error', (e) => errors.push(String(e.error)));
  await sleep(60);
  check('clinker panel: loads without error', errors.length === 0, errors.join('; '));
  const frame = () => d.querySelector('#lv-clinker-frame iframe');
  check('clinker panel: L001 loads clinker/L001.html',
        frame() && frame().getAttribute('src') === 'clinker/L001.html',
        frame() && frame().getAttribute('src'));
  check('clinker panel: lists the strains shown with a reason',
        /S1: locus exemplar/.test(d.getElementById('lv-clinker-list').textContent),
        d.getElementById('lv-clinker-list').textContent);
  [...d.querySelectorAll('#lv-list .isv-item')][1].dispatchEvent(ev(w, 'click'));
  check('clinker panel: a locus without a file says so',
        d.getElementById('lv-clinker-note').textContent === 'No synteny figure for this locus.',
        d.getElementById('lv-clinker-note').textContent);
  check('clinker panel: and shows no iframe', !frame());
}

// ------------------------------- island synteny: no loci keeps the old page
```


- [ ] **Step 2: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus_template.py tests/test_pangenome_island_synteny.py
```

Expected: FAIL: 3 template tests (`ValueError: substring not found` for `clinkerPanelState`) and 2 CLI tests (`unrecognized arguments: --slices_tsv`)


- [ ] **Step 3: Make the change**

In `lib/island_locus_template.py`, replace this text (it occurs once):

```python
  .lv-confound { margin: 0 0 12px; font-size: 12px; color: var(--text-secondary); }
```

with:

```python
  .lv-confound { margin: 0 0 12px; font-size: 12px; color: var(--text-secondary); }
  .lv-clinker { margin-top: 16px; }
  .lv-clinker h3 { margin: 0 0 6px; font-size: 14px; }
  .lv-clinker-list { margin: 0 0 10px; padding-left: 18px; font-size: 12px; color: var(--text-secondary); }
  .lv-clinker-iframe { width: 100%; height: 640px; border: 1px solid var(--border); border-radius: 8px; background: var(--surface-1); }
```

In `lib/island_locus_template.py`, replace this text (it occurs once):

```python
          <div class="isv-vscroll"><canvas class="isv-grid" id="lv-grid"></canvas></div>
        </div>
```

with:

```python
          <div class="isv-vscroll"><canvas class="isv-grid" id="lv-grid"></canvas></div>
        </div>
        <section class="lv-clinker" id="lv-clinker" aria-labelledby="lv-clinker-title">
          <h3 id="lv-clinker-title">Synteny (clinker)</h3>
          <p class="isv-main-note" id="lv-clinker-note"></p>
          <ul class="lv-clinker-list" id="lv-clinker-list"></ul>
          <div id="lv-clinker-frame"></div>
        </section>
```

In `lib/island_locus_template.py`, replace this text (it occurs once):

```python
  function lcolX(i) {
```

with:

```python
  // ---- clinker panel (spec section 8) ----
  var CLINKER = LMETA.clinker || { enabled: false, keys: [] };
  function clinkerPanelState(locus, clinker) {
    if (!clinker || !clinker.enabled) return { mode: "off", src: "" };
    if (!/^L[0-9]+$/.test(locus.key) || (clinker.keys || []).indexOf(locus.key) === -1) {
      return { mode: "missing", src: "" };
    }
    return { mode: "ok", src: "clinker/" + locus.key + ".html" };
  }
  function clinkerReason(pick) {
    if (pick.reason === "exemplar") return "locus exemplar";
    var cls = locusClassLabel(pick.row_class);
    if (pick.reason === "fill") return "more " + cls + " strains, best assembly first";
    return "best-assembled " + cls + " strain" + (pick.species ? " of " + pick.species : "");
  }
  function clinkerRegionText(pick) {
    var where = pick.bp_start
      ? pick.contig + ":" + pick.bp_start + "-" + pick.bp_end
      : pick.contig + " gene ranks " + pick.rank_lo + "-" + pick.rank_hi;
    return where + (pick.n_genes !== undefined ? ", " + pick.n_genes + " genes" : "");
  }
  function renderClinkerPanel(locus) {
    var st = clinkerPanelState(locus, CLINKER);
    var note = document.getElementById("lv-clinker-note");
    var list = document.getElementById("lv-clinker-list");
    var frame = document.getElementById("lv-clinker-frame");
    list.textContent = "";
    frame.textContent = "";
    if (st.mode === "off") {
      note.textContent = "The clinker step was not run for this page (--pangenome_clinker false).";
      return;
    }
    if (st.mode === "missing") {
      note.textContent = "No synteny figure for this locus.";
      return;
    }
    note.textContent = "clinker 0.0.32: each strain's region, genes linked by similarity and " +
      "grouped by tier-1 family (the column IDs above). Strains shown:";
    (locus.clinker_strains || []).forEach(function (pick) {
      list.appendChild(el("li", null, pick.strain + ": " + clinkerReason(pick) + "; " + clinkerRegionText(pick)));
    });
    var f = document.createElement("iframe");
    f.className = "lv-clinker-iframe";
    f.title = "clinker synteny figure for " + locus.locus_id;
    f.src = st.src;
    frame.appendChild(f);
  }

  function lcolX(i) {
```

In `bin/pangenome_island_synteny.py`, replace this text (it occurs once):

```python
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
```

with:

```python
    ap.add_argument("--slices_tsv", default=None,
                    help="island_slices.tsv from bin/pangenome_island_gbk_slice.py: bp "
                    "range and gene count of each clinker strain's region.")
    ap.add_argument("--clinker_keys", default="",
                    help="comma list of locus keys (L001,...) that have clinker/<key>.html")
    ap.add_argument("--clinker_enabled", default="false", choices=["true", "false"],
                    help="whether the clinker step ran (--pangenome_clinker)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
```

In `bin/pangenome_island_synteny.py`, replace this text (it occurs once):

```python
    add_locus_view(payload, args.loci_json, args.diagnostics_tsv)
```

with:

```python
    add_locus_view(payload, args.loci_json, args.diagnostics_tsv)
    add_clinker(payload, args.slices_tsv, args.clinker_keys, args.clinker_enabled == "true")
```

In `bin/pangenome_island_synteny.py`, replace this text (it occurs once):

```python
def main() -> int:
```

with:

```python
def add_clinker(payload: dict, slices_tsv: str | None, keys: str, enabled: bool) -> None:
    """locus_meta.clinker = {enabled, keys}; each locus's clinker_strains gain
    bp_start, bp_end and n_genes from island_slices.tsv."""
    if "loci" not in payload:
        return
    payload.setdefault("locus_meta", {})["clinker"] = {
        "enabled": enabled, "keys": sorted(k for k in keys.split(",") if k)}
    if not slices_tsv or not Path(slices_tsv).is_file() or Path(slices_tsv).stat().st_size == 0:
        return
    with open_maybe_compressed(slices_tsv) as fh:
        by_key = {(r["locus_key"], r["strain"]): r for r in csv.DictReader(fh, delimiter="\t")}
    for locus in payload["loci"]:
        for pick in locus.get("clinker_strains", []):
            row = by_key.get((locus["key"], pick["strain"]))
            if row:
                pick["bp_start"] = int(row["bp_start"])
                pick["bp_end"] = int(row["bp_end"])
                pick["n_genes"] = int(row["n_genes"])


def main() -> int:
```


- [ ] **Step 4: Run**

```bash
pixi run python -m pytest -q tests/test_island_locus_template.py tests/test_pangenome_island_synteny.py tests/test_report_templates.py tests/test_report_js_behaviour.py
```

Expected: all pass; with `NOVINVENIO_JSDOM` set the driver prints 7 `PASS clinker panel:` lines


- [ ] **Step 5: Commit**

```bash
git add lib/island_locus_template.py bin/pangenome_island_synteny.py tests/test_island_locus_template.py tests/test_pangenome_island_synteny.py tests/test_report_js_behaviour.py tests/js/drive_reports.mjs
git commit -m "clinker panel: Synteny (clinker) panel below the grid (#${ISSUE_B})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 22: Nextflow: `ISLAND_GBK_SLICE`, `ISLAND_CLINKER`, `--pangenome_clinker`

**Files:**
- Create: `modules/pangenome/island_gbk_slice.nf`, `modules/pangenome/island_clinker.nf`
- Modify: `modules/pangenome/island_loci.nf`, `modules/pangenome/island_synteny.nf`, `workflows/pangenome_profile.nf`, `nextflow.config` (params block only), `pangenome.nf` (help)
- Modify: `docs/pangenome-assumptions.md`, `CHANGES.md`, `.living/decisions.md`, `.living/learnings.md`

**Interfaces:**
- Consumes: `ISLAND_LOCI.out.regions`, `samplesheet`, `data_dir_abs`, `gff3_dir_abs`, `GENE_POSITIONS.out.positions`, `rescue_positions`, `CLUSTER_TIER1.out.cluster_tsv`.
- Produces: `ISLAND_LOCI` `emit: regions`; `ISLAND_GBK_SLICE` (`emit: locus_dirs`, optional; `emit: slices`); `ISLAND_CLINKER` (`emit: html`, optional; publishes to `pangenome/clinker/`); `ISLAND_SYNTENY` inputs `path(slices_tsv)`, `val(clinker_keys)`, `val(clinker_enabled)`; params `pangenome_clinker=true`, `pangenome_clinker_max_strains=12`, `pangenome_clinker_slim=true`, `pangenome_clinker_batch=50`.


- [ ] **Step 1: Read this first**

Nextflow 26's strict parser rejects `collate` ("Missing process or function collate") and `path(..., type: 'dir'), emit:`. The code below uses `buffer(size:, remainder: true)` and a plain `path("gbk/L*")`. During planning a toy workflow with the same channel chain gave `keys=[L001,L003]` for 3 loci with one failing, `[L001]` for 1 locus and `[]` for none, and the page process ran in all three cases.


- [ ] **Step 2: Make the change**

Create `modules/pangenome/island_gbk_slice.nf`:

```groovy
// ISLAND_GBK_SLICE -- one GenBank file per (locus, strain) region in
// island_regions.tsv, for ISLAND_CLINKER (spec section 8). One task for all
// drawn loci. Reads each chosen strain's genome FASTA, GFF3 and protein
// FASTA from the study data_dir (passed as `val`, the GENE_POSITIONS
// convention -- needs `--bind /bigdata` under singularity).
// Measured 2026-09-26 on 3 loci x 12 strains of the 529-strain run: 35 s,
// 60 MB peak RSS.
process ISLAND_GBK_SLICE {
    label 'low_cpu'
    tag "island_gbk_slice"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy',
        pattern: 'gbk/island_slices.tsv'

    input:
    path(regions)
    path(samplesheet)
    val(data_dir_abs)
    val(gff3_dir_abs)
    path(gene_positions)
    path(rescue_positions)
    path(cluster_tsv)

    output:
    path("gbk/L*"), optional: true, emit: locus_dirs
    path("gbk/island_slices.tsv"), emit: slices

    script:
    """
    pangenome_island_gbk_slice.py \
        --regions ${regions} \
        --config ${samplesheet} \
        --data_dir ${data_dir_abs} \
        --gff3_dir ${gff3_dir_abs} \
        --gene_positions ${gene_positions} \
        --rescue_positions ${rescue_positions} \
        --cluster_tsv ${cluster_tsv} \
        --id_sep '${params.pangenome_id_sep}' \
        --out_dir gbk
    """
}
```

Create `modules/pangenome/island_clinker.nf`:

```groovy
// ISLAND_CLINKER -- clinker figure per locus (spec section 8), batched
// --pangenome_clinker_batch loci per task (plan Ruling R10: one locus took
// 56 s wall on 4 cores for 12 strains x ~20 genes, 2026-09-26, so 50 loci
// fit one task of under an hour instead of 50 two-minute SLURM jobs).
//
// clinker is gamcil/clinker 0.0.32 from PyPI (pixi.toml
// [pypi-dependencies]). Bioconda's "clinker" 1.33 is an unrelated RNA-seq
// tool -- do not add it. -gf groups genes by tier-1 family, so clinker's
// groups are the grid's column IDs. bin/pangenome_island_clinker.py drops the
// embedded sequences (5.35 MB -> 1.01 MB on a real locus) unless
// --pangenome_clinker_slim false, and skips a locus whose clinker run fails
// (plan Ruling R15); the page then says "No synteny figure for this locus".
process ISLAND_CLINKER {
    label 'med_cpu'
    tag "island_clinker"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome/clinker" }, mode: 'copy'

    input:
    path(locus_dirs)

    output:
    path("L*.html"), optional: true, emit: html

    script:
    def keep = Helpers.asBool(params.pangenome_clinker_slim) ? '' : '--keep_sequences'
    """
    pangenome_island_clinker.py --locus_dirs ${locus_dirs} --cpus ${task.cpus} --out_dir . ${keep}
    """
}
```

In `modules/pangenome/island_loci.nf`, replace this text (it occurs once):

```groovy
// island_loci.json for ISLAND_SYNTENY.
//
```

with:

```groovy
// island_loci.json for ISLAND_SYNTENY. With --pangenome_clinker it also
// writes island_regions.tsv (clinker strains and rank ranges, spec section 8).
//
```

In `modules/pangenome/island_loci.nf`, replace this text (it occurs once):

```groovy
    output:
    path("island_loci.json"), emit: loci

    script:
    """
```

with:

```groovy
    output:
    path("island_loci.json"), emit: loci
    path("island_regions.tsv"), emit: regions

    script:
    def clinker_n = Helpers.asBool(params.pangenome_clinker) ? params.pangenome_clinker_max_strains : 0
    """
```

In `modules/pangenome/island_loci.nf`, replace this text (it occurs once):

```groovy
        --min_strains ${params.pangenome_top_islands_min_strains} \
        --project
```

with:

```groovy
        --min_strains ${params.pangenome_top_islands_min_strains} \
        --clinker_max_strains ${clinker_n} \
        --regions_out island_regions.tsv \
        --project
```

In `modules/pangenome/island_synteny.nf`, replace this text (it occurs once):

```groovy
    path(diagnostics_tsv)

    output:
```

with:

```groovy
    path(diagnostics_tsv)
    path(slices_tsv)
    val(clinker_keys)
    val(clinker_enabled)

    output:
```

In `modules/pangenome/island_synteny.nf`, replace this text (it occurs once):

```groovy
        --diagnostics_tsv ${diagnostics_tsv} \
        --output island_synteny.html
```

with:

```groovy
        --diagnostics_tsv ${diagnostics_tsv} \
        --slices_tsv ${slices_tsv} \
        --clinker_keys '${clinker_keys}' \
        --clinker_enabled ${clinker_enabled} \
        --output island_synteny.html
```

In `workflows/pangenome_profile.nf`, replace this text (it occurs once):

```groovy
include { ISLAND_LOCI }                                                    from '../modules/pangenome/island_loci'
```

with:

```groovy
include { ISLAND_LOCI }                                                    from '../modules/pangenome/island_loci'
include { ISLAND_GBK_SLICE }                                               from '../modules/pangenome/island_gbk_slice'
include { ISLAND_CLINKER }                                                 from '../modules/pangenome/island_clinker'
```

In `workflows/pangenome_profile.nf`, replace this text (it occurs once):

```groovy
include { EMPTY_EVALUES_STUB as EMPTY_DOMTBLOUT_STUB }           from '../modules/empty_evalues_stub'
```

with:

```groovy
include { EMPTY_EVALUES_STUB as EMPTY_DOMTBLOUT_STUB }           from '../modules/empty_evalues_stub'
include { EMPTY_EVALUES_STUB as EMPTY_ISLAND_SLICES_STUB }       from '../modules/empty_evalues_stub'
```

In `workflows/pangenome_profile.nf`, replace this text (it occurs once):

```groovy
        ISLAND_SYNTENY(
            REPORT_TABLES.out.islands_with_domains,
            rescued_matrix,
            FAMILY_POSITIONS.out.positions,
            FAMILY_PFAM_SCAN.out.domtblout,
            samplesheet,
            DIAGNOSTICS.out.banner_html,
            GENE_POSITIONS.out.positions,
            CLUSTER_TIER1.out.cluster_tsv,
            rescue_positions,
            ISLAND_LOCI.out.loci,
            DIAGNOSTICS.out.tsv,
        )
```

with:

```groovy
        // Clinker synteny panel (spec section 8). ISLAND_CLINKER tasks take
        // --pangenome_clinker_batch loci each (plan Ruling R10).
        if (Helpers.asBool(params.pangenome_clinker)) {
            ISLAND_GBK_SLICE(
                ISLAND_LOCI.out.regions, samplesheet, data_dir_abs, gff3_dir_abs,
                GENE_POSITIONS.out.positions, rescue_positions, CLUSTER_TIER1.out.cluster_tsv,
            )
            ISLAND_CLINKER(
                ISLAND_GBK_SLICE.out.locus_dirs.flatten()
                    .buffer(size: params.pangenome_clinker_batch as int, remainder: true)
            )
            clinker_keys = ISLAND_CLINKER.out.html.flatten()
                .map { html -> html.baseName }
                .collect()
                .map { keys -> keys.sort().join(',') }
                .ifEmpty('')
            island_slices = ISLAND_GBK_SLICE.out.slices
            clinker_enabled = 'true'
        }
        else {
            EMPTY_ISLAND_SLICES_STUB()
            island_slices = EMPTY_ISLAND_SLICES_STUB.out.evalues
            clinker_keys = channel.value('')
            clinker_enabled = 'false'
        }
        ISLAND_SYNTENY(
            REPORT_TABLES.out.islands_with_domains,
            rescued_matrix,
            FAMILY_POSITIONS.out.positions,
            FAMILY_PFAM_SCAN.out.domtblout,
            samplesheet,
            DIAGNOSTICS.out.banner_html,
            GENE_POSITIONS.out.positions,
            CLUSTER_TIER1.out.cluster_tsv,
            rescue_positions,
            ISLAND_LOCI.out.loci,
            DIAGNOSTICS.out.tsv,
            island_slices,
            clinker_keys,
            clinker_enabled,
        )
```

In `nextflow.config`, replace this text (it occurs once):

```groovy
    pangenome_locus_candidates  = 200            // loci whose states are computed before ranking
}
```

with:

```groovy
    pangenome_locus_candidates  = 200            // loci whose states are computed before ranking

    // Clinker synteny panel (spec section 8; ISLAND_GBK_SLICE + ISLAND_CLINKER).
    pangenome_clinker             = true
    pangenome_clinker_max_strains = 12           // strains per clinker figure (<= 66 pairs)
    pangenome_clinker_slim        = true         // drop embedded sequences from the clinker HTML
    pangenome_clinker_batch       = 50           // loci per ISLAND_CLINKER task
}
```

In `pangenome.nf`, replace this text (it occurs once):

```groovy
      --pangenome_locus_candidates     Loci scored before ranking (default: 200).
```

with:

```groovy
      --pangenome_locus_candidates     Loci scored before ranking (default: 200).
      --pangenome_clinker              Build the clinker synteny panel (default: true).
      --pangenome_clinker_max_strains  Strains per clinker figure (default: 12).
      --pangenome_clinker_slim         Remove embedded sequences from clinker pages
                                       (default: true).
      --pangenome_clinker_batch        Loci per ISLAND_CLINKER task (default: 50).
```


- [ ] **Step 3: Run**

```bash
pixi run nextflow lint modules/pangenome/island_loci.nf modules/pangenome/island_gbk_slice.nf modules/pangenome/island_clinker.nf modules/pangenome/island_synteny.nf workflows/pangenome_profile.nf pangenome.nf
```

Expected: `Nextflow linting complete!`, no errors


- [ ] **Step 4: Run**

```bash
pixi run nextflow run pangenome.nf --help | grep -E 'pangenome_clinker'
```

Expected: the four `--pangenome_clinker*` help lines


- [ ] **Step 5: Make the change**

In `docs/pangenome-assumptions.md`, replace this text (it occurs once):

```markdown
| `pangenome_locus_candidates` | `200` | `theory+practical` | Compute bound: states for 200 candidates took 1 min 41 s and 1.36 GB (2026-09-26); 1301 loci on that run. | same | Informative loci outside the top 200 by carrier proxy | Rerun with 400 and compare the drawn set |
```

with:

```markdown
| `pangenome_locus_candidates` | `200` | `theory+practical` | Compute bound: states for 200 candidates took 1 min 41 s and 1.36 GB (2026-09-26); 1301 loci on that run. | same | Informative loci outside the top 200 by carrier proxy | Rerun with 400 and compare the drawn set |
| `pangenome_clinker_max_strains` | `12` | `theory+practical` | Spec section 8: clinker compares every pair, so 12 strains = 66 pairs; one 12-strain locus took 56 s at 4 cores, 128 MB. | same | — | Time one locus at the new cap |
| `pangenome_clinker_slim` | `true` | `validated-empirically` | Removing `sequence`/`translation` shrinks a page 5.35 -> 1.01 MB; headless Chromium draws the same genes, clusters and link paths (`tests/test_clinker_render.py`). | spike + L001 | A clinker version whose scripts read these fields | Rerun `tests/test_clinker_render.py` after any clinker upgrade |
| `pangenome_clinker_batch` | `50` | `theory+practical` | HPCC job sizing (about 1 h per job): 50 loci x ~1 min at 4 cores. | same | Much larger regions or strain caps | Check ISLAND_CLINKER task time in the trace |
```

In `CHANGES.md`, replace this text (it occurs once):

```markdown
## Unreleased

```

with:

```markdown
## Unreleased

### New: clinker synteny panel

- **`lib/genbank_slice.py`**, **`bin/pangenome_island_gbk_slice.py`** (`ISLAND_GBK_SLICE`),
  **`lib/clinker_html.py`**, **`bin/pangenome_island_clinker.py`** (`ISLAND_CLINKER`) -- for each drawn
  locus, up to 12 strains' regions are written as GenBank files and drawn with gamcil/clinker 0.0.32
  (PyPI; bioconda's `clinker` is an unrelated tool), grouped by tier-1 family. Embedded sequences are
  removed (5.35 -> 1.01 MB per page). `island_synteny.html` shows the figure in a "Synteny (clinker)"
  panel. Pages publish to `pangenome/clinker/<key>.html`.
- New params: `--pangenome_clinker` (true), `--pangenome_clinker_max_strains` (12),
  `--pangenome_clinker_slim` (true), `--pangenome_clinker_batch` (50).
- The published container `0.5.0` has no clinker: rebuild and push the image before container runs.

```

Append to the end of `.living/decisions.md`:

```markdown


## 2026-09-26 — Clinker panel: batching, slimming, failure handling

**Context**: Spec section 8 asks for one ISLAND_CLINKER task per locus; one locus measured 56 s at 4 cores.
**Decision**: `--pangenome_clinker_batch` loci per task (default 50; plan Ruling R10); slim pages by default (R16); a failed locus is skipped with a warning and the page says so (R15).
**Alternatives**: one task per locus (50 two-minute SLURM jobs); fail the run on any clinker error.
**Rationale**: the HPCC job-sizing rule (about 1 h per job); one bad locus must not cost the other figures. Headless Chromium drew identical figures for slimmed and full pages.
```

Append to the end of `.living/learnings.md`:

```markdown


## 2026-09-26 — Island locus view: "empty site" at the top loci is a gene-model split

**Context**: Sequence spot check of empty-site calls on the Coccidioides run (plan Task 14), blastn of each strain's region DNA.
**Finding**: 0 of 6 empty-site calls on L001-L003 were deletions. The "empty" strains carry the locus DNA at 97.0-99.9% identity. At L001 one strain has two gene models (277 aa + 446 aa, two shell families) where the other has one 705 aa model in a third family, so the locus families are absent by clustering, not by DNA.
**Why it matters**: The informative ranking picks loci with many empty-site and full strains, which on this run selects gene-model differences. Read empty-site counts as annotation states until a DNA-level presence test exists.
**Tags**: pangenome, island-locus-view, annotation, gene-model, validation

## 2026-09-26 — Nextflow 26 strict parser: no `collate`, no `type: 'dir'` before `emit:`

**Finding**: `channel.collate(n)` fails with "Missing process or function collate"; use `buffer(size: n, remainder: true)`. `path("x*", type: 'dir'), optional: true, emit: x` fails `nextflow lint` ("Unexpected input: ','"); use `path("x*"), optional: true, emit: x`.
**Tags**: nextflow, strict-syntax, channels
```


- [ ] **Step 6: Commit**

```bash
git add modules/pangenome/island_gbk_slice.nf modules/pangenome/island_clinker.nf modules/pangenome/island_loci.nf modules/pangenome/island_synteny.nf workflows/pangenome_profile.nf nextflow.config pangenome.nf docs/pangenome-assumptions.md CHANGES.md .living/decisions.md .living/learnings.md
git commit -m "clinker panel: ISLAND_GBK_SLICE and ISLAND_CLINKER, --pangenome_clinker (#${ISSUE_B})" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 23: NovInvenio_Investigations: publish `clinker/` with the synteny page (DIFFERENT REPO)

**Files:**
- Repo: `NovInvenio_Investigations` (not NovInvenio). Use a worktree: `git -C /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations worktree add /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/.worktrees/clinker-sync -b clinker-sync origin/main`
- Modify: `bin/sync_pangenome_report.py`, `lib/pangenome_site.py` (docstring), `.gitignore`, `bin/publish_report_release.sh`, `.github/workflows/static.yml`, `tests/test_sync_pangenome_report.py`

**Interfaces:**
- Consumes: `<pangenome dir>/clinker/*.html` from Task 22.
- Produces: `docs/<domain>/<set>/<run>/clinker/` staged, gitignored (class 3, release asset), packed by `publish_report_release.sh` and merged by `static.yml`.


- [ ] **Step 1: Read this first**

All commands in this task run from the NII worktree root. NII's own CLAUDE.md applies: `clinker/` grows with locus count, so it is class 3 (release asset only, never committed).


- [ ] **Step 2: Write the failing test**

In `tests/test_sync_pangenome_report.py`, replace this text (it occurs once):

```python
        (pg / "island_synteny.html").write_text("<html>synteny</html>")
```

with:

```python
        (pg / "island_synteny.html").write_text("<html>synteny</html>")
        (pg / "clinker").mkdir()
        (pg / "clinker" / "L001.html").write_text("<html>clinker</html>")
```

In `tests/test_sync_pangenome_report.py`, replace this text (it occurs once):

```python
    ("docs/fungi/demo/run_a/island_synteny.html", True),
```

with:

```python
    ("docs/fungi/demo/run_a/island_synteny.html", True),
    ("docs/fungi/demo/run_a/clinker/L001.html", True),
```

Append to the end of `tests/test_sync_pangenome_report.py`:

```python


def test_clinker_pages_are_staged_next_to_the_synteny_viewer(tmp_path):
    _fake_study(tmp_path)
    assert _run(tmp_path, "--run", "run_a") == 0
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    assert (run_docs / "clinker" / "L001.html").read_text() == "<html>clinker</html>"


def test_resync_replaces_the_clinker_dir(tmp_path):
    study = _fake_study(tmp_path)
    _run(tmp_path, "--run", "run_a")
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    (run_docs / "clinker" / "L999.html").write_text("stale")
    pg = study / "results" / "run_a" / "output" / "pangenome"
    (pg / "clinker" / "L001.html").unlink()
    _run(tmp_path, "--run", "run_a")
    assert not (run_docs / "clinker" / "L999.html").exists()
    assert not (run_docs / "clinker" / "L001.html").exists()


def test_publish_script_and_pages_deploy_ship_clinker():
    publish = (REPO_ROOT / "bin" / "publish_report_release.sh").read_text()
    static = (REPO_ROOT / ".github" / "workflows" / "static.yml").read_text()
    assert "archive clinker island_synteny.html" in publish
    assert "for d in figures figures_pdf archive clinker; do" in static
```


- [ ] **Step 3: Run**

```bash
pixi run python -m pytest -q tests/test_sync_pangenome_report.py
```

Expected: FAIL: `test_clinker_pages_are_staged...` (FileNotFoundError), `test_gitignore_matches_data_classes[docs/fungi/demo/run_a/clinker/L001.html-True]`, `test_publish_script_and_pages_deploy_ship_clinker`


- [ ] **Step 4: Make the change**

In `bin/sync_pangenome_report.py`, replace this text (it occurs once):

```python
  2. Replace docs/<domain>/<set>/<run>/{figures,figures_pdf,archive}/ and
     island_synteny.html / assembly_quality.html with this run's copies.
```

with:

```python
  2. Replace docs/<domain>/<set>/<run>/{figures,figures_pdf,archive,clinker}/
     and island_synteny.html / assembly_quality.html with this run's copies.
     clinker/ holds the per-locus clinker pages island_synteny.html loads in
     its "Synteny (clinker)" panel (NovInvenio spec 2026-09-24 section 8).
```

In `bin/sync_pangenome_report.py`, replace this text (it occurs once):

```python
ASSET_DIRS = [("report/figures", "figures"), ("report/figures_pdf", "figures_pdf")]
RELEASE_ONLY = ["figures", "figures_pdf", "archive", "island_synteny.html", "assembly_quality.html"]
```

with:

```python
ASSET_DIRS = [("report/figures", "figures"), ("report/figures_pdf", "figures_pdf"),
              ("clinker", "clinker")]
RELEASE_ONLY = ["figures", "figures_pdf", "archive", "clinker", "island_synteny.html",
                "assembly_quality.html"]
```

In `.gitignore`, replace this text (it occurs once):

```
docs/*/*/*/island_synteny.html
```

with:

```
docs/*/*/*/island_synteny.html
docs/*/*/*/clinker/
```

In `bin/publish_report_release.sh`, replace this text (it occurs once):

```bash
             figures figures_pdf archive island_synteny.html assembly_quality.html; do
```

with:

```bash
             figures figures_pdf archive clinker island_synteny.html assembly_quality.html; do
```

In `bin/publish_report_release.sh`, replace this text (it occurs once):

```bash
# (docs/<domain>/<set>/<run>/{figures,figures_pdf,archive}/, island_synteny.html,
```

with:

```bash
# (docs/<domain>/<set>/<run>/{figures,figures_pdf,archive,clinker}/, island_synteny.html,
```

In `.github/workflows/static.yml`, replace this text (it occurs once):

```yaml
            for d in figures figures_pdf archive; do
```

with:

```yaml
            for d in figures figures_pdf archive clinker; do
```

In `lib/pangenome_site.py`, replace this text (it occurs once):

```python
  docs/<domain>/<set>/<run>/island_synteny.html  time by static.yml)
```

with:

```python
  docs/<domain>/<set>/<run>/island_synteny.html  time by static.yml)
  docs/<domain>/<set>/<run>/clinker/       } (clinker pages the synteny viewer loads)
```


- [ ] **Step 5: Run**

```bash
pixi run python -m pytest -q tests/test_sync_pangenome_report.py
```

Expected: 28 passed


- [ ] **Step 6: Commit (in the NII worktree)**

```bash
git add bin/sync_pangenome_report.py lib/pangenome_site.py .gitignore bin/publish_report_release.sh .github/workflows/static.yml tests/test_sync_pangenome_report.py
git commit -m "pangenome site: stage and publish clinker/ with island_synteny.html" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```


### Task 24: Real-data run of the clinker chain and visual checks (controller task, spec validation item 4)

**Files:**
- No repo files. Outputs go to `$SCRATCH/lv_b/`.

**Interfaces:**
- Consumes: Tasks 15-22 on the real run.
- Produces: measured cost per locus, three rendered clinker pages, the visual check notes.


- [ ] **Step 1: Run**

```bash
export PATH="$PWD/.pixi/envs/default/bin:$PATH"
ST=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome
P=$ST/results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome
O="${SCRATCH:?}/lv_b"; mkdir -p "$O/site"
python bin/pangenome_island_loci.py \
  --islands_with_domains $P/report_tables/islands_with_domains.tsv \
  --presence_matrix $P/presence_matrix.rescued.tsv --family_positions $P/family_positions.tsv.zst \
  --frequency_table $P/frequency_table.tsv --assembly_quality $P/assembly_quality_vs_content.tsv \
  --config $ST/config_immitis_in_posadasii_out.csv --top_loci 3 \
  --regions_out "$O/island_regions.tsv" --project cocci --output "$O/island_loci.json"
/usr/bin/time -v python bin/pangenome_island_gbk_slice.py --regions "$O/island_regions.tsv" \
  --config $ST/config_immitis_in_posadasii_out.csv --data_dir $ST/data_dir --gff3_dir $ST/data_dir/gff3 \
  --gene_positions $P/gene_positions.tsv.zst --rescue_positions $P/rescue_positions.tsv \
  --cluster_tsv $P/cluster/tier1_cluster.tsv --out_dir "$O/gbk" 2> "$O/gbk.time"
/usr/bin/time -v python bin/pangenome_island_clinker.py --locus_dirs "$O"/gbk/L* --cpus 4 \
  --out_dir "$O/site/clinker" 2> "$O/clinker.time"
grep -E 'slices|drawn|Elapsed|Maximum|WARNING|ERROR' "$O/gbk.time" "$O/clinker.time"
python bin/pangenome_island_synteny.py \
  --islands_with_domains $P/report_tables/islands_with_domains.tsv \
  --presence_matrix $P/presence_matrix.rescued.tsv --family_positions $P/family_positions.tsv.zst \
  --config $ST/config_immitis_in_posadasii_out.csv --loci_json "$O/island_loci.json" \
  --slices_tsv "$O/gbk/island_slices.tsv" --clinker_keys L001,L002,L003 --clinker_enabled true \
  --project cocci --output "$O/site/island_synteny.html"
ls -la "$O/site" "$O/site/clinker"
```


- [ ] **Step 2: Read this first**

Expected from planning (same code, 2026-09-26): 36 regions (12 per locus), 36 slices, no anchor mismatch, GenBank step 35 s and 60 MB; clinker on L001 alone 56 s wall at 4 cores and 128 MB, 5.35 MB page slimmed to 1.01 MB, and no "Could not find parent gene" warnings (the gene features stop them).


- [ ] **Step 3: Read this first**

Visual checks (spec validation item 4; the controller does these by eye and writes the result into the PR). Open `$O/site/island_synteny.html` in a desktop browser from `file://`. For each of L001, L002, L003:
1. The exemplar's clinker track shows the grid's columns as the same genes in the same order (the grid's column IDs are clinker's group labels; clinker may draw the track reversed, because gene_positions has no strand and the spec defers orientation to Phase 2).
2. An empty-site strain shows its left and right flanks joined with no locus genes between them. Planning note: at L001 the empty-site strains carry one merged gene model in the locus position (Task 14), so expect a gene between the flanks there; record it, do not change code.
3. The strain list under the figure names each strain, its reason and its region (contig:start-end, gene count).
Spike follow-up: for the spike island `407-0_S_OLD_CPA0002:scaffold_30:125276-129229` (3 of 5 member families co-located), rerun `pangenome_island_loci.py` with `--top_loci 200` and find the locus whose `families` contain `B11057|CCF9097_006184-T1`; report which locus columns are in place in its full-locus strains. The two near-singleton members should be `absent` or `elsewhere` there, so they do not widen any strain's region.


## Self-review (done while writing the plan)

**Spec coverage.**

| Spec item | Task |
|---|---|
| 1 Loci, not islands (50% containment, variants listed) | 1, 6 (`n_variants`), 9 (sidebar) |
| 2 Exemplar choice, `F_min` fallback, "short flanks" / "exemplar at contig end" | 2, 9 (badge) |
| 3 Columns: `F` + locus + `F`, bin, Pfam class, exemplar location, anchor marks | 3, 7, 9 |
| 4 Cell states incl. rescue hatched and contig break | 4, 9 |
| 5 Row classes, grouping by class then species, breakpoint track by species | 5, 6, 9 |
| 6 Informative ranking, alternative sorts, `--top_loci` 50 | 5, 7, 9, 12 |
| 7 Page layout, popups, diagnostics banner + confound line | 9, 10, 11 |
| 8 Clinker strains, regions, GenBank, clinker 0.0.32 PyPI, `-gf`, slimming, panel, switch, release assets | 15-24 |
| Data and wiring (inputs, streaming, payload codes) | 7, 12, 22 |
| Validation 1 (synthetic unit tests) | 1-6, 15, 17 |
| Validation 2 (regression vs feasibility numbers) | 8 |
| Validation 3 (sequence spot checks, 5 empty + 5 partial) | 14 |
| Validation 4 (clinker visual checks on 3 loci) | 24 |
| Spike: check that non-co-located members drop out of regions | 24 |
| Spike: slimmed page renders the same | 20 |

**Placeholder scan.** No "TBD", "TODO" or "similar to Task N". `ISSUE_A`/`ISSUE_B` are shell variables set in the setup steps. Every code step shows the full code.

**Type consistency.** Names were checked against the code that passed: `Placement`, `Columns`, `StrainCells`, `compute_locus` keys (`_cells`, `_pairs`, `_row_class`, `_placement`), `locus_payload`, `island_loci.json` keys, `island_regions.tsv` and `island_slices.tsv` columns, CLI flags, DOM ids and JS helper names are the same in every task that uses them.

**Review Focus.** The five items above each have a test in their owning task. Checked and left to existing tests: empty regions and `--top_loci 0` (Tasks 7, 18), missing N50 file (Task 7), unsafe strain names in file names (Task 17).

**Replay.** The code in this plan was generated from files that passed the full suite in a scratch copy (1062 passed, 8 skipped; ruff clean; `nextflow lint` clean; jsdom driver and clinker render check passed). A task-by-task replay onto a clean `origin/main` copy confirmed the red/green steps of Tasks 1-10 (test counts above are the replayed counts). The replay of Tasks 11-23 was still running when the plan was committed, so their stated pass counts come from the scratch run, not the replay; the executor should treat any count mismatch there as a plan error to report, not a code error.

## Execution

Recommended: subagent-driven. The tasks depend on each other's names and payload keys, a shipped mistake in the page or the rank rebuild is hard to see without a fresh reviewer, and Tasks 13, 14, 20 and 24 need the controller on a compute node. Tasks 13, 14 and 24 are controller tasks (real data, `$SCRATCH`); Task 23 runs in the other repo.
