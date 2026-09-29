# Ingroup-outgroup overlap figures (#212, PR 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `report_tables/group_class_overlap.tsv` and three figures (class x class heatmap, UpSet plot, three-way bar) plus a report section showing how the ingroup's and outgroup's pangenomes overlap.

**Architecture:** REPORT_TABLES collapses each family's `bin` / `bin_out` into 7 side values and writes a 49-cell long table (empty file when the outgroup is not binned). REPORT_RENDER reads that table, draws the three figures with matplotlib and adds one report section. Nextflow passes the table from REPORT_TABLES to REPORT_RENDER.

**Tech Stack:** Python 3 (csv, numpy, matplotlib), pytest via `pixi run python -m pytest`, Nextflow DSL2.

**Spec:** `docs/superpowers/specs/2026-09-28-pangenome-group-overlap-figures-design.md`

## Global Constraints

- Side values, in this order: `core`, `soft_core`, `shell`, `cloud`, `singleton`, `nonrep_only`, `absent`.
- Collapse: ingroup `outgroup_only` / `absent` -> `absent`; outgroup `ingroup_only` / `absent` -> `absent`.
- Table header exactly `ingroup_class\toutgroup_class\tn_families`, then 49 data rows in row-major order of the list above.
- Skip rule: no `bin_out` column, or every `bin_out` is `-` -> write the table as a 0-byte file; render skips the section with the reason.
- Figure names: `group_class_overlap_heatmap`, `group_class_overlap_upset`, `group_class_overlap_shared` (PNG + PDF via `_savefig_both`).
- No new dependency (matplotlib only).

## Review Focus

1. Old table (no `bin_out`) or outgroup not binned: REPORT_TABLES writes a 0-byte file and REPORT_RENDER must not crash on it. Pinned in Task 1 and Task 4.
2. Samplesheet with several species in one group (FSSC_ambrosia ingroup has 16): axis label falls back to "Ingroup (IN)". Pinned in Task 2.
3. A cell of 0 in the heatmap: blank, not "0" text on a coloured square, and log colour must not warn on log10(0). Pinned in Task 2.
4. REPORT_RENDER on a run without the new input (manual re-render, old pipeline outputs): `--group_class_overlap` omitted -> section skipped silently. Pinned in Task 4.
5. Nextflow: the new REPORT_RENDER input must be a real file on every run (REPORT_TABLES always emits it, possibly empty). Checked in Task 5 by lint and by reading the emit.

---

### Task 1: Overlap table in REPORT_TABLES

**Files:** Modify `bin/pangenome_report_tables.py`; Test `tests/test_pangenome_report_tables.py`

**Interfaces — Produces:**
- `OVERLAP_CLASSES = ["core", "soft_core", "shell", "cloud", "singleton", "nonrep_only", "absent"]`
- `overlap_side(label: str) -> str`
- `group_class_overlap(frequency_rows: list[dict]) -> list[tuple[str, str, int]] | None` (None = skip rule)
- CLI: REPORT_TABLES always writes `group_class_overlap.tsv` (header + 49 rows, or 0 bytes).

- [ ] **Step 1: failing tests** (append; add `group_class_overlap, overlap_side` to the import list)

```python
def test_overlap_side_collapses_other_group_labels():
    assert overlap_side("outgroup_only") == "absent"
    assert overlap_side("ingroup_only") == "absent"
    assert overlap_side("absent") == "absent"
    assert overlap_side("nonrep_only") == "nonrep_only"
    assert overlap_side("core") == "core"


def test_group_class_overlap_counts_all_49_cells():
    rows = [{"bin": "core", "bin_out": "core"}, {"bin": "core", "bin_out": "ingroup_only"},
            {"bin": "outgroup_only", "bin_out": "cloud"}, {"bin": "nonrep_only", "bin_out": "nonrep_only"}]
    cells = group_class_overlap(rows)
    assert len(cells) == 49
    d = {(a, b): n for a, b, n in cells}
    assert d[("core", "core")] == 1 and d[("core", "absent")] == 1
    assert d[("absent", "cloud")] == 1 and d[("nonrep_only", "nonrep_only")] == 1
    assert sum(d.values()) == 4
    assert cells[0][:2] == ("core", "core") and cells[-1][:2] == ("absent", "absent")


def test_group_class_overlap_skip_rule():
    assert group_class_overlap([{"bin": "core"}]) is None
    assert group_class_overlap([{"bin": "core", "bin_out": "-"}]) is None
```

- [ ] **Step 2: run, expect FAIL** — `pixi run python -m pytest -q tests/test_pangenome_report_tables.py` → ImportError.

- [ ] **Step 3: implement** (near `PER_STRAIN_FIELDS`)

```python
OVERLAP_CLASSES = ["core", "soft_core", "shell", "cloud", "singleton", "nonrep_only", "absent"]


def overlap_side(label: str) -> str:
    """Collapse one group's class for the overlap table (PR 2 spec): the
    other-group-only labels and `absent` all mean 'not in this group'."""
    return "absent" if label in ("outgroup_only", "ingroup_only", "absent") else label


def group_class_overlap(frequency_rows: list[dict]) -> list[tuple[str, str, int]] | None:
    """(ingroup value, outgroup value, n_families) for all 49 cells in
    OVERLAP_CLASSES row-major order; None when the outgroup is not binned
    (no `bin_out` column, or every value is `-`)."""
    if not frequency_rows or "bin_out" not in frequency_rows[0]:
        return None
    if all(r["bin_out"] == "-" for r in frequency_rows):
        return None
    counts: dict[tuple[str, str], int] = {}
    for r in frequency_rows:
        key = (overlap_side(r["bin"]), overlap_side(r["bin_out"]))
        counts[key] = counts.get(key, 0) + 1
    return [(a, b, counts.get((a, b), 0)) for a in OVERLAP_CLASSES for b in OVERLAP_CLASSES]
```

In `main()`, keep the frequency-table rows while building `family_bin` (collect `freq_rows = list(csv.DictReader(...))`, then build the two dicts from `freq_rows`), and after writing `per_strain_summary.tsv`:

```python
    overlap = group_class_overlap(freq_rows)
    with open(out_dir / "group_class_overlap.tsv", "w") as out:
        if overlap is not None:
            out.write("ingroup_class\toutgroup_class\tn_families\n")
            for a, b, n in overlap:
                out.write(f"{a}\t{b}\t{n}\n")
```

- [ ] **Step 4: run, expect PASS.**
- [ ] **Step 5: commit** `report tables: group_class_overlap.tsv (#212 PR 2)`.

### Task 2: Heatmap and axis labels

**Files:** Modify `bin/pangenome_report_render.py`; Test `tests/test_pangenome_report_render.py`

**Interfaces — Produces:**
- `read_group_class_overlap(path: str | None) -> dict[tuple[str, str], int] | None` (None for missing path or 0-byte file)
- `group_axis_label(samplesheet: dict[str, tuple[str, str]] | None, label: str, fallback: str) -> str`
- `plot_group_class_overlap_heatmap(cells: dict[tuple[str, str], int], in_label: str, out_label: str, out_dir: Path) -> None`
- `OVERLAP_CLASSES` (same list as Task 1, defined in render too — the scripts do not import each other)

- [ ] **Step 1: failing tests**

```python
def test_read_group_class_overlap(tmp_path):
    p = tmp_path / "o.tsv"
    p.write_text("ingroup_class\toutgroup_class\tn_families\ncore\tcore\t5\ncore\tabsent\t0\n")
    assert pangenome_report_render.read_group_class_overlap(str(p)) == {("core", "core"): 5, ("core", "absent"): 0}
    (tmp_path / "e.tsv").write_text("")
    assert pangenome_report_render.read_group_class_overlap(str(tmp_path / "e.tsv")) is None
    assert pangenome_report_render.read_group_class_overlap(None) is None


def test_group_axis_label_single_vs_multi_species():
    ss = {"a": ("IN", "Sp one"), "b": ("IN", "Sp one"), "c": ("OUT", "Sp two"), "d": ("OUT", "Sp three")}
    f = pangenome_report_render.group_axis_label
    assert f(ss, "IN", "Ingroup") == "Sp one (IN)"
    assert f(ss, "OUT", "Outgroup") == "Outgroup (OUT)"
    assert f(None, "IN", "Ingroup") == "Ingroup (IN)"


def test_heatmap_writes_files_and_handles_zero_cells(tmp_path):
    import warnings
    C = pangenome_report_render.OVERLAP_CLASSES
    cells = {(a, b): 0 for a in C for b in C}
    cells[("core", "core")] = 10
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        pangenome_report_render.plot_group_class_overlap_heatmap(cells, "Sp one (IN)", "Sp two (OUT)", tmp_path)
    assert (tmp_path / "figures" / "group_class_overlap_heatmap.png").stat().st_size > 0
    assert (tmp_path / "figures_pdf" / "group_class_overlap_heatmap.pdf").exists()
```

- [ ] **Step 2: run, expect FAIL.**
- [ ] **Step 3: implement**

```python
OVERLAP_CLASSES = ["core", "soft_core", "shell", "cloud", "singleton", "nonrep_only", "absent"]


def read_group_class_overlap(path: str | None) -> dict[tuple[str, str], int] | None:
    if not path or Path(path).stat().st_size == 0:
        return None
    with open(path, newline="") as fh:
        return {(r["ingroup_class"], r["outgroup_class"]): int(r["n_families"])
                for r in csv.DictReader(fh, delimiter="\t")}


def group_axis_label(samplesheet, label: str, fallback: str) -> str:
    """'<Species> (<GROUP>)' when every samplesheet strain of the group has
    one Species, else '<fallback> (<GROUP>)'."""
    species = {sp for grp, sp in (samplesheet or {}).values() if grp == label}
    name = species.pop() if len(species) == 1 else fallback
    return f"{name} ({label})"


def plot_group_class_overlap_heatmap(cells, in_label: str, out_label: str, out_dir: Path) -> None:
    n = len(OVERLAP_CLASSES)
    grid = np.array([[cells.get((a, b), 0) for b in OVERLAP_CLASSES] for a in OVERLAP_CLASSES], dtype=float)
    shown = np.ma.masked_where(grid == 0, np.log10(grid + 1))
    fig, ax = plt.subplots(figsize=(8, 6.5))
    im = ax.imshow(shown, cmap="Blues")
    for i in range(n):
        for j in range(n):
            if grid[i, j]:
                dark = shown[i, j] > 0.6 * shown.max()
                ax.text(j, i, f"{int(grid[i, j])}", ha="center", va="center", fontsize=7,
                        color="white" if dark else "black")
    ax.set_xticks(range(n)); ax.set_xticklabels(OVERLAP_CLASSES, rotation=45, ha="right")
    ax.set_yticks(range(n)); ax.set_yticklabels(OVERLAP_CLASSES)
    ax.set_xlabel(out_label); ax.set_ylabel(in_label)
    ax.set_title("Gene families by class in each group")
    fig.colorbar(im, ax=ax, label="log10(families + 1)")
    fig.tight_layout()
    _savefig_both(fig, out_dir, "group_class_overlap_heatmap")
    plt.close(fig)
```


- [ ] **Step 4: run, expect PASS.** - [ ] **Step 5: commit** `report: group class overlap heatmap (#212 PR 2)`.

### Task 3: Three-way bar and UpSet plot

**Files:** Modify `bin/pangenome_report_render.py`; Test `tests/test_pangenome_report_render.py`

**Interfaces — Produces:**
- `overlap_intersections(cells) -> list[tuple[str | None, str | None, int]]`: `(in_class, out_class, n)` for every non-zero cell except (absent, absent); `None` for an absent side; sorted by `n` descending, ties by OVERLAP_CLASSES order.
- `shared_split(cells) -> dict[str, dict[str, int]]`: keys `shared`, `ingroup_only`, `outgroup_only`; inner keys are classes (ingroup class for shared / ingroup_only, outgroup class for outgroup_only).
- `plot_group_class_overlap_shared(cells, in_label, out_label, out_dir)`, `plot_group_class_overlap_upset(cells, in_label, out_label, out_dir)`.

- [ ] **Step 1: failing tests**

```python
def _cells(**kw):
    C = pangenome_report_render.OVERLAP_CLASSES
    d = {(a, b): 0 for a in C for b in C}
    for k, v in kw.items():
        a, b = k.split("__")
        d[(a, b)] = v
    return d


def test_overlap_intersections_pairs_and_singles():
    cells = _cells(core__core=5, core__absent=3, absent__cloud=7, absent__absent=0)
    assert pangenome_report_render.overlap_intersections(cells) == [
        (None, "cloud", 7), ("core", "core", 5), ("core", None, 3)]


def test_shared_split_totals():
    cells = _cells(core__core=5, shell__cloud=2, core__absent=3, absent__cloud=7, nonrep_only__absent=1)
    s = pangenome_report_render.shared_split(cells)
    assert s == {"shared": {"core": 5, "shell": 2}, "ingroup_only": {"core": 3, "nonrep_only": 1},
                 "outgroup_only": {"cloud": 7}}


def test_shared_and_upset_write_files(tmp_path):
    cells = _cells(core__core=5, core__absent=3, absent__cloud=7, nonrep_only__nonrep_only=1)
    pangenome_report_render.plot_group_class_overlap_shared(cells, "A (IN)", "B (OUT)", tmp_path)
    pangenome_report_render.plot_group_class_overlap_upset(cells, "A (IN)", "B (OUT)", tmp_path)
    for name in ("group_class_overlap_shared", "group_class_overlap_upset"):
        assert (tmp_path / "figures" / f"{name}.png").stat().st_size > 0
```

- [ ] **Step 2: run, expect FAIL.**
- [ ] **Step 3: implement**

```python
def overlap_intersections(cells):
    order = {c: i for i, c in enumerate(OVERLAP_CLASSES)}
    out = []
    for (a, b), n in cells.items():
        if n and (a, b) != ("absent", "absent"):
            out.append((None if a == "absent" else a, None if b == "absent" else b, n))
    return sorted(out, key=lambda t: (-t[2], order.get(t[0] or "absent"), order.get(t[1] or "absent")))


def shared_split(cells):
    out = {"shared": {}, "ingroup_only": {}, "outgroup_only": {}}
    for (a, b), n in cells.items():
        if not n or (a, b) == ("absent", "absent"):
            continue
        if a != "absent" and b != "absent":
            key, cls = "shared", a
        elif b == "absent":
            key, cls = "ingroup_only", a
        else:
            key, cls = "outgroup_only", b
        out[key][cls] = out[key].get(cls, 0) + n
    return out


def plot_group_class_overlap_shared(cells, in_label: str, out_label: str, out_dir: Path) -> None:
    """Spec 2c: shared (by ingroup class), ingroup only, outgroup only."""
    split = shared_split(cells)
    rows = [("shared", f"shared (classes: {in_label})"), ("ingroup_only", f"only {in_label}"),
            ("outgroup_only", f"only {out_label}")]
    fig, ax = plt.subplots(figsize=(9, 3.2))
    for i, (key, _) in enumerate(rows):
        left = 0
        for cls in OVERLAP_CLASSES[:-1]:
            v = split[key].get(cls, 0)
            if v:
                ax.barh(i, v, left=left, color=BAND_COLORS[cls],
                        label=cls if cls not in ax.get_legend_handles_labels()[1] else None)
                left += v
        ax.text(left, i, f" {left}", va="center", fontsize=8)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[1] for r in rows])
    ax.invert_yaxis()
    ax.set_xlabel("Gene families")
    ax.set_title("Shared and group-specific gene families")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=6, frameon=False)
    fig.tight_layout()
    _savefig_both(fig, out_dir, "group_class_overlap_shared")
    plt.close(fig)


def plot_group_class_overlap_upset(cells, in_label: str, out_label: str, out_dir: Path) -> None:
    """Spec 2b: intersection sizes over a dot matrix of the 12 group-class
    sets (6 per group; `absent` is not a set)."""
    inter = overlap_intersections(cells)
    sets = [("IN", c) for c in OVERLAP_CLASSES[:-1]] + [("OUT", c) for c in OVERLAP_CLASSES[:-1]]
    set_names = [f"{in_label if g == 'IN' else out_label}: {c}" for g, c in sets]
    fig = plt.figure(figsize=(max(8, 0.28 * len(inter) + 4), 7))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 4], height_ratios=[2, 1.6], hspace=0.05, wspace=0.02)
    ax_bar = fig.add_subplot(gs[0, 1]); ax_dot = fig.add_subplot(gs[1, 1], sharex=ax_bar)
    ax_set = fig.add_subplot(gs[1, 0], sharey=ax_dot)
    x = np.arange(len(inter))
    ax_bar.bar(x, [n for _, _, n in inter], color="#444444")
    ax_bar.set_ylabel("Families"); ax_bar.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_bar.set_title("Gene families by group-class intersection")
    ys = {s: i for i, s in enumerate(sets)}
    ax_dot.scatter(np.repeat(x, len(sets)), np.tile(range(len(sets)), len(x)), s=10, color="#dddddd")
    for xi, (a, b, _) in enumerate(inter):
        members = [ys[("IN", a)]] if a else []
        members += [ys[("OUT", b)]] if b else []
        ax_dot.plot([xi] * len(members), members, "-o", color="#222222", markersize=4)
    ax_dot.set_yticks(range(len(sets))); ax_dot.set_yticklabels([])
    ax_dot.invert_yaxis(); ax_dot.set_xticks([])
    size = [sum(n for (a, b), n in cells.items() if (g == "IN" and a == c) or (g == "OUT" and b == c))
            for g, c in sets]
    ax_set.barh(range(len(sets)), size, color=[BAND_COLORS[c] for _, c in sets])
    ax_set.invert_xaxis(); ax_set.set_yticks(range(len(sets))); ax_set.set_yticklabels(set_names, fontsize=7)
    ax_set.set_xlabel("Set size")
    _savefig_both(fig, out_dir, "group_class_overlap_upset")
    plt.close(fig)
```

- [ ] **Step 4: run, expect PASS.** - [ ] **Step 5: commit** `report: overlap UpSet and shared/specific bar (#212 PR 2)`.

### Task 4: Report section and render CLI

**Files:** Modify `bin/pangenome_report_render.py`; Test `tests/test_pangenome_report_render.py`

**Interfaces — Produces:** `render_report_markdown(..., overlap_section: str = "")` where `overlap_section` is `""` (no section), `"ok"` or a skip reason `"no_out_columns"` / `"not_binned"`; CLI `--group_class_overlap` (optional).

- [ ] **Step 1: failing tests**

```python
def test_overlap_section_ok_and_skip():
    md = render_report_markdown({}, {}, {}, [], 0, None, None, [], overlap_section="ok")
    assert "## Ingroup vs outgroup content" in md
    for name in ("group_class_overlap_heatmap", "group_class_overlap_shared", "group_class_overlap_upset"):
        assert f"figures/{name}.png" in md
    assert "report_tables/group_class_overlap.tsv" in md
    md = render_report_markdown({}, {}, {}, [], 0, None, None, [], overlap_section="not_binned")
    assert "## Ingroup vs outgroup content" in md and "not binned" in md
    assert "group_class_overlap_heatmap" not in md
    md = render_report_markdown({}, {}, {}, [], 0, None, None, [])
    assert "Ingroup vs outgroup content" not in md
```

- [ ] **Step 2: run, expect FAIL.**
- [ ] **Step 3: implement.** After the per-genome block in `render_report_markdown`:

```python
    if overlap_section:
        lines += ["## Ingroup vs outgroup content", ""]
        if overlap_section == "ok":
            lines += ["Each gene family's class in the ingroup against its class in the outgroup "
                      "(counts; colour is log-scaled).", "",
                      "![Class overlap](figures/group_class_overlap_heatmap.png)", "",
                      "Families in both groups, split by ingroup class, and families found in one group only.", "",
                      "![Shared and group-specific families](figures/group_class_overlap_shared.png)", "",
                      "Every non-empty group-class combination, largest first.", "",
                      "![Class intersections](figures/group_class_overlap_upset.png)", "",
                      "Counts: `report_tables/group_class_overlap.tsv`.", ""]
        elif overlap_section == "no_out_columns":
            lines += ["*Not shown: this run's frequency table predates per-group bins.*", ""]
        else:
            lines += ["*Not shown: the outgroup is not binned (too few representative genomes).*", ""]
```

In `main()`: add `ap.add_argument("--group_class_overlap", default=None)`. After the per-genome figure:

```python
    overlap_cells = read_group_class_overlap(args.group_class_overlap)
    if overlap_cells is not None:
        ss = read_samplesheet_groups(args.samplesheet) if args.samplesheet else None
        in_lab = group_axis_label(ss, args.ingroup_label, "Ingroup")
        out_lab = group_axis_label(ss, args.outgroup_label, "Outgroup")
        plot_group_class_overlap_heatmap(overlap_cells, in_lab, out_lab, out_dir)
        plot_group_class_overlap_shared(overlap_cells, in_lab, out_lab, out_dir)
        plot_group_class_overlap_upset(overlap_cells, in_lab, out_lab, out_dir)
        overlap_section = "ok"
    elif args.group_class_overlap is None:
        overlap_section = ""
    elif frequency_table_rows and "bin_out" not in frequency_table_rows[0]:
        overlap_section = "no_out_columns"
    else:
        overlap_section = "not_binned"
```

and pass `overlap_section=overlap_section` to `render_report_markdown`.

- [ ] **Step 4: run the three test files, expect PASS.** - [ ] **Step 5: commit** `report: ingroup vs outgroup content section (#212 PR 2)`.

### Task 5: Nextflow wiring

**Files:** `modules/pangenome/report.nf`, `workflows/pangenome_profile.nf`

- [ ] REPORT_TABLES `output:` add `path("group_class_overlap.tsv"), emit: group_class_overlap`.
- [ ] REPORT_RENDER `input:` add `path(group_class_overlap)` after `path(samplesheet)`; script: `--group_class_overlap ${group_class_overlap} \` before `--out_dir report`.
- [ ] Workflow REPORT_RENDER call: add `REPORT_TABLES.out.group_class_overlap,` after `samplesheet,`.
- [ ] `nextflow lint modules/pangenome/report.nf workflows/pangenome_profile.nf` → no errors, warning count unchanged from `main`.
- [ ] Commit `REPORT_TABLES -> REPORT_RENDER: group_class_overlap (#212 PR 2)`.

### Task 6: Real-data check, docs

- [ ] Run REPORT_TABLES + render offline for NII `rescue_freqpol_immitis_in_posadasii_out` (current outputs are PR 1's) into scratch; the 49-cell table must equal the spec's example table exactly.
- [ ] Look at the three PNGs.
- [ ] `CHANGES.md` Unreleased: one entry for `group_class_overlap.tsv` + three figures + section.
- [ ] Full three test files; commit.
