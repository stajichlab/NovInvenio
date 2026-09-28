# Per-group frequency bins (#212, PR 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bin every family per `GROUP` (ingroup, outgroup) with correct `nonrep_only` / `outgroup_only` / `ingroup_only` / `absent` classes, and make the per-genome table, figure, outlier flags and composition section use each genome's own group's classes.

**Architecture:** `bin/pangenome_frequency_bins.py` gains a per-group class function and writes three `*_out` columns. `bin/pangenome_report_tables.py` reads them plus the samplesheet and strain inventory to build a per-group `per_strain_summary.tsv`. `bin/pangenome_report_render.py` draws the per-genome figure from the new columns and replaces the pie with per-group bars. Nextflow passes the new inputs and one new param.

**Tech Stack:** Python 3 (stdlib csv, numpy, matplotlib), pytest via `pixi run python -m pytest`, Nextflow DSL2.

**Spec:** `docs/superpowers/specs/2026-09-28-pangenome-group-bins-design.md`

## Global Constraints

- Class labels, exactly: `core`, `soft_core`, `shell`, `cloud`, `singleton`, `nonrep_only`, `outgroup_only`, `ingroup_only`, `absent`.
- Bin cutoffs unchanged: core 0.95, soft_core 0.90, shell 0.15 (params `pangenome_core_cutoff` / `pangenome_softcore_cutoff` / `pangenome_shell_cutoff`).
- New param `pangenome_outgroup_min_bin_strains = 3`.
- Not-binned sentinel in `*_out` columns: `-`.
- `frequency_table.tsv` keeps `family`, `frequency`, `strain_count`, `bin` with their ingroup meaning; new columns `frequency_out`, `strain_count_out`, `bin_out` are appended after them.
- Presence = `present` or `genome_only` (`PresenceMatrix.is_present`).
- Tables written before this change (no `*_out` columns) must still be read everywhere, as "outgroup not binned".
- Tests run with `pixi run python -m pytest -q <files>` from the repo root.

## Review Focus

1. Dereplication off: `strain_inventory` is the 0-byte stub `empty_evalues.tsv`. Expected: every group strain is a representative, `nonrep_only` never appears, nothing raises (`read_representative_shorts` raises on an empty file). Pinned in Task 1 and Task 3.
2. Re-render of an old run: `frequency_table.tsv` without `*_out` columns reaches REPORT_TABLES and the render. Expected: outgroup rows use `bin`, one composition bar, no crash. Pinned in Task 3 and Task 5.
3. A matrix column that is not in the samplesheet, or has a third `GROUP` value. Expected: `group` = its value (or empty), counted with `bin`, not scored with either group's median. Pinned in Task 3.
4. Outgroup present but with fewer than 3 representatives (e.g. 4 outgroup strains of which 2 are representatives). Expected: `*_out` = `-`, fallback note. Pinned in Task 1.
5. Stub-name collision (#203): REPORT_TABLES already receives several `empty_evalues.tsv` stubs; the new inventory input must not collide. Expected: staged under its own name. Pinned in Task 4 by `stageAs`, checked with `nextflow lint` and by reading the task command.

---

### Task 1: Per-group class rules and frequency table rows

**Files:**
- Modify: `bin/pangenome_frequency_bins.py`
- Test: `tests/test_pangenome_frequency_bins.py` (new)

**Interfaces:**
- Produces:
  - `group_class(count: int, n_reps: int, in_nonrep: bool, in_other_group: bool, is_ingroup: bool, core_cutoff=0.95, softcore_cutoff=0.90, shell_cutoff=0.15) -> str`
  - `GroupStrains` dataclass: `in_reps: list[str]`, `in_nonreps: list[str]`, `out_reps: list[str]`, `out_nonreps: list[str]`
  - `compute_group_frequency_table(matrix: PresenceMatrix, groups: GroupStrains, outgroup_min_bin_strains: int = 3, core_cutoff=0.95, softcore_cutoff=0.90, shell_cutoff=0.15) -> list[dict]`; each dict has `family, frequency (float), strain_count (int), bin, frequency_out (float | None), strain_count_out (int | None), bin_out (str)`; `None` / `"-"` when the outgroup is not binned.
  - `resolve_group_strains(matrix, config_path: str, inventory_path: str | None, ingroup_label: str, outgroup_label: str) -> GroupStrains`

- [ ] **Step 1: Write the failing tests**

```python
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_frequency_bins import (
    group_class, GroupStrains, compute_group_frequency_table, resolve_group_strains,
)
from pangenome_matrix import PresenceMatrix


def _matrix(tmp_path, text):
    p = tmp_path / "m.tsv"
    p.write_text(text)
    return PresenceMatrix.from_tsv(str(p))


def test_group_class_rules_in_order():
    assert group_class(10, 10, False, False, True) == "core"
    assert group_class(1, 10, False, False, True) == "singleton"
    assert group_class(0, 10, True, True, True) == "nonrep_only"
    assert group_class(0, 10, False, True, True) == "outgroup_only"
    assert group_class(0, 10, False, True, False) == "ingroup_only"
    assert group_class(0, 10, False, False, True) == "absent"


def test_group_frequency_table_classes_both_sides(tmp_path):
    # i1-i3 ingroup reps, i4 ingroup non-rep; o1-o3 outgroup reps, o4 outgroup non-rep.
    m = _matrix(tmp_path,
        "family\ti1\ti2\ti3\ti4\to1\to2\to3\to4\n"
        "fA\tpresent\tpresent\tpresent\tpresent\tpresent\tpresent\tpresent\tpresent\n"
        "fB\tabsent\tabsent\tabsent\tabsent\tpresent\tpresent\tpresent\tabsent\n"
        "fC\tabsent\tabsent\tabsent\tpresent\tabsent\tabsent\tabsent\tabsent\n"
        "fD\tabsent\tabsent\tabsent\tpresent\tpresent\tabsent\tabsent\tabsent\n"
        "fE\tabsent\tabsent\tabsent\tabsent\tabsent\tabsent\tabsent\tgenome_only\n"
        "fF\tpresent\tabsent\tabsent\tabsent\tabsent\tabsent\tabsent\tabsent\n")
    g = GroupStrains(["i1", "i2", "i3"], ["i4"], ["o1", "o2", "o3"], ["o4"])
    rows = {r["family"]: r for r in compute_group_frequency_table(m, g)}
    assert (rows["fA"]["bin"], rows["fA"]["bin_out"]) == ("core", "core")
    assert (rows["fB"]["bin"], rows["fB"]["bin_out"]) == ("outgroup_only", "core")
    assert (rows["fC"]["bin"], rows["fC"]["bin_out"]) == ("nonrep_only", "ingroup_only")
    # in an ingroup non-rep AND the outgroup: ingroup side is nonrep_only
    assert (rows["fD"]["bin"], rows["fD"]["bin_out"]) == ("nonrep_only", "singleton")
    assert (rows["fE"]["bin"], rows["fE"]["bin_out"]) == ("outgroup_only", "nonrep_only")
    assert (rows["fF"]["bin"], rows["fF"]["bin_out"]) == ("singleton", "ingroup_only")
    assert rows["fB"]["strain_count"] == 0 and rows["fB"]["strain_count_out"] == 3


def test_outgroup_below_minimum_is_not_binned(tmp_path):
    m = _matrix(tmp_path, "family\ti1\ti2\ti3\to1\to2\nfA\tpresent\tpresent\tpresent\tpresent\tabsent\n")
    g = GroupStrains(["i1", "i2", "i3"], [], ["o1", "o2"], [])
    (row,) = compute_group_frequency_table(m, g, outgroup_min_bin_strains=3)
    assert (row["frequency_out"], row["strain_count_out"], row["bin_out"]) == (None, None, "-")
    assert row["bin"] == "core"


def test_no_outgroup_is_not_binned(tmp_path):
    m = _matrix(tmp_path, "family\ti1\ti2\ti3\nfA\tabsent\tabsent\tabsent\n")
    g = GroupStrains(["i1", "i2", "i3"], [], [], [])
    (row,) = compute_group_frequency_table(m, g)
    assert row["bin_out"] == "-"
    assert row["bin"] == "absent"


def _config(tmp_path, rows):
    p = tmp_path / "c.csv"
    p.write_text("GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup\n"
                 + "".join(f"{g},sp,{s},{s}.fa,{s}.dna,{s}.gff3,{s},\n" for g, s in rows))
    return str(p)


def test_resolve_group_strains_with_inventory(tmp_path):
    m = _matrix(tmp_path, "family\ti1\ti2\to1\to2\tx1\nfA\tpresent\tpresent\tpresent\tpresent\tpresent\n")
    cfg = _config(tmp_path, [("IN", "i1"), ("IN", "i2"), ("OUT", "o1"), ("OUT", "o2"), ("NEAR_INGROUP", "x1")])
    inv = tmp_path / "inv.tsv"
    inv.write_text("Short\tdedup_group\tis_representative\ni1\t0\t1\ni2\t0\t0\no1\t1\t1\no2\t2\t1\nx1\t3\t1\n")
    g = resolve_group_strains(m, cfg, str(inv), "IN", "OUT")
    assert g == GroupStrains(["i1"], ["i2"], ["o1", "o2"], [])


def test_resolve_group_strains_empty_inventory_stub_means_no_dereplication(tmp_path):
    m = _matrix(tmp_path, "family\ti1\ti2\to1\nfA\tpresent\tpresent\tpresent\n")
    cfg = _config(tmp_path, [("IN", "i1"), ("IN", "i2"), ("OUT", "o1")])
    stub = tmp_path / "empty_evalues.tsv"
    stub.write_text("")
    g = resolve_group_strains(m, cfg, str(stub), "IN", "OUT")
    assert g == GroupStrains(["i1", "i2"], [], ["o1"], [])
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run python -m pytest -q tests/test_pangenome_frequency_bins.py`
Expected: FAIL at import (`cannot import name 'group_class'`).

- [ ] **Step 3: Implement**

In `bin/pangenome_frequency_bins.py`, add after `assign_bin` (keep `assign_bin`; replace `compute_frequency_table` and `resolve_strains`, which have no other callers — confirm with `grep -rn "compute_frequency_table\|resolve_strains" bin lib tests`):

```python
from dataclasses import dataclass, field


def group_class(
    count: int, n_reps: int, in_nonrep: bool, in_other_group: bool, is_ingroup: bool,
    core_cutoff: float = 0.95, softcore_cutoff: float = 0.90, shell_cutoff: float = 0.15,
) -> str:
    """Class of one family within one group (spec section 1), rules in order:
    counted in >= 1 representative -> the frequency bin; else present in a
    non-representative of this group -> nonrep_only; else present in the
    other group -> outgroup_only / ingroup_only; else absent."""
    if count >= 1:
        return assign_bin(count / n_reps, count, core_cutoff, softcore_cutoff, shell_cutoff)
    if in_nonrep:
        return "nonrep_only"
    if in_other_group:
        return "outgroup_only" if is_ingroup else "ingroup_only"
    return "absent"


@dataclass
class GroupStrains:
    in_reps: list[str]
    in_nonreps: list[str] = field(default_factory=list)
    out_reps: list[str] = field(default_factory=list)
    out_nonreps: list[str] = field(default_factory=list)


def compute_group_frequency_table(
    matrix: PresenceMatrix, groups: GroupStrains, outgroup_min_bin_strains: int = 3,
    core_cutoff: float = 0.95, softcore_cutoff: float = 0.90, shell_cutoff: float = 0.15,
) -> list[dict]:
    """One row per family: ingroup frequency/count/class (`frequency`,
    `strain_count`, `bin`) and, when the outgroup has at least
    `outgroup_min_bin_strains` representatives, the outgroup ones
    (`frequency_out`, `strain_count_out`, `bin_out`; None/None/"-" otherwise)."""
    cut = (core_cutoff, softcore_cutoff, shell_cutoff)
    n_in, n_out = len(groups.in_reps), len(groups.out_reps)
    bin_out_on = n_out >= outgroup_min_bin_strains
    rows = []
    for fam in matrix.families:
        def present(strains):
            return [s for s in strains if matrix.is_present(fam, s)]
        c_in = len(present(groups.in_reps))
        c_out = len(present(groups.out_reps))
        in_nonrep = bool(present(groups.in_nonreps))
        out_nonrep = bool(present(groups.out_nonreps))
        in_any = c_in > 0 or in_nonrep
        out_any = c_out > 0 or out_nonrep
        row = {
            "family": fam,
            "frequency": c_in / n_in if n_in else 0.0,
            "strain_count": c_in,
            "bin": group_class(c_in, n_in, in_nonrep, out_any, True, *cut),
            "frequency_out": None, "strain_count_out": None, "bin_out": "-",
        }
        if bin_out_on:
            row["frequency_out"] = c_out / n_out
            row["strain_count_out"] = c_out
            row["bin_out"] = group_class(c_out, n_out, out_nonrep, in_any, False, *cut)
        rows.append(row)
    return rows


def _representatives(inventory_path: str | None) -> set[str] | None:
    """Representative Shorts, or None when dereplication is off (no path, or
    the 0-byte EMPTY_INVENTORY_STUB file)."""
    if not inventory_path or Path(inventory_path).stat().st_size == 0:
        return None
    return set(read_representative_shorts(inventory_path))


def resolve_group_strains(
    matrix: PresenceMatrix, config_path: str, inventory_path: str | None,
    ingroup_label: str, outgroup_label: str,
) -> GroupStrains:
    """Ingroup/outgroup strains from the samplesheet that are matrix columns,
    split into representatives and non-representatives (spec section 1)."""
    from config_parser import parse_config  # noqa: E402

    samples = parse_config(config_path)
    in_matrix = set(matrix.strains)
    reps = _representatives(inventory_path)
    missing = [s.short for s in samples
               if s.group in (ingroup_label, outgroup_label) and s.short not in in_matrix]
    if missing:
        print(f"WARNING: {len(missing)} ingroup/outgroup strains are absent from the matrix "
              f"columns and are ignored (first few: {missing[:5]})", file=sys.stderr)

    def split(label):
        members = [s.short for s in samples if s.group == label and s.short in in_matrix]
        if reps is None:
            return members, []
        return [s for s in members if s in reps], [s for s in members if s not in reps]

    in_reps, in_nonreps = split(ingroup_label)
    out_reps, out_nonreps = split(outgroup_label)
    return GroupStrains(in_reps, in_nonreps, out_reps, out_nonreps)
```

Update the module docstring's first paragraph to say bins are computed per group and to name the new classes.

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi run python -m pytest -q tests/test_pangenome_frequency_bins.py`
Expected: 6 passed.

- [ ] **Step 5: Commit**

```bash
git add bin/pangenome_frequency_bins.py tests/test_pangenome_frequency_bins.py
git commit -m "frequency bins: per-group classes with nonrep_only/outgroup_only/ingroup_only (#212)"
```

### Task 2: Frequency-bins CLI, new columns, Nextflow param

**Files:**
- Modify: `bin/pangenome_frequency_bins.py` (`main`)
- Modify: `modules/pangenome/frequency_cooccurrence.nf` (FREQUENCY_BINS script)
- Modify: `nextflow.config` (params block, next to `pangenome_ingroup_label`)
- Test: `tests/test_pangenome_frequency_bins.py`

**Interfaces:**
- Consumes: Task 1's `resolve_group_strains`, `compute_group_frequency_table`.
- Produces: `frequency_table.tsv` header `family\tfrequency\tstrain_count\tbin\tfrequency_out\tstrain_count_out\tbin_out`; CLI flags `--outgroup_label` (default `OUT`), `--outgroup_min_bin_strains` (default 3).

- [ ] **Step 1: Write the failing test**

```python
import subprocess


def test_cli_writes_out_columns(tmp_path):
    (tmp_path / "m.tsv").write_text(
        "family\ti1\ti2\ti3\to1\to2\to3\nfA\tpresent\tpresent\tpresent\tabsent\tabsent\tpresent\n")
    cfg = _config(tmp_path, [("IN", "i1"), ("IN", "i2"), ("IN", "i3"),
                             ("OUT", "o1"), ("OUT", "o2"), ("OUT", "o3")])
    out = tmp_path / "ft.tsv"
    subprocess.run([sys.executable, str(Path(__file__).parent.parent / "bin" / "pangenome_frequency_bins.py"),
                    "--matrix", str(tmp_path / "m.tsv"), "--config", cfg,
                    "--outgroup_label", "OUT", "--outgroup_min_bin_strains", "3",
                    "--output", str(out)], check=True)
    lines = out.read_text().splitlines()
    assert lines[0] == "family\tfrequency\tstrain_count\tbin\tfrequency_out\tstrain_count_out\tbin_out"
    assert lines[1] == "fA\t1.0000\t3\tcore\t0.3333\t1\tshell"
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run python -m pytest -q tests/test_pangenome_frequency_bins.py -k cli`
Expected: FAIL (`unrecognized arguments: --outgroup_label`).

- [ ] **Step 3: Implement**

In `main()` add the two arguments and replace the strain resolution / table / writer:

```python
    ap.add_argument("--outgroup_label", default="OUT",
                    help="samplesheet GROUP value identifying outgroup strains (default: OUT)")
    ap.add_argument("--outgroup_min_bin_strains", type=int, default=3,
                    help="bin the outgroup only with at least this many representatives (default: 3)")
```

```python
    matrix = PresenceMatrix.from_tsv(args.matrix)
    groups = resolve_group_strains(matrix, args.config, args.inventory,
                                   args.ingroup_label, args.outgroup_label)
    if not groups.in_reps:
        print("ERROR: no ingroup strains left to bin over (check --config/--inventory "
              "against the matrix's columns)", file=sys.stderr)
        sys.exit(1)
    binned_out = len(groups.out_reps) >= args.outgroup_min_bin_strains
    print(f"Binning over {len(groups.in_reps)} ingroup representatives "
          f"({len(groups.in_nonreps)} non-representatives); outgroup: "
          f"{len(groups.out_reps)} representatives, "
          f"{'binned' if binned_out else 'not binned'}", file=sys.stderr)

    table = compute_group_frequency_table(
        matrix, groups, args.outgroup_min_bin_strains,
        args.core_cutoff, args.softcore_cutoff, args.shell_cutoff,
    )
    with open(args.output, "w") as fh:
        fh.write("family\tfrequency\tstrain_count\tbin\tfrequency_out\tstrain_count_out\tbin_out\n")
        for row in table:
            f_out = "-" if row["frequency_out"] is None else f"{row['frequency_out']:.4f}"
            c_out = "-" if row["strain_count_out"] is None else str(row["strain_count_out"])
            fh.write(f"{row['family']}\t{row['frequency']:.4f}\t{row['strain_count']}\t{row['bin']}\t"
                     f"{f_out}\t{c_out}\t{row['bin_out']}\n")
```

In `modules/pangenome/frequency_cooccurrence.nf` FREQUENCY_BINS script, after the `--ingroup_label` line add:

```
        --outgroup_label '${params.pangenome_outgroup_label}' \
        --outgroup_min_bin_strains ${params.pangenome_outgroup_min_bin_strains} \
```

In `nextflow.config`, after the `pangenome_outgroup_label` line:

```
    pangenome_outgroup_min_bin_strains = 3   // bin the outgroup (frequency_table *_out columns) only with >= this many representatives (#212)
```

Also add the param to `pangenome.nf`'s help text next to `--pangenome_outgroup_label` (one line, same format).

- [ ] **Step 4: Run tests and lint**

Run: `pixi run python -m pytest -q tests/test_pangenome_frequency_bins.py` → 7 passed.
Run: `nextflow lint modules/pangenome/frequency_cooccurrence.nf` → no errors.

- [ ] **Step 5: Commit**

```bash
git add bin/pangenome_frequency_bins.py modules/pangenome/frequency_cooccurrence.nf nextflow.config pangenome.nf tests/test_pangenome_frequency_bins.py
git commit -m "frequency bins: write frequency_out/strain_count_out/bin_out; pangenome_outgroup_min_bin_strains (#212)"
```

### Task 3: Per-group per-strain summary and outlier flags

**Files:**
- Modify: `bin/pangenome_report_tables.py` (`per_strain_summary`, `add_outlier_flags`, `main`)
- Test: `tests/test_pangenome_report_tables.py`

**Interfaces:**
- Consumes: frequency table columns from Task 2 (optional `bin_out`).
- Produces:
  - `PER_STRAIN_FIELDS = ["Short", "group", "is_representative", "bins_from", "n_families", "core", "soft_core", "shell", "cloud", "singleton", "nonrep_only", "outgroup_only", "singleton_z", "is_outlier"]`
  - `per_strain_summary(presence_matrix_path: str, family_bin: dict[str, str], family_bin_out: dict[str, str] | None = None, strain_group: dict[str, str] | None = None, representatives: set[str] | None = None, ingroup_label: str = "IN", outgroup_label: str = "OUT") -> list[dict]`
  - `add_outlier_flags(totals, mad_multiplier=0.6745, threshold=3.5, reference=None)`: `reference` = the rows whose `singleton` values define median/MAD (default: all of `totals`).
  - `add_group_outlier_flags(totals: list[dict]) -> list[dict]`: per `group`, reference = rows with `is_representative == "Y"`.
  - `read_strain_groups(samplesheet_path: str) -> dict[str, str]` (Short → GROUP, via `config_parser.parse_config`).
  - `read_representatives(inventory_path: str | None) -> set[str] | None` (None when no path or 0-byte file).
  - CLI: `--samplesheet`, `--strain_inventory`, `--ingroup_label` (IN), `--outgroup_label` (OUT), all optional.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_pangenome_report_tables.py`; add `add_group_outlier_flags, read_representatives` to its import list)

```python
def test_per_strain_summary_uses_own_group_classes(tmp_path):
    pm = tmp_path / "pm.tsv"
    pm.write_text("family\ti1\ti2\to1\tx1\n"
                  "fA\tpresent\tpresent\tpresent\tpresent\n"
                  "fB\tabsent\tabsent\tpresent\tabsent\n"
                  "fC\tabsent\tpresent\tabsent\tabsent\n")
    family_bin = {"fA": "core", "fB": "outgroup_only", "fC": "nonrep_only"}
    family_bin_out = {"fA": "core", "fB": "cloud", "fC": "ingroup_only"}
    rows = {r["Short"]: r for r in per_strain_summary(
        str(pm), family_bin, family_bin_out,
        strain_group={"i1": "IN", "i2": "IN", "o1": "OUT", "x1": "NEAR_INGROUP"},
        representatives={"i1", "o1", "x1"})}
    assert rows["i1"]["bins_from"] == "in" and rows["i1"]["core"] == 1
    assert rows["i2"]["is_representative"] == "N" and rows["i2"]["nonrep_only"] == 1
    assert rows["o1"]["bins_from"] == "out"
    assert (rows["o1"]["core"], rows["o1"]["cloud"], rows["o1"]["outgroup_only"]) == (1, 1, 0)
    assert rows["x1"]["group"] == "NEAR_INGROUP" and rows["x1"]["bins_from"] == "in"


def test_per_strain_summary_old_table_falls_back_to_ingroup_bins(tmp_path):
    pm = tmp_path / "pm.tsv"
    pm.write_text("family\ti1\to1\nfA\tpresent\tpresent\nfB\tabsent\tpresent\n")
    rows = {r["Short"]: r for r in per_strain_summary(
        str(pm), {"fA": "core", "fB": "outgroup_only"}, None,
        strain_group={"i1": "IN", "o1": "OUT"}, representatives=None)}
    assert rows["o1"]["bins_from"] == "in"
    assert rows["o1"]["outgroup_only"] == 1
    assert rows["o1"]["is_representative"] == "Y"


def test_per_strain_summary_absent_class_counts_only_in_total(tmp_path):
    pm = tmp_path / "pm.tsv"
    pm.write_text("family\tx1\nfZ\tpresent\n")
    (row,) = per_strain_summary(str(pm), {"fZ": "absent"}, None, strain_group={"x1": "NEAR"})
    assert row["n_families"] == 1
    assert sum(row[k] for k in ("core", "soft_core", "shell", "cloud", "singleton",
                                "nonrep_only", "outgroup_only")) == 0


def test_group_outlier_flags_use_group_representatives_only():
    rows = [{"Short": f"r{i}", "group": "IN", "is_representative": "Y", "singleton": v}
            for i, v in enumerate([10, 11, 9, 10, 12])]
    rows += [{"Short": "dup", "group": "IN", "is_representative": "N", "singleton": 400},
             {"Short": "o1", "group": "OUT", "is_representative": "Y", "singleton": 1500},
             {"Short": "o2", "group": "OUT", "is_representative": "Y", "singleton": 1400},
             {"Short": "o3", "group": "OUT", "is_representative": "Y", "singleton": 1600}]
    out = {r["Short"]: r for r in add_group_outlier_flags(rows)}
    assert out["dup"]["is_outlier"] == "Y"          # scored against IN reps
    assert out["r0"]["is_outlier"] == "N"
    assert out["o1"]["is_outlier"] == "N"           # OUT judged against OUT, not IN


def test_read_representatives_empty_stub_is_none(tmp_path):
    stub = tmp_path / "empty_evalues.tsv"
    stub.write_text("")
    assert read_representatives(str(stub)) is None
    assert read_representatives(None) is None
```

Also update the existing `test_per_strain_summary_counts_genes_and_bins` expected dicts to the new field set: add `"group": "", "is_representative": "Y", "bins_from": "in", "nonrep_only": 0, "outgroup_only": 0` to both expected rows.

- [ ] **Step 2: Run to verify they fail**

Run: `pixi run python -m pytest -q tests/test_pangenome_report_tables.py`
Expected: FAIL at import (`add_group_outlier_flags`).

- [ ] **Step 3: Implement** in `bin/pangenome_report_tables.py`

```python
PER_STRAIN_FIELDS = ["Short", "group", "is_representative", "bins_from", "n_families",
                     "core", "soft_core", "shell", "cloud", "singleton",
                     "nonrep_only", "outgroup_only", "singleton_z", "is_outlier"]
_COUNTED = ("core", "soft_core", "shell", "cloud", "singleton", "nonrep_only", "outgroup_only")


def read_strain_groups(samplesheet_path: str) -> dict[str, str]:
    from config_parser import parse_config  # noqa: E402
    return {s.short: s.group for s in parse_config(samplesheet_path)}


def read_representatives(inventory_path: str | None) -> set[str] | None:
    """Representative Shorts, or None when dereplication is off (no path or
    the 0-byte stub)."""
    if not inventory_path or Path(inventory_path).stat().st_size == 0:
        return None
    from pangenome_strain_inventory import read_representative_shorts  # noqa: E402
    return set(read_representative_shorts(inventory_path))


def per_strain_summary(
    presence_matrix_path: str, family_bin: dict[str, str],
    family_bin_out: dict[str, str] | None = None,
    strain_group: dict[str, str] | None = None,
    representatives: set[str] | None = None,
    ingroup_label: str = "IN", outgroup_label: str = "OUT",
) -> list[dict]:
    """One row per strain (spec section 3): families present, tallied by the
    strain's own group's class -- `bin_out` for outgroup strains when the
    outgroup is binned (`family_bin_out` given), `bin` otherwise.
    `genome_only` counts as present. Classes outside PER_STRAIN_FIELDS
    (`absent`) count only toward n_families."""
    strain_group = strain_group or {}
    with open(presence_matrix_path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        strains = next(reader)[1:]
        totals, lookup = {}, {}
        for s in strains:
            group = strain_group.get(s, "")
            use_out = group == outgroup_label and family_bin_out is not None
            lookup[s] = family_bin_out if use_out else family_bin
            totals[s] = {"Short": s, "group": group,
                         "is_representative": "Y" if representatives is None or s in representatives else "N",
                         "bins_from": "out" if use_out else "in", "n_families": 0,
                         **{k: 0 for k in _COUNTED}}
        for row in reader:
            family = row[0]
            for strain, call in zip(strains, row[1:]):
                if call != "absent":
                    t = totals[strain]
                    t["n_families"] += 1
                    b = lookup[strain].get(family)
                    if b in _COUNTED:
                        t[b] += 1
    return add_group_outlier_flags(list(totals.values()))
```

Change `add_outlier_flags` to take `reference`:

```python
def add_outlier_flags(totals: list[dict], mad_multiplier: float = 0.6745, threshold: float = 3.5,
                      reference: list[dict] | None = None) -> list[dict]:
    """... (keep the existing docstring) ... `reference` (default: `totals`)
    supplies the singleton values the median/MAD are computed from; every
    row of `totals` is scored against them."""
    values = [t["singleton"] for t in (totals if reference is None else reference)]
```

(The rest of the body is unchanged: it already scores `out`, built from `totals`.)

```python
def add_group_outlier_flags(totals: list[dict]) -> list[dict]:
    """Spec section 3: per `group`, median/MAD over that group's
    representatives; every strain of the group is scored."""
    out: list[dict] = []
    by_group: dict[str, list[dict]] = {}
    for t in totals:
        by_group.setdefault(t.get("group", ""), []).append(t)
    for rows in by_group.values():
        reps = [r for r in rows if r.get("is_representative", "Y") == "Y"]
        out += add_outlier_flags(rows, reference=reps)
    order = {t["Short"]: i for i, t in enumerate(totals)}
    return sorted(out, key=lambda r: order[r["Short"]])
```

In `main()`: add the four optional arguments, then replace the frequency-table read and per-strain block:

```python
    ap.add_argument("--samplesheet", default=None)
    ap.add_argument("--strain_inventory", default=None)
    ap.add_argument("--ingroup_label", default="IN")
    ap.add_argument("--outgroup_label", default="OUT")
```

```python
    family_bin: dict[str, str] = {}
    family_bin_out: dict[str, str] = {}
    with open(args.frequency_table, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            family_bin[row["family"]] = row["bin"]
            if row.get("bin_out", "-") != "-":
                family_bin_out[row["family"]] = row["bin_out"]
    strain_rows = per_strain_summary(
        args.presence_matrix, family_bin, family_bin_out or None,
        strain_group=read_strain_groups(args.samplesheet) if args.samplesheet else None,
        representatives=read_representatives(args.strain_inventory),
        ingroup_label=args.ingroup_label, outgroup_label=args.outgroup_label,
    )
    with open(out_dir / "per_strain_summary.tsv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=PER_STRAIN_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in sorted(strain_rows, key=lambda r: r["n_families"]):
            writer.writerow(row)
```

`sys.path` for `lib/` is already set in this script if it imports from lib; check the top of the file and add `sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))` only if missing.

- [ ] **Step 4: Run tests**

Run: `pixi run python -m pytest -q tests/test_pangenome_report_tables.py`
Expected: all pass (existing + 5 new).

- [ ] **Step 5: Commit**

```bash
git add bin/pangenome_report_tables.py tests/test_pangenome_report_tables.py
git commit -m "report tables: per-group per_strain_summary, within-group outlier flags (#212)"
```

### Task 4: REPORT_TABLES wiring

**Files:**
- Modify: `modules/pangenome/report.nf` (REPORT_TABLES input + script)
- Modify: `workflows/pangenome_profile.nf` (REPORT_TABLES call)

**Interfaces:**
- Consumes: Task 3 CLI flags.

- [ ] **Step 1: Edit the module.** In REPORT_TABLES `input:` after `path(gene_positions)` add:

```
    path(samplesheet)
    path(strain_inventory, stageAs: 'strain_inventory.tsv')   // 0-byte stub when dereplication is off; own name avoids the #203 stub collision
```

In its script, before `--out_dir .` add:

```
        --samplesheet ${samplesheet} \
        --strain_inventory ${strain_inventory} \
        --ingroup_label '${params.pangenome_ingroup_label}' \
        --outgroup_label '${params.pangenome_outgroup_label}' \
```

- [ ] **Step 2: Edit the workflow.** In the `REPORT_TABLES(` call add two arguments after `GENE_POSITIONS.out.positions,`:

```
        effective_samplesheet,
        strain_inventory,
```

- [ ] **Step 3: Lint**

Run: `nextflow lint modules/pangenome/report.nf workflows/pangenome_profile.nf`
Expected: no errors; warning count equal to `main` (23 on 2026-09-28 across the two runs of both files — compare with `git stash` as before).

- [ ] **Step 4: Commit**

```bash
git add modules/pangenome/report.nf workflows/pangenome_profile.nf
git commit -m "REPORT_TABLES: pass samplesheet and strain inventory (#212)"
```

### Task 5: Per-genome figure from the new columns

**Files:**
- Modify: `bin/pangenome_report_render.py` (`BAND_COLORS`, `_band_total`, `plot_per_genome_class_composition`, `render_report_markdown` note, `main` guard)
- Test: `tests/test_pangenome_report_render.py`

**Interfaces:**
- Consumes: `per_strain_summary.tsv` columns (Task 3); old files without `nonrep_only` / `is_representative` must still work.
- Produces: `GENOME_SEGMENTS = BAND_ORDER + ["nonrep_only", "outgroup_only"]`; `render_report_markdown(..., outgroup_fallback: bool = False)`.

- [ ] **Step 1: Write the failing tests**

```python
def test_per_genome_figure_marks_non_representatives(tmp_path, monkeypatch):
    labels = []
    real = pangenome_report_render.plt.Axes.set_yticklabels
    monkeypatch.setattr(pangenome_report_render.plt.Axes, "set_yticklabels",
                        lambda self, l, **k: (labels.extend(l), real(self, l, **k))[1])
    blocks = [("Sp a (IN)", [
        dict(_strain("a", core=5), is_representative="Y"),
        dict(_strain("b", core=5, nonrep_only=2), is_representative="N")])]
    plot_per_genome_class_composition(blocks, tmp_path)
    assert labels == ["a", "b †"]


def test_band_total_includes_new_segments():
    assert pangenome_report_render._band_total(
        {"core": "1", "nonrep_only": "2", "outgroup_only": "3"}) == 6


def test_report_notes_outgroup_fallback():
    md = render_report_markdown({}, {}, {}, [], 0, None, None, [100],
                                per_genome_figure=True, outgroup_fallback=True)
    assert "outgroup genomes use ingroup classes" in md
```

(`_strain` in this test file already builds rows; rows without the new keys must keep working, which the existing per-genome tests cover.)

- [ ] **Step 2: Run to verify they fail**

Run: `pixi run python -m pytest -q tests/test_pangenome_report_render.py -k "non_representatives or new_segments or fallback"`
Expected: 3 FAIL.

- [ ] **Step 3: Implement**

```python
BAND_COLORS = {
    "core": "#2c7fb8", "soft_core": "#7fcdbb", "shell": "#fed976",
    "cloud": "#fd8d3c", "singleton": "#bdbdbd",
    "nonrep_only": "#9e9ac8", "outgroup_only": "#636363",
}
GENOME_SEGMENTS = BAND_ORDER + ["nonrep_only", "outgroup_only"]
```

`_band_total` sums `GENOME_SEGMENTS`. In `plot_per_genome_class_composition`: loop over `GENOME_SEGMENTS` instead of `BAND_ORDER`, skipping a segment whose values are all 0 (so old files and runs without these classes show no empty legend entries); `legend ncol` = number of drawn segments. Genome labels:

```python
    nonrep = [r.get("is_representative", "Y") == "N" for r in rows]
    if show_labels:
        ax.set_yticks(y)
        ax.set_yticklabels([r["Short"] + (" †" if n else "") for r, n in zip(rows, nonrep)], fontsize=6)
    else:
        ax.set_yticks([])
        ax.set_ylabel(f"{len(rows)} genomes (per-genome values in per_strain_summary.tsv)")
        ys = [i for i, n in enumerate(nonrep) if n]
        if ys:
            ax.scatter([0] * len(ys), ys, marker="|", s=12, color="black",
                       transform=ax.get_yaxis_transform(), clip_on=False, zorder=3)
```

Add a legend note when any row is non-representative: append a text entry `"† / left tick: non-representative (dereplicated near-duplicate)"` via `ax.text` below the legend (`transform=ax.transAxes`, `y` just under the legend anchor).

`render_report_markdown` gets `outgroup_fallback: bool = False`; when `per_genome_figure and outgroup_fallback`, add after the figure line:

```python
        if outgroup_fallback:
            lines += ["*The outgroup has too few representative genomes to bin on its own; "
                      "outgroup genomes use ingroup classes.*", ""]
```

In `main()`: the per-genome guard becomes `all(b in per_strain_rows[0] for b in BAND_ORDER)` (unchanged; new columns are optional), and `outgroup_fallback = any(r.get("group") == args.outgroup_label and r.get("bins_from") == "in" for r in per_strain_rows)`; pass it to `render_report_markdown`.

- [ ] **Step 4: Run tests**

Run: `pixi run python -m pytest -q tests/test_pangenome_report_render.py`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add bin/pangenome_report_render.py tests/test_pangenome_report_render.py
git commit -m "per-genome figure: nonrep_only/outgroup_only segments, non-representative marks (#212)"
```

### Task 6: Per-group composition bars replace the pie

**Files:**
- Modify: `bin/pangenome_report_render.py` (`plot_frequency_bins` → `plot_group_composition`, composition section, `main` counts)
- Test: `tests/test_pangenome_report_render.py`

**Interfaces:**
- Produces:
  - `group_class_counts(frequency_table_rows: list[dict]) -> dict[str, dict[str, int]]`: keys `"ingroup"` always, `"outgroup"` when any row has `bin_out` not `-`/missing; values count every class label seen.
  - `plot_group_composition(group_counts, out_dir)` → `figures/group_composition.png` / `.pdf`.
  - `render_report_markdown(counts: dict[str, dict[str, int]], ...)`: `counts` is now the `group_class_counts` result (first positional argument; same position as before).

- [ ] **Step 1: Write the failing tests**

```python
def test_group_class_counts_both_groups():
    rows = [{"bin": "core", "bin_out": "core"}, {"bin": "outgroup_only", "bin_out": "cloud"},
            {"bin": "singleton", "bin_out": "ingroup_only"}]
    assert pangenome_report_render.group_class_counts(rows) == {
        "ingroup": {"core": 1, "outgroup_only": 1, "singleton": 1},
        "outgroup": {"core": 1, "cloud": 1, "ingroup_only": 1},
    }


def test_group_class_counts_old_table_is_ingroup_only():
    assert pangenome_report_render.group_class_counts([{"bin": "core"}]) == {"ingroup": {"core": 1}}


def test_composition_section_per_group():
    md = render_report_markdown(
        {"ingroup": {"core": 8, "singleton": 2, "outgroup_only": 5},
         "outgroup": {"core": 6, "cloud": 4, "ingroup_only": 3}},
        {}, {}, [], 0, None, None, [])
    assert "figures/group_composition.png" in md
    assert "core_shell_cloud_pie" not in md
    assert "- **core**: 8 (80.0%)" in md          # of the ingroup's 10 binned families
    assert "- **outgroup_only**: 5" in md
    assert "- **core**: 6 (60.0%)" in md
```

Update the existing tests that pass `counts={"core": 100, ...}` (e.g. `test_render_report_markdown_includes_key_sections`) to `counts={"ingroup": {...same dict...}}`; tests passing `{}` keep working.

- [ ] **Step 2: Run to verify they fail**

Run: `pixi run python -m pytest -q tests/test_pangenome_report_render.py -k "group_class_counts or per_group"`
Expected: FAIL.

- [ ] **Step 3: Implement**

```python
EXTRA_CLASSES = ["nonrep_only", "outgroup_only", "ingroup_only", "absent"]


def group_class_counts(frequency_table_rows: list[dict]) -> dict[str, dict[str, int]]:
    out: dict[str, dict[str, int]] = {"ingroup": {}}
    for r in frequency_table_rows:
        out["ingroup"][r["bin"]] = out["ingroup"].get(r["bin"], 0) + 1
        b = r.get("bin_out", "-")
        if b not in ("-", None, ""):
            out.setdefault("outgroup", {})
            out["outgroup"][b] = out["outgroup"].get(b, 0) + 1
    return out


def plot_group_composition(group_counts: dict[str, dict[str, int]], out_dir: Path) -> None:
    """One horizontal stacked bar per binned group on a shared axis (spec
    section 4); only core ... singleton are drawn."""
    groups = list(group_counts)
    fig, ax = plt.subplots(figsize=(9, 1.2 + 0.8 * len(groups)))
    left = np.zeros(len(groups))
    for band in BAND_ORDER:
        vals = np.array([group_counts[g].get(band, 0) for g in groups], dtype=float)
        ax.barh(range(len(groups)), vals, left=left, color=BAND_COLORS[band], label=band)
        left += vals
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(groups)
    ax.invert_yaxis()
    ax.set_xlabel("Gene families")
    ax.set_title("Pangenome composition by group")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.35), ncol=len(BAND_ORDER), frameon=False)
    fig.tight_layout()
    _savefig_both(fig, out_dir, "group_composition")
    plt.close(fig)
```

Delete `plot_frequency_bins`. Composition section in `render_report_markdown`:

```python
    lines += ["## Pangenome composition", ""]
    for group, gc in counts.items():
        binned = sum(gc.get(b, 0) for b in BAND_ORDER)
        lines += [f"**{group.capitalize()}** ({binned} gene families)", ""]
        for label in BAND_ORDER:
            if gc.get(label, 0):
                pct = 100 * gc[label] / binned if binned else 0
                lines.append(f"- **{label}**: {gc[label]} ({pct:.1f}%)")
        extra = [f"- **{label}**: {gc[label]}" for label in EXTRA_CLASSES if gc.get(label, 0)]
        if extra:
            lines += ["", "Not in this group's pangenome:"] + extra
        lines.append("")
    if counts:
        lines += ["![Composition](figures/group_composition.png)", ""]
    lines += ["![Frequency distribution](figures/frequency_distribution.png)", ""]
```

Remove the old `total_families` lines and the old pie line. In `main()`: replace the `counts` loop with `counts = group_class_counts(frequency_table_rows)` and `plot_frequency_bins(counts, out_dir)` with `plot_group_composition(counts, out_dir)`. Grep for any other `total_families` use in the function before deleting it.

- [ ] **Step 4: Run tests**

Run: `pixi run python -m pytest -q tests/test_pangenome_report_render.py tests/test_pangenome_report_tables.py tests/test_pangenome_frequency_bins.py`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add bin/pangenome_report_render.py tests/test_pangenome_report_render.py
git commit -m "report: per-group composition bars replace the pie (#212)"
```

### Task 7: Real-data verification and PR

**Files:** none changed unless a check fails. Scratch outputs go to the session scratchpad.

- [ ] **Step 1: Rebuild the frequency table offline** for NII `studies/fungi/coccidioides_pangenome/results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome` (`P`):

```bash
pixi run python bin/pangenome_frequency_bins.py --matrix $P/presence_matrix.rescued.tsv \
  --config $P/samplesheet.with_clades.csv --inventory $P/strain_inventory.tsv \
  --ingroup_label IN --outgroup_label OUT --outgroup_min_bin_strains 3 --output $SP/ft_new.tsv
```

Check with a short script:
- `bin` counts: `singleton` == 7,808; `outgroup_only + nonrep_only + absent` == 21,694; `nonrep_only >= 909`; `outgroup_only <= 20,785`.
- For every family with `strain_count >= 1`: `bin` equals the old table's `bin`.
- The sets of families with `bin` in {shell, cloud} and in {core, soft_core} equal the old table's sets.

- [ ] **Step 2: REPORT_TABLES + render offline** on the same run into `$SP/rt_new/` and `$SP/rr_new/` (inputs as in the 2026-09-28 re-render script, plus `--samplesheet` / `--strain_inventory` for REPORT_TABLES). Look at `per_genome_class_composition.png` and `group_composition.png`. Check: outgroup median `singleton` in `per_strain_summary.tsv` is now in the same order of magnitude as the ingroup's (report both numbers; no threshold assumed).

- [ ] **Step 3: Full test files + lint**

Run: `pixi run python -m pytest -q tests/test_pangenome_frequency_bins.py tests/test_pangenome_report_tables.py tests/test_pangenome_report_render.py`
Run: `nextflow lint modules/pangenome/*.nf workflows/pangenome_profile.nf`

- [ ] **Step 4: Re-grep consumers** of `bin` / `frequency_table` (`grep -rn '\["bin"\]\|row\["bin"\]\|frequency_table' bin lib workflows modules`) and confirm the spec's consumer table is complete.

- [ ] **Step 5: Push and open the PR** referencing #212 and the spec, listing the real-data numbers from Steps 1–2 and that no Nextflow test run was done unless one was.
