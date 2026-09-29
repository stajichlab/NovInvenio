#!/usr/bin/env python3
"""Figures (matplotlib, PNG+PDF per figure) + templated Markdown report for
the pangenome island+Pfam enrichment step. Fully automated -- no
hand-written narrative -- so it works for any future pangenome study, not
just Afumigatus (whose REPORT.md was hand-assembled prose; this is not
that). Figure functions generalize
studies/fungi/Afumigatus_pangenome/bin/plot_pangenome_summary.py's five
figures (frequency histogram, composition pie, presence/absence raster,
accumulation/rarefaction curve, classification bar), plus two new figures
(top enriched Pfam domains, which Afumigatus embedded by hand with no
checked-in generator) and a Heaps'-law openness fit + core-genome
asymptote fit that Afumigatus's own report never computed (it only
asserted openness from the curve's visual shape).

--islands_with_domains/--island_size_distribution/--island_pfam_enrichment/
--marker_summary are optional (issue #135): the core figures/sections
(frequency distribution, presence/absence heatmap, accumulation curve,
classification counts) render on every run regardless, using only
--frequency_table/--presence_matrix/--classification_counts/
--per_strain_summary, none of which need the `--pangenome_island_pfam_hmm`
accessory-island branch. Omitting the islands-only args (or letting them
point at an empty/0-byte file, as REPORT_RENDER's Nextflow wiring does when
that branch didn't run) simply skips the "Accessory islands"/"Pfam domain
enrichment"/"Marker co-occurrence" sections instead of erroring.

Usage (full report, islands + Pfam branch enabled):
  pangenome_report_render.py --frequency_table frequency_table.tsv \\
      --presence_matrix presence_matrix.tsv \\
      --islands_with_domains islands_with_domains.tsv \\
      --island_size_distribution island_size_distribution.tsv \\
      --classification_counts classification_counts.tsv \\
      --island_pfam_enrichment island_pfam_enrichment.tsv \\
      --per_strain_summary per_strain_summary.tsv \\
      --n_permutations 20 --seed 0 \\
      --out_dir report/

Usage (core report only, no islands/Pfam branch):
  pangenome_report_render.py --frequency_table frequency_table.tsv \\
      --presence_matrix presence_matrix.tsv \\
      --classification_counts classification_counts.tsv \\
      --per_strain_summary per_strain_summary.tsv \\
      --n_permutations 20 --seed 0 \\
      --out_dir report/
"""
from __future__ import annotations

import argparse
import csv
import random
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import PresenceMatrix, PRESENT, GENOME_ONLY  # noqa: E402

BAND_ORDER = ["core", "soft_core", "shell", "cloud", "singleton"]
BAND_COLORS = {
    "core": "#2c7fb8", "soft_core": "#7fcdbb", "shell": "#fed976",
    "cloud": "#fd8d3c", "singleton": "#bdbdbd",
    "nonrep_only": "#9e9ac8", "outgroup_only": "#636363",
}
GENOME_SEGMENTS = BAND_ORDER + ["nonrep_only", "outgroup_only"]
"""Per-genome figure segments: the five bands plus the #212 classes a genome
can carry that are outside its group's bands."""


def _savefig_both(fig, out_dir: Path, name: str) -> None:
    """Writes both <name>.png (into out_dir/figures/) and <name>.pdf
    (into out_dir/figures_pdf/) -- matches Afumigatus's convention of a
    PNG for the Markdown embed and a parallel PDF vector copy for print."""
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)
    (out_dir / "figures_pdf").mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "figures" / f"{name}.png", dpi=150)
    fig.savefig(out_dir / "figures_pdf" / f"{name}.pdf")


def plot_frequency_distribution(frequency_table_rows: list[dict], out_dir: Path) -> None:
    """Ported from plot_pangenome_summary.py's plot_frequency_distribution --
    family-frequency histogram, one series per band."""
    fig, ax = plt.subplots(figsize=(8, 5))
    for band in BAND_ORDER:
        freqs = [float(r["frequency"]) for r in frequency_table_rows if r["bin"] == band]
        if freqs:
            ax.hist(freqs, bins=40, range=(0, 1), color=BAND_COLORS[band],
                     label=f"{band} (n={len(freqs)})", alpha=0.85)
    ax.set_xlabel("Fraction of strains carrying the family")
    ax.set_ylabel("Number of gene families")
    ax.set_title("Family-frequency distribution")
    ax.legend()
    fig.tight_layout()
    _savefig_both(fig, out_dir, "frequency_distribution")
    plt.close(fig)


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


def outgroup_fallback_reason(frequency_table_rows: list[dict], per_strain_rows: list[dict],
                             outgroup_label: str) -> str:
    """Why outgroup genomes are drawn with ingroup classes: "" (they are not),
    "no_out_columns" (the frequency table predates per-group bins) or
    "too_few" (the outgroup had too few representatives to bin)."""
    if not any(r.get("group") == outgroup_label and r.get("bins_from") == "in"
               for r in per_strain_rows):
        return ""
    if frequency_table_rows and "bin_out" not in frequency_table_rows[0]:
        return "no_out_columns"
    return "too_few"


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
    species = {sp.strip() for grp, sp in (samplesheet or {}).values()
               if grp.strip() == label and sp.strip()}
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
    ax.set_xticks(range(n))
    ax.set_xticklabels(OVERLAP_CLASSES, rotation=45, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(OVERLAP_CLASSES)
    ax.set_xlabel(out_label)
    ax.set_ylabel(in_label)
    ax.set_title("Gene families by class in each group")
    fig.colorbar(im, ax=ax, label="log10(families + 1)")
    fig.tight_layout()
    _savefig_both(fig, out_dir, "group_class_overlap_heatmap")
    plt.close(fig)


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
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[1] for r in rows])
    ax.invert_yaxis()
    ax.set_xlabel("Gene families")
    ax.set_title("Shared and group-specific gene families")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=3, frameon=False)
    fig.tight_layout()
    _savefig_both(fig, out_dir, "group_class_overlap_shared")
    plt.close(fig)


def plot_group_class_overlap_upset(cells, in_label: str, out_label: str, out_dir: Path) -> None:
    """Spec 2b: intersection sizes over a dot matrix of the 12 group-class
    sets (6 per group
    `absent` is not a set)."""
    inter = overlap_intersections(cells)
    sets = [("IN", c) for c in OVERLAP_CLASSES[:-1]] + [("OUT", c) for c in OVERLAP_CLASSES[:-1]]
    set_names = [f"{in_label if g == 'IN' else out_label}: {c}" for g, c in sets]
    fig = plt.figure(figsize=(max(8, 0.28 * len(inter) + 6), 7))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 4], height_ratios=[2, 1.6], hspace=0.05, wspace=0.02,
                          left=0.22, right=0.98)
    ax_bar = fig.add_subplot(gs[0, 1])
    # Not sharey with ax_set: a shared axis leaks ax_set's tick labels into the dot panel.
    ax_dot = fig.add_subplot(gs[1, 1], sharex=ax_bar)
    ax_set = fig.add_subplot(gs[1, 0])
    x = np.arange(len(inter))
    ax_bar.bar(x, [n for _, _, n in inter], color="#444444")
    ax_bar.set_ylabel("Families")
    ax_bar.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_bar.set_title("Gene families by group-class intersection")
    ys = {s: i for i, s in enumerate(sets)}
    ax_dot.scatter(np.repeat(x, len(sets)), np.tile(range(len(sets)), len(x)), s=10, color="#dddddd")
    for xi, (a, b, _) in enumerate(inter):
        members = [ys[("IN", a)]] if a else []
        members += [ys[("OUT", b)]] if b else []
        ax_dot.plot([xi] * len(members), members, "-o", color="#222222", markersize=4)
    ax_dot.set_yticks([])
    ax_dot.set_ylim(len(sets) - 0.5, -0.5)
    ax_dot.set_xticks([])
    size = [sum(n for (a, b), n in cells.items() if (g == "IN" and a == c) or (g == "OUT" and b == c))
            for g, c in sets]
    ax_set.barh(range(len(sets)), size, color=[BAND_COLORS[c] for _, c in sets])
    ax_set.invert_xaxis()
    ax_set.set_ylim(len(sets) - 0.5, -0.5)
    ax_set.set_yticks(range(len(sets)))
    ax_set.set_yticklabels(set_names, fontsize=7)
    ax_set.set_xlabel("Set size")
    _savefig_both(fig, out_dir, "group_class_overlap_upset")
    plt.close(fig)


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


def build_presence_bool_array(
    matrix: PresenceMatrix, family_order: list[str], strain_order: list[str],
) -> np.ndarray:
    """Materialize a `len(family_order) x len(strain_order)` boolean presence
    array in ONE pass over `matrix.calls` (O(populated cells)), instead of
    the double loop this replaced -- `for fam in family_order: for strain in
    strain_order: arr[i, j] = matrix.is_present(fam, strain)` -- which was
    O(families * strains) Python-level dict lookups + function calls. At real
    study scale (47,743 tier-1 families x 529 strains) that's ~25M calls,
    and `plot_presence_absence_matrix` and `accumulation_curve` each built
    this same array redundantly from scratch.

    `numpy` is not otherwise a `lib/` dependency (only `bin/` scripts import
    it), so this stays in `bin/pangenome_report_render.py` rather than
    becoming a `PresenceMatrix` method, to avoid adding a new dependency to
    `lib/pangenome_matrix.py` for this alone.

    `family_order`/`strain_order` set both the row/column ordering AND which
    rows/columns are included in the result -- a family or strain omitted
    from these lists is simply not represented, same as never calling
    `is_present` for it. A `(family, strain)` pair in `matrix.calls` that
    isn't covered by `family_order`/`strain_order` is skipped. A family or
    strain in `family_order`/`strain_order` with no entry in `matrix.calls`
    at all defaults to False (ABSENT), matching `PresenceMatrix.call`'s own
    default.
    """
    family_idx = {fam: i for i, fam in enumerate(family_order)}
    strain_idx = {strain: j for j, strain in enumerate(strain_order)}
    arr = np.zeros((len(family_order), len(strain_order)), dtype=bool)
    for (fam, strain), state in matrix.calls.items():
        i = family_idx.get(fam)
        if i is None:
            continue
        j = strain_idx.get(strain)
        if j is None:
            continue
        if state in (PRESENT, GENOME_ONLY):
            arr[i, j] = True
    return arr


def plot_presence_absence_matrix(matrix: PresenceMatrix, frequency_table_rows: list[dict], out_dir: Path) -> None:
    """Ported from plot_pangenome_summary.py's plot_presence_absence_matrix
    -- families (rows, frequency-sorted) x strains raster via imshow, kept
    fast at tens-of-thousands-of-families scale by rasterizing rather than
    drawing one element per cell."""
    freq_by_family = {r["family"]: float(r["frequency"]) for r in frequency_table_rows}
    family_order = sorted((f for f in matrix.families if f in freq_by_family), key=lambda f: -freq_by_family[f])
    strain_order = list(matrix.strains)
    arr = build_presence_bool_array(matrix, family_order, strain_order)

    fig, ax = plt.subplots(figsize=(8, 10))
    ax.imshow(arr, aspect="auto", cmap="Greys", interpolation="nearest")
    ax.set_xlabel(f"Strains (n={arr.shape[1]})")
    ax.set_ylabel(f"Gene families (n={arr.shape[0]}), sorted by frequency")
    ax.set_title("Presence/absence matrix")
    ax.set_xticks([])
    ax.set_yticks([])
    # Bug fix: this is a large raster with no other annotation of what
    # black/white mean -- a reader has no way to know which color is
    # "present" without this legend. "Greys" maps the boolean True (present)
    # to black and False (absent) to white, so the patch colors here must
    # mirror that mapping exactly.
    legend_handles = [
        Patch(facecolor="black", edgecolor="black", label="Present"),
        Patch(facecolor="white", edgecolor="black", label="Absent"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", bbox_to_anchor=(1.0, 1.08),
              ncol=2, frameon=True, fontsize="small")
    fig.tight_layout()
    _savefig_both(fig, out_dir, "presence_absence_matrix")
    plt.close(fig)


def accumulation_curve(matrix: PresenceMatrix, n_permutations: int = 20, seed: int = 0) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Ported from plot_pangenome_summary.py's accumulation_curve --
    (pan_mean, pan_std, core_mean, core_std) averaged over `n_permutations`
    random strain orderings."""
    n_strains = len(matrix.strains)
    n_families = len(matrix.families)
    arr = build_presence_bool_array(matrix, matrix.families, matrix.strains)

    rng = random.Random(seed)
    pan_runs = np.zeros((n_permutations, n_strains), dtype=int)
    core_runs = np.zeros((n_permutations, n_strains), dtype=int)
    strain_indices = list(range(n_strains))
    for p in range(n_permutations):
        order = strain_indices[:]
        rng.shuffle(order)
        cum_any = np.zeros(n_families, dtype=bool)
        cum_all = np.ones(n_families, dtype=bool)
        for k, s in enumerate(order):
            cum_any |= arr[:, s]
            cum_all &= arr[:, s]
            pan_runs[p, k] = cum_any.sum()
            core_runs[p, k] = cum_all.sum()
    return (
        pan_runs.mean(axis=0), pan_runs.std(axis=0),
        core_runs.mean(axis=0), core_runs.std(axis=0),
    )


def fit_heaps_law(pan_mean: np.ndarray) -> dict:
    """Fit Heaps' law n = kappa * N^gamma to the pangenome accumulation
    curve via log-log linear regression. gamma < 1 conventionally indicates
    an OPEN pangenome (new families keep appearing as more genomes are
    added); gamma >= 1 indicates a closed one. This is the fitted openness
    statistic Afumigatus's own report never computed -- it only asserted
    openness from the curve's visual shape.

    Below 3 strains, np.polyfit's log-log linear regression degenerates
    (2 points fit a line trivially with a meaningless/undefined R^2, 1 point
    can't fit at all) -- short-circuits to the same not-available sentinel
    fit_core_decay's own except-path returns, rather than attempting the
    fit."""
    if len(pan_mean) < 3:
        return {"kappa": float("nan"), "gamma": float("nan"), "r_squared": float("nan"),
                "is_open": False, "fit_ok": False}
    n_strains = np.arange(1, len(pan_mean) + 1)
    log_n = np.log(n_strains)
    log_pan = np.log(pan_mean)
    gamma, log_kappa = np.polyfit(log_n, log_pan, 1)
    kappa = np.exp(log_kappa)
    predicted = log_kappa + gamma * log_n
    ss_res = np.sum((log_pan - predicted) ** 2)
    ss_tot = np.sum((log_pan - log_pan.mean()) ** 2)
    r_squared = float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
    return {"kappa": float(kappa), "gamma": float(gamma), "r_squared": r_squared,
            "is_open": bool(gamma < 1), "fit_ok": True}


def fit_core_decay(core_mean: np.ndarray) -> dict:
    """Fit an exponential-decay model core(n) = core_inf + amplitude *
    exp(-n / tau) to the core-genome accumulation curve, reporting the
    extrapolated asymptotic core-genome size (core_inf) -- the standard
    companion statistic to Heaps' law openness."""
    from scipy.optimize import curve_fit

    n_strains = np.arange(1, len(core_mean) + 1).astype(float)

    def model(n, core_inf, amplitude, tau):
        return core_inf + amplitude * np.exp(-n / tau)

    try:
        popt, _ = curve_fit(
            model, n_strains, core_mean,
            p0=[core_mean[-1], core_mean[0] - core_mean[-1], max(len(core_mean) / 4, 1)],
            maxfev=10000,
        )
        core_inf, _amplitude, tau = popt
        return {"core_inf": float(core_inf), "tau": float(tau), "fit_ok": True}
    except Exception:
        # curve_fit raises RuntimeError when the fit doesn't converge, but
        # also TypeError when there are fewer data points than free
        # parameters (a real failure mode for a very small strain cohort,
        # e.g. 2 strains) -- both are non-fatal, fall-back-gracefully cases
        # for this report's last pipeline step, not something that should
        # crash REPORT_RENDER.
        return {"core_inf": float(core_mean[-1]), "tau": float("nan"), "fit_ok": False}


def plot_accumulation_curve(pan_mean, pan_std, core_mean, core_std, heaps_fit: dict, out_dir: Path) -> None:
    n = len(pan_mean)
    x = np.arange(1, n + 1)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, pan_mean, color="#2c7fb8", label="Pangenome (union)")
    ax.fill_between(x, pan_mean - pan_std, pan_mean + pan_std, color="#2c7fb8", alpha=0.2)
    ax.plot(x, core_mean, color="#d95f02", label="Core (intersection)")
    ax.fill_between(x, core_mean - core_std, core_mean + core_std, color="#d95f02", alpha=0.2)
    ax.set_xlabel("Number of strains sampled")
    ax.set_ylabel("Number of gene families")
    openness = "open" if heaps_fit["is_open"] else "closed"
    ax.set_title(f"Accumulation curve (Heaps' γ={heaps_fit['gamma']:.2f}, {openness})")
    ax.legend()
    fig.tight_layout()
    _savefig_both(fig, out_dir, "accumulation_curve")
    plt.close(fig)


def plot_classification_counts(classification_counts_dict: dict[str, int], out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    labels = list(classification_counts_dict.keys())
    values = [classification_counts_dict[k] for k in labels]
    ax.bar(labels, values)
    ax.set_ylabel("Number of pairs")
    ax.set_title("Pair classification breakdown")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    _savefig_both(fig, out_dir, "pair_classification_summary")
    plt.close(fig)


def plot_island_size_distribution(size_dist: dict[int, int], out_dir: Path) -> None:
    sizes = sorted(size_dist)
    counts = [size_dist[s] for s in sizes]
    # Numeric x axis so missing sizes show as gaps. Real runs reach sizes of
    # ~70-80, so widen the figure with the largest size, label every size up
    # to 10, then every 5 above it to keep the labels from overlapping.
    max_size = sizes[-1] if sizes else 0
    fig, ax = plt.subplots(figsize=(max(8.0, 0.15 * max_size), 5))
    ax.bar(sizes, counts, width=0.8)
    ax.set_xticks(list(range(min(sizes, default=1), min(max_size, 10) + 1))
                  + list(range(15, max_size + 1, 5)))
    ax.tick_params(axis="x", labelsize=8)
    ax.set_xlabel("Island size (genes)")
    ax.set_ylabel("Number of islands")
    ax.set_title("Accessory island size distribution")
    fig.tight_layout()
    _savefig_both(fig, out_dir, "island_size_distribution")
    plt.close(fig)


def _band_total(row: dict) -> int:
    return sum(_as_int(row.get(band)) for band in GENOME_SEGMENTS)


def read_samplesheet_groups(path: str) -> dict[str, tuple[str, str]]:
    """{Short: (GROUP, Species)} from the pangenome samplesheet CSV."""
    with open(path, newline="") as fh:
        return {r["Short"]: (r.get("GROUP", ""), r.get("Species", ""))
                for r in csv.DictReader(fh)}


def order_genome_blocks(
    per_strain_rows: list[dict],
    samplesheet: dict[str, tuple[str, str]] | None,
    ingroup_label: str,
    outgroup_label: str,
) -> list[tuple[str, list[dict]]]:
    """Genomes for the per-genome class figure, as (block label, rows)
    blocks: ingroup, then outgroup, then every other GROUP value or strain
    missing from the samplesheet (one "other" block); one block per Species
    inside ingroup/outgroup; largest total family count first in a block.
    Without a samplesheet, one unlabelled block sorted by total."""
    def by_total(rows):
        return sorted(rows, key=lambda r: (-_band_total(r), r["Short"]))

    if not samplesheet:
        return [("", by_total(per_strain_rows))]
    grouped: dict[tuple[int, str], list[dict]] = {}
    for row in per_strain_rows:
        group, species = samplesheet.get(row["Short"], ("", ""))
        if group == ingroup_label:
            key = (0, f"{species} ({group})")
        elif group == outgroup_label:
            key = (1, f"{species} ({group})")
        else:
            key = (2, "other")
        grouped.setdefault(key, []).append(row)
    return [(label, by_total(rows)) for (_, label), rows in sorted(grouped.items())]


PER_GENOME_LABEL_MAX = 150
"""Above this many genomes the per-genome figure drops its genome-name
labels: they no longer fit at a readable size. Values stay in
per_strain_summary.tsv."""


def per_genome_figure_height(n_genomes: int) -> tuple[float, bool]:
    """(figure height in inches, draw genome labels?) for the per-genome
    class figure. Labelled rows get 0.15in each; above
    PER_GENOME_LABEL_MAX the rows shrink to 0.04in, capped at 30in."""
    if n_genomes <= PER_GENOME_LABEL_MAX:
        return 1.5 + 0.15 * n_genomes, True
    return min(30.0, 1.5 + 0.04 * n_genomes), False


def plot_per_genome_class_composition(blocks: list[tuple[str, list[dict]]], out_dir: Path) -> None:
    """One horizontal stacked bar per genome: family counts per frequency
    band (after PPanGGOLiN's Fig 3, without the tree). Top to bottom in
    `blocks` order; a divider and label mark each block. Segments with no
    families in any genome are not drawn (older per_strain_summary.tsv files
    have no nonrep_only/outgroup_only columns). Non-representative genomes
    (`is_representative == "N"`) get a dagger after the name, or a tick at
    the left edge when names are not drawn."""
    rows = [r for _, block in blocks for r in block]
    height, show_labels = per_genome_figure_height(len(rows))
    fig, ax = plt.subplots(figsize=(10, height))
    y = np.arange(len(rows))
    left = np.zeros(len(rows))
    drawn = []
    for band in GENOME_SEGMENTS:
        values = np.array([_as_int(r.get(band)) for r in rows], dtype=float)
        if not values.any():
            continue
        drawn.append(band)
        # Unlabelled (many-genome) rows touch: thin gaps between hundreds of
        # bars alias into false stripes in the PNG.
        ax.barh(y, values, left=left, height=0.85 if show_labels else 1.0,
                color=BAND_COLORS[band], label=band, linewidth=0,
                antialiased=show_labels)
        left += values
    nonrep = [r.get("is_representative", "Y") == "N" for r in rows]
    if show_labels:
        ax.set_yticks(y)
        ax.set_yticklabels([r["Short"] + (" †" if n else "") for r, n in zip(rows, nonrep)], fontsize=6)
    else:
        ax.set_yticks([])
        ax.set_ylabel(f"{len(rows)} genomes (per-genome values in per_strain_summary.tsv)")
        ys = [i for i, n in enumerate(nonrep) if n]
        if ys:
            ax.scatter([0] * len(ys), ys, marker="_", s=40, linewidths=1.2, color="black",
                       transform=ax.get_yaxis_transform(), clip_on=False, zorder=3)
    start = 0
    for i, (label, block) in enumerate(blocks):
        if i:
            ax.axhline(start - 0.5, color="black", linewidth=0.8)
        if label:
            ax.text(1.005, start - 0.4, label, transform=ax.get_yaxis_transform(),
                    va="top", ha="left", fontsize=8)
        start += len(block)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xlabel("Gene families")
    ax.set_title("Gene families per genome by frequency class")
    if any(nonrep):
        # Marker-less legend entry: explains the dagger / tick without extra
        # figure space (short figures have none to spare).
        mark = "†" if show_labels else "left tick"
        ax.plot([], [], " ", label=f"{mark} = non-representative genome")
    handles, labels = ax.get_legend_handles_labels()
    order = sorted(range(len(labels)), key=lambda i: labels[i] not in drawn)  # bands first, note last
    ax.legend([handles[i] for i in order], [labels[i] for i in order],
              loc="upper center", bbox_to_anchor=(0.5, -0.04 * 20 / height),
              ncol=max(1, min(4, len(labels))), frameon=False)
    fig.tight_layout()
    _savefig_both(fig, out_dir, "per_genome_class_composition")
    plt.close(fig)


_ZERO_Q_FALLBACK_NEG_LOG_Q = 10.0
"""Fixed fallback for cap_zero_q_sentinel's all-zero-q edge case: chosen as
a value clearly visible on a -log10(q) axis (an FDR q of 1e-10 would be an
extraordinarily strong result) without depending on any finite value in the
data (there is none to scale off when every plotted domain has fdr_q == 0)."""


def cap_zero_q_sentinel(
    neg_log_q: list[float], is_zero_q: list[bool], multiple: float = 1.3,
) -> list[float]:
    """Replaces each zero-q placeholder value in `neg_log_q` (positions
    flagged by `is_zero_q`) with a value scaled off the largest FINITE
    (non-zero-q) value actually being plotted, capped at `multiple` times
    that max -- instead of a fixed sentinel (the previous 300.0) that
    dwarfs every real, meaningful bar on the same linear axis regardless of
    what scale the rest of the chart is at.

    A q=0 domain is the most statistically significant result on the chart
    by construction (q=0 beats any q>0), so it should still read as the
    tallest bar -- just not so tall the others become invisible by
    comparison. `multiple` (default 1.3x the max finite value) keeps that
    "still tallest" property while staying on a comparable visual scale.

    Falls back to `_ZERO_Q_FALLBACK_NEG_LOG_Q` when every plotted value is
    a zero-q placeholder (`is_zero_q` all True) -- there is no finite value
    to scale off in that edge case.
    """
    finite_values = [v for v, is_zero in zip(neg_log_q, is_zero_q) if not is_zero]
    if finite_values:
        capped_value = max(finite_values) * multiple
    else:
        capped_value = _ZERO_Q_FALLBACK_NEG_LOG_Q
    return [capped_value if is_zero else v for v, is_zero in zip(neg_log_q, is_zero_q)]


def plot_domain_enrichment(top_domains: list[dict], out_dir: Path, top_n: int = 20) -> None:
    """Horizontal bar chart of the top-N (by fisher_p; `top_domains` arrives
    already sorted that way from main()) significantly enriched
    Pfam domains among island-member families -- the figure
    (`island_domain_enrichment.png`) Afumigatus's own REPORT.md embedded by
    hand with NO checked-in generator; this makes it reproducible."""
    domains = top_domains[:top_n]
    if not domains:
        return
    fig, ax = plt.subplots(figsize=(8, max(4, 0.3 * len(domains))))
    labels = [d["domain"] for d in reversed(domains)]
    fdr_qs = [float(d["fdr_q"]) for d in reversed(domains)]
    is_zero_q = [q == 0 for q in fdr_qs]
    raw_neg_log_q = [-np.log10(q) if q > 0 else 0.0 for q in fdr_qs]
    neg_log_q = cap_zero_q_sentinel(raw_neg_log_q, is_zero_q)
    ax.barh(labels, neg_log_q, color="#2c7fb8")
    ax.set_xlabel("-log10(FDR q)")
    ax.set_title(f"Top {len(domains)} enriched Pfam domains in accessory islands")
    fig.tight_layout()
    _savefig_both(fig, out_dir, "island_domain_enrichment")
    plt.close(fig)


def _as_int(value, default: int = 0) -> int:
    """Parse a TSV cell to int, tolerating '' and the '-' missing-value sentinel."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _fmt(value, spec: str) -> str:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "-"
    return "-" if v != v else format(v, spec)


def _neighborhood_section(rows: list[dict]) -> list[str]:
    """View B1 table (issue #182): per-module genomic clustering in each strain's
    own assembly, against a within-strain permutation null. The scale and null
    pool are printed with the table (spec: never left implicit)."""
    lines = ["## Trans-module genomic clustering", ""]
    if not rows:
        return lines + ["No trans modules to score.", ""]
    r0 = rows[0]
    lines += [
        f"Scale: same contig, at least {r0['min_gene_gap']} genes apart, gene starts within "
        f"max {_fmt(r0['max_kb'], 'g')} kb. Null: {r0['null_pool']} genes of the same strain, "
        f"{r0['n_perm']} permutations, seed {r0['seed']}. Pairs on different contigs are "
        "excluded, not counted as distant. With millions of pairs a tiny effect still "
        "reaches the smallest possible p, so read the effect ratio first.",
        "",
        "| Module | Families | Strains | Obs. frac (same contig) | Null | Effect | p "
        "| Effect (all pairs) | p (all pairs) | Cross-contig excluded |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        xc = _fmt(float(r["cross_contig_frac"]) * 100, ".1f") if r.get("cross_contig_frac") not in (None, "", "nan") else "-"
        lines.append(
            f"| {r['module_id']} | {r['module_size']} | {r['strains_scored']} "
            f"| {_fmt(r['obs_frac'], '.3f')} | {_fmt(r['null_mean_frac'], '.3f')} "
            f"| {_fmt(r['effect_ratio'], '.2f')} | {_fmt(r['p_empirical'], '.3f')} "
            f"| {_fmt(r['effect_ratio_total'], '.2f')} | {_fmt(r['p_empirical_total'], '.3f')} "
            f"| {xc}% |")
    return lines + [""]


def _linked_pfam_domains(pfam_domains: str, pfam_urls: dict[str, str] | None) -> str:
    """Comma-joined Pfam domain names -> markdown links to each domain's Pfam
    page. A name with no URL ('-', or absent from `pfam_urls`) stays plain."""
    if not pfam_domains or pfam_domains == "-":
        return "-"
    cells = []
    for name in pfam_domains.split(","):
        url = (pfam_urls or {}).get(name, "-")
        cells.append(f"[{name}]({url})" if url and url != "-" else name)
    return ", ".join(cells)


def render_report_markdown(
    counts: dict[str, dict[str, int]],
    size_dist: dict[int, int],
    classification_counts_dict: dict[str, int],
    top_domains: list[dict],
    n_islands: int,
    heaps_fit: dict | None,
    core_decay: dict | None,
    strain_family_counts: list[int],
    marker_rows: list[dict] | None = None,
    islands_with_domains_rows: list[dict] | None = None,
    per_strain_rows: list[dict] | None = None,
    top_islands_min_strains: int = 2,
    diagnostics_banner: str | None = None,
    islands_available: bool = True,
    neighborhood_rows: list[dict] | None = None,
    pfam_urls: dict[str, str] | None = None,
    per_genome_figure: bool = False,
    outgroup_fallback: str = "",
    overlap_section: str = "",
) -> str:
    lines: list[str] = []
    if diagnostics_banner:
        # Issue #134: pipeline diagnostics land at the TOP of report.md,
        # before any results -- a reader must see them before the first
        # figure, not scroll past them or find them only in stderr.
        lines += [diagnostics_banner.rstrip("\n"), ""]
    lines += ["# Pangenome Island + Pfam Enrichment Report", ""]
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
            lines += ["", "Not counted in this group's classes:"] + extra
        lines.append("")
    if counts:
        lines += ["![Composition](figures/group_composition.png)", ""]
    lines += ["![Frequency distribution](figures/frequency_distribution.png)", ""]

    if heaps_fit is not None:
        lines += ["## Pangenome openness", ""]
        openness = "open" if heaps_fit["is_open"] else "closed"
        lines += [f"Heaps' law fit: κ={heaps_fit['kappa']:.1f}, γ={heaps_fit['gamma']:.3f} "
                  f"(R²={heaps_fit['r_squared']:.3f}) -- pangenome is **{openness}** "
                  f"(γ {'<' if heaps_fit['is_open'] else '>='} 1).", ""]
        if core_decay is not None and core_decay["fit_ok"]:
            lines += [f"Extrapolated asymptotic core-genome size: {core_decay['core_inf']:.0f} families.", ""]
        lines += ["![Accumulation curve](figures/accumulation_curve.png)", ""]
        lines += ["![Presence/absence matrix](figures/presence_absence_matrix.png)", ""]

    if strain_family_counts:
        lines += ["## Per-strain summary", ""]
        lines += [f"Families per strain: min={min(strain_family_counts)}, "
                  f"median={sorted(strain_family_counts)[len(strain_family_counts)//2]}, "
                  f"max={max(strain_family_counts)} (n={len(strain_family_counts)} strains). "
                  "See per_strain_summary.tsv for outliers.", ""]
        flagged = [r["Short"] for r in (per_strain_rows or []) if r.get("is_outlier") == "Y"]
        if flagged:
            lines += [f"**Outlier strains (singleton-count modified z-score beyond threshold):** "
                      f"{', '.join(flagged)}", ""]
    if per_genome_figure:
        lines += ["![Gene families per genome by frequency class]"
                  "(figures/per_genome_class_composition.png)", ""]
        if outgroup_fallback == "too_few":
            lines += ["*The outgroup has too few representative genomes to bin on its own; "
                      "outgroup genomes use ingroup classes.*", ""]
        elif outgroup_fallback == "no_out_columns":
            lines += ["*This run's frequency table predates per-group bins (no outgroup "
                      "columns); outgroup genomes use ingroup classes.*", ""]

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
                      "Counts: the run's `group_class_overlap.tsv` table.", ""]
        elif overlap_section == "no_out_columns":
            lines += ["*Not shown: this run's frequency table predates per-group bins.*", ""]
        else:
            lines += ["*Not shown: the outgroup is not binned (too few representative genomes).*", ""]

    if islands_available:
        lines += ["## Accessory islands", ""]
        lines += [f"{n_islands} statistically significant accessory islands found "
                  "(built from adjacency of non-core genes, gated by containing "
                  "at least one FDR-significant physically-linked pair).", ""]
        if size_dist:
            lines += ["![Island sizes](figures/island_size_distribution.png)", ""]

        located = [r for r in (islands_with_domains_rows or []) if r.get("locus_id", "-") != "-"]
        top_islands = [r for r in located
                       if _as_int(r.get("n_strains")) >= top_islands_min_strains]
        n_excluded = len(located) - len(top_islands)
        if top_islands:
            top_islands.sort(key=lambda r: -int(r.get("island_size", 0)))
            lines += ["", "**Top islands (by size):**", ""]
            if n_excluded:
                if top_islands_min_strains == 2:
                    note = (f"*{n_excluded} single-strain islands excluded "
                            "(present in one strain only -- strain-private content, "
                            "which dominates the size ranking).*")
                else:
                    note = (f"*{n_excluded} islands present in fewer than "
                            f"{top_islands_min_strains} strains excluded.*")
                lines += [note, ""]
            lines += ["| Locus (strain:contig:start-end) | Families (#) | Span (kb) | "
                      "Strains (#) | Pfam domains |",
                      "|---|---|---|---|---|"]
            for row in top_islands[:20]:
                start, end = _as_int(row.get("locus_start"), -1), _as_int(row.get("locus_end"), -1)
                span = f"{(end - start) / 1000:.1f}" if start >= 0 and end >= 0 else "-"
                lines.append(f"| {row.get('locus_id', '-')} | {row.get('island_size', '-')} | "
                             f"{span} | {row.get('n_strains', '-')} | "
                             f"{_linked_pfam_domains(row.get('pfam_domains', '-'), pfam_urls)} |")
            lines.append("")

    if marker_rows:
        lines += ["## Marker co-occurrence", ""]
        lines += ["| Marker | Islands with marker | Islands total | % with marker |",
                  "|---|---|---|---|"]
        for row in marker_rows:
            lines.append(
                f"| {row['marker_name']} | {row['n_islands_with_marker']} | "
                f"{row['n_islands_total']} | {float(row['pct_islands_with_marker']):.1f}% |"
            )
        lines.append("")

    lines += ["## Pair classification breakdown", ""]
    for classification, count in sorted(classification_counts_dict.items(), key=lambda kv: -kv[1]):
        lines.append(f"- **{classification}**: {count}")
    if classification_counts_dict:
        lines += ["", "![Classification breakdown](figures/pair_classification_summary.png)", ""]

    if neighborhood_rows is not None:
        lines += _neighborhood_section(neighborhood_rows)

    if islands_available:
        lines += ["## Pfam domain enrichment", ""]
        if not top_domains:
            lines += ["No significantly enriched Pfam domains found.", ""]
        else:
            lines += ["![Top enriched domains](figures/island_domain_enrichment.png)", ""]
            has_go = any(row.get("go_terms") for row in top_domains)
            if has_go:
                lines += ["| Domain | Fisher p | FDR q | GO terms |", "|---|---|---|---|"]
            else:
                lines += ["| Domain | Fisher p | FDR q |", "|---|---|---|"]
            for row in top_domains:
                pfam_url = row.get("pfam_url")
                domain_cell = f"[{row['domain']}]({pfam_url})" if pfam_url and pfam_url != "-" else row["domain"]
                cells = f"| {domain_cell} | {float(row['fisher_p']):.2e} | {float(row['fdr_q']):.2e} |"
                if has_go:
                    cells += f" {row.get('go_terms', '-')} |"
                lines.append(cells)
            lines.append("")

    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--presence_matrix", required=True)
    ap.add_argument(
        "--islands_with_domains", default=None,
        help="Optional (issue #135 core-report ungating): islands_with_domains.tsv "
        "from REPORT_TABLES. Omitted (or an empty file) on a run with no "
        "--pangenome_island_pfam_hmm -- the 'Accessory islands' section is then "
        "skipped rather than reporting '0 islands found'.",
    )
    ap.add_argument("--island_size_distribution", default=None)
    ap.add_argument("--classification_counts", required=True)
    ap.add_argument(
        "--island_pfam_enrichment", default=None,
        help="Optional, same as --islands_with_domains -- omitted skips the "
        "'Pfam domain enrichment' section.",
    )
    ap.add_argument("--marker_summary", default=None)
    ap.add_argument("--per_strain_summary", required=True)
    ap.add_argument(
        "--samplesheet", default=None,
        help="Optional pangenome samplesheet CSV (GROUP,Species,...,Short). Orders the "
        "per-genome class figure by GROUP, then Species; without it, by total only.")
    ap.add_argument("--ingroup_label", default="IN")
    ap.add_argument("--outgroup_label", default="OUT")
    ap.add_argument("--top_islands_min_strains", type=int, default=2,
                    help="Minimum carrying strains for an island to appear in the "
                         "report's 'Top islands' table (default 2: exclude "
                         "strain-private islands). 1 = no filter.")
    ap.add_argument("--fdr_threshold", type=float, default=0.05)
    ap.add_argument("--n_permutations", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--diagnostics_banner", default=None,
        help="Optional pangenome_diagnostics.py diagnostics_banner.md file "
        "(issue #134) -- prepended to report.md, before any results.",
    )
    ap.add_argument("--group_class_overlap", default=None,
                    help="report_tables/group_class_overlap.tsv (#212 PR 2); omitted -> no "
                    "'Ingroup vs outgroup content' section, empty -> section with the skip reason.")
    ap.add_argument("--module_neighborhood", default=None,
                    help="module_neighborhood.tsv (View B1, issue #182); optional")
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.frequency_table, newline="") as fh:
        frequency_table_rows = list(csv.DictReader(fh, delimiter="\t"))
    counts = group_class_counts(frequency_table_rows)

    size_dist: dict[int, int] = {}
    if args.island_size_distribution:
        with open(args.island_size_distribution, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                size_dist[int(row["island_size"])] = int(row["count"])

    classification_counts_dict: dict[str, int] = {}
    with open(args.classification_counts, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            classification_counts_dict[row["classification"]] = int(row["count"])

    # islands_available drives render_report_markdown's "Accessory islands"/
    # "Pfam domain enrichment" sections -- see --islands_with_domains's help
    # text. Only --islands_with_domains (not --island_size_distribution,
    # which is empty whenever size_dist is empty regardless of reason) is the
    # signal, matching how REPORT_TABLES's islands_available flag is decided.
    islands_available = args.islands_with_domains is not None
    islands_with_domains_rows: list[dict] = []
    n_islands = 0
    if islands_available:
        with open(args.islands_with_domains, newline="") as fh:
            islands_with_domains_rows = list(csv.DictReader(fh, delimiter="\t"))
        n_islands = len(islands_with_domains_rows)

    top_domains: list[dict] = []
    # Every tested domain (not only the FDR-significant ones) has a row in
    # the enrichment table, so this covers all names in the Top islands table.
    pfam_urls: dict[str, str] = {}
    if args.island_pfam_enrichment:
        with open(args.island_pfam_enrichment, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                pfam_urls[row["domain"]] = row.get("pfam_url", "-")
                if float(row["fdr_q"]) < args.fdr_threshold:
                    top_domains.append(row)
        top_domains.sort(key=lambda r: float(r["fisher_p"]))

    strain_family_counts: list[int] = []
    per_strain_rows: list[dict] = []
    with open(args.per_strain_summary, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            strain_family_counts.append(int(row["n_families"]))
            per_strain_rows.append(row)

    marker_rows: list[dict] = []
    if args.marker_summary:
        with open(args.marker_summary, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                marker_rows.append(row)

    plot_frequency_distribution(frequency_table_rows, out_dir)
    plot_group_composition(counts, out_dir)
    if size_dist:
        plot_island_size_distribution(size_dist, out_dir)
    if classification_counts_dict:
        plot_classification_counts(classification_counts_dict, out_dir)
    plot_domain_enrichment(top_domains, out_dir)
    # Older per_strain_summary.tsv files carry no per-band columns.
    per_genome_figure = bool(per_strain_rows) and all(b in per_strain_rows[0] for b in BAND_ORDER)
    if per_genome_figure:
        samplesheet = read_samplesheet_groups(args.samplesheet) if args.samplesheet else None
        plot_per_genome_class_composition(
            order_genome_blocks(per_strain_rows, samplesheet,
                                args.ingroup_label, args.outgroup_label),
            out_dir)
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

    matrix = PresenceMatrix.from_tsv(args.presence_matrix)
    plot_presence_absence_matrix(matrix, frequency_table_rows, out_dir)
    pan_mean, pan_std, core_mean, core_std = accumulation_curve(
        matrix, n_permutations=args.n_permutations, seed=args.seed,
    )
    heaps_fit = fit_heaps_law(pan_mean)
    core_decay = fit_core_decay(core_mean)
    plot_accumulation_curve(pan_mean, pan_std, core_mean, core_std, heaps_fit, out_dir)

    with open(out_dir / "pangenome_openness.tsv", "w") as out:
        out.write("kappa\tgamma\tr_squared\tis_open\tcore_inf\ttau\n")
        out.write(f"{heaps_fit['kappa']:.4f}\t{heaps_fit['gamma']:.4f}\t{heaps_fit['r_squared']:.4f}\t"
                  f"{heaps_fit['is_open']}\t{core_decay['core_inf']:.4f}\t{core_decay['tau']:.4f}\n")

    diagnostics_banner = (
        Path(args.diagnostics_banner).read_text() if args.diagnostics_banner else None
    )
    neighborhood_rows = None
    if args.module_neighborhood:
        with open(args.module_neighborhood, newline="") as fh:
            neighborhood_rows = list(csv.DictReader(fh, delimiter="\t"))
    markdown = render_report_markdown(
        counts, size_dist, classification_counts_dict, top_domains, n_islands,
        heaps_fit, core_decay, strain_family_counts, marker_rows,
        islands_with_domains_rows=islands_with_domains_rows,
        per_strain_rows=per_strain_rows,
        top_islands_min_strains=args.top_islands_min_strains,
        diagnostics_banner=diagnostics_banner,
        islands_available=islands_available,
        neighborhood_rows=neighborhood_rows,
        pfam_urls=pfam_urls,
        per_genome_figure=per_genome_figure,
        outgroup_fallback=outgroup_fallback_reason(frequency_table_rows, per_strain_rows,
                                                   args.outgroup_label),
        overlap_section=overlap_section,
    )
    (out_dir / "report.md").write_text(markdown)

    print(f"pangenome_report_render: wrote {out_dir}/report.md, pangenome_openness.tsv, "
          f"and figures/, figures_pdf/", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
