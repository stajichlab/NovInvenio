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

Usage:
  pangenome_report_render.py --frequency_table frequency_table.tsv \\
      --presence_matrix presence_matrix.tsv \\
      --islands_with_domains islands_with_domains.tsv \\
      --island_size_distribution island_size_distribution.tsv \\
      --classification_counts classification_counts.tsv \\
      --island_pfam_enrichment island_pfam_enrichment.tsv \\
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
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import PresenceMatrix  # noqa: E402

BAND_ORDER = ["core", "soft_core", "shell", "cloud", "singleton"]
BAND_COLORS = {
    "core": "#2c7fb8", "soft_core": "#7fcdbb", "shell": "#fed976",
    "cloud": "#fd8d3c", "singleton": "#bdbdbd",
}


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


def plot_frequency_bins(counts: dict[str, int], out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    labels = [k for k in BAND_ORDER if counts.get(k, 0) > 0]
    values = [counts[k] for k in labels]
    ax.pie(values, labels=labels, autopct="%1.1f%%", colors=[BAND_COLORS[b] for b in labels])
    ax.set_title("Pangenome composition")
    fig.tight_layout()
    _savefig_both(fig, out_dir, "core_shell_cloud_pie")
    plt.close(fig)


def plot_presence_absence_matrix(matrix: PresenceMatrix, frequency_table_rows: list[dict], out_dir: Path) -> None:
    """Ported from plot_pangenome_summary.py's plot_presence_absence_matrix
    -- families (rows, frequency-sorted) x strains raster via imshow, kept
    fast at tens-of-thousands-of-families scale by rasterizing rather than
    drawing one element per cell."""
    freq_by_family = {r["family"]: float(r["frequency"]) for r in frequency_table_rows}
    family_order = sorted((f for f in matrix.families if f in freq_by_family), key=lambda f: -freq_by_family[f])
    strain_order = list(matrix.strains)
    arr = np.zeros((len(family_order), len(strain_order)), dtype=bool)
    for i, fam in enumerate(family_order):
        for j, strain in enumerate(strain_order):
            arr[i, j] = matrix.is_present(fam, strain)

    fig, ax = plt.subplots(figsize=(8, 10))
    ax.imshow(arr, aspect="auto", cmap="Greys", interpolation="nearest")
    ax.set_xlabel(f"Strains (n={arr.shape[1]})")
    ax.set_ylabel(f"Gene families (n={arr.shape[0]}), sorted by frequency")
    ax.set_title("Presence/absence matrix")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    _savefig_both(fig, out_dir, "presence_absence_matrix")
    plt.close(fig)


def accumulation_curve(matrix: PresenceMatrix, n_permutations: int = 20, seed: int = 0) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Ported from plot_pangenome_summary.py's accumulation_curve --
    (pan_mean, pan_std, core_mean, core_std) averaged over `n_permutations`
    random strain orderings."""
    n_strains = len(matrix.strains)
    n_families = len(matrix.families)
    arr = np.zeros((n_families, n_strains), dtype=bool)
    for i, fam in enumerate(matrix.families):
        for j, strain in enumerate(matrix.strains):
            arr[i, j] = matrix.is_present(fam, strain)

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
    fig.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    _savefig_both(fig, out_dir, "pair_classification_summary")
    plt.close(fig)


def plot_island_size_distribution(size_dist: dict[int, int], out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    sizes = sorted(size_dist)
    counts = [size_dist[s] for s in sizes]
    ax.bar([str(s) for s in sizes], counts)
    ax.set_xlabel("Island size (genes)")
    ax.set_ylabel("Number of islands")
    ax.set_title("Accessory island size distribution")
    fig.tight_layout()
    _savefig_both(fig, out_dir, "island_size_distribution")
    plt.close(fig)


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
    neg_log_q = [-np.log10(float(d["fdr_q"])) if float(d["fdr_q"]) > 0 else 300.0 for d in reversed(domains)]
    ax.barh(labels, neg_log_q, color="#2c7fb8")
    ax.set_xlabel("-log10(FDR q)")
    ax.set_title(f"Top {len(domains)} enriched Pfam domains in accessory islands")
    fig.tight_layout()
    _savefig_both(fig, out_dir, "island_domain_enrichment")
    plt.close(fig)


def render_report_markdown(
    counts: dict[str, int],
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
) -> str:
    total_families = sum(counts.values())
    lines = ["# Pangenome Island + Pfam Enrichment Report", ""]
    lines += ["## Pangenome composition", ""]
    lines += [f"Total families: {total_families}", ""]
    for label in BAND_ORDER:
        if counts.get(label, 0):
            pct = 100 * counts[label] / total_families if total_families else 0
            lines.append(f"- **{label}**: {counts[label]} ({pct:.1f}%)")
    lines += ["", "![Composition](figures/core_shell_cloud_pie.png)", ""]
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

    lines += ["## Accessory islands", ""]
    lines += [f"{n_islands} statistically significant accessory islands found "
              "(built from adjacency of non-core genes, gated by containing "
              "at least one FDR-significant physically-linked pair).", ""]
    if size_dist:
        lines += ["![Island sizes](figures/island_size_distribution.png)", ""]

    top_islands = [r for r in (islands_with_domains_rows or []) if r.get("locus_id", "-") != "-"]
    if top_islands:
        top_islands.sort(key=lambda r: -int(r.get("island_size", 0)))
        lines += ["", "**Top islands (by size):**", "",
                  "| Locus | Size | Strains | Pfam domains |", "|---|---|---|---|"]
        for row in top_islands[:20]:
            lines.append(f"| {row.get('locus_id', '-')} | {row.get('island_size', '-')} | "
                          f"{row.get('n_strains', '-')} | {row.get('pfam_domains', '-')} |")
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
    ap.add_argument("--islands_with_domains", required=True)
    ap.add_argument("--island_size_distribution", required=True)
    ap.add_argument("--classification_counts", required=True)
    ap.add_argument("--island_pfam_enrichment", required=True)
    ap.add_argument("--marker_summary", required=True)
    ap.add_argument("--per_strain_summary", required=True)
    ap.add_argument("--fdr_threshold", type=float, default=0.05)
    ap.add_argument("--n_permutations", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.frequency_table, newline="") as fh:
        frequency_table_rows = list(csv.DictReader(fh, delimiter="\t"))
    counts: dict[str, int] = {}
    for row in frequency_table_rows:
        counts[row["bin"]] = counts.get(row["bin"], 0) + 1

    size_dist: dict[int, int] = {}
    with open(args.island_size_distribution, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            size_dist[int(row["island_size"])] = int(row["count"])

    classification_counts_dict: dict[str, int] = {}
    with open(args.classification_counts, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            classification_counts_dict[row["classification"]] = int(row["count"])

    with open(args.islands_with_domains, newline="") as fh:
        islands_with_domains_rows = list(csv.DictReader(fh, delimiter="\t"))
    n_islands = len(islands_with_domains_rows)

    top_domains: list[dict] = []
    with open(args.island_pfam_enrichment, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
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
    with open(args.marker_summary, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            marker_rows.append(row)

    plot_frequency_distribution(frequency_table_rows, out_dir)
    plot_frequency_bins(counts, out_dir)
    if size_dist:
        plot_island_size_distribution(size_dist, out_dir)
    if classification_counts_dict:
        plot_classification_counts(classification_counts_dict, out_dir)
    plot_domain_enrichment(top_domains, out_dir)

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

    markdown = render_report_markdown(
        counts, size_dist, classification_counts_dict, top_domains, n_islands,
        heaps_fit, core_decay, strain_family_counts, marker_rows,
        islands_with_domains_rows=islands_with_domains_rows,
        per_strain_rows=per_strain_rows,
    )
    (out_dir / "report.md").write_text(markdown)

    print(f"pangenome_report_render: wrote {out_dir}/report.md, pangenome_openness.tsv, "
          f"and figures/, figures_pdf/", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
