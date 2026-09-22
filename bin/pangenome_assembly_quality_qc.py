#!/usr/bin/env python3
"""Assembly-quality vs pangenome-content QC (issue #130).

A real 529-strain Coccidioides pangenome study found `rho(accessory ~ N50)
= -0.53` and `rho(accessory ~ n_contigs) = +0.53`: more-fragmented
assemblies show MORE accessory (variable) gene families, not because they
carry more real biology, but because a gene broken across a contig
boundary is predicted as two partial proteins, neither of which clusters
with the intact family -- each becomes a spurious "strain-private" family.
Confirmed mechanism: proteins in strain-private (singleton) families sit
at a contig terminus far more often than proteins genome-wide.

This makes that confound VISIBLE on every pangenome run, automatically,
rather than requiring a one-off manual analysis to rediscover it. Per
ingroup strain, it reports assembly-quality metrics (n_contigs, N50,
total_length -- computed directly from that strain's own DNA FASTA, so
this runs regardless of whether `--pangenome_dereplicate` is enabled)
alongside that strain's accessory-family-present count, core-family-missing
count, and the fraction of that strain's private-family proteins sitting
within `--terminus_window_bp` of a contig end (vs. the same fraction over
all of that strain's proteins). It also reports the genome-wide Spearman
rho of accessory-present/core-missing against n_contigs/N50 -- both raw and
partial (controlling for total_length) -- and warns when |rho| for
accessory-vs-assembly-quality exceeds `--rho_warn_threshold`.

This is a DIAGNOSTIC, not a filter: it never excludes a strain or corrects
a presence call. See DESIGN.md / issue #130 for the full reasoning; the
`--terminus_window_bp` default (1000 bp) is this script's own choice of
"near a contig end", not a reproduction of a specific number from the
originating analysis (which did not check in an exact window value).

Usage:
  pangenome_assembly_quality_qc.py --config config.csv --data_dir data_dir \\
      --matrix presence_matrix.rescued.tsv --frequency_table frequency_table.tsv \\
      --gene_positions gene_positions.tsv.zst --cluster_tsv cluster.tsv \\
      --ingroup_label IN --out_dir .
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import open_maybe_compressed  # noqa: E402
from pangenome_matrix import PresenceMatrix, read_cluster_tsv  # noqa: E402
from config_parser import parse_config  # noqa: E402


def parse_contig_lengths(fasta_path: str | Path) -> dict[str, int]:
    """{contig_name: length} from a DNA FASTA (header first token, before
    any whitespace, matching this repo's standard FASTA-header-as-ID
    convention)."""
    lengths: dict[str, int] = {}
    name = None
    current = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = current
                name = line[1:].split()[0].strip()
                current = 0
            else:
                current += len(line.strip())
        if name is not None:
            lengths[name] = current
    return lengths


def assembly_stats_from_lengths(lengths: dict[str, int]) -> dict:
    """{n_contigs, n50, total_length} from a {contig: length} map."""
    values = sorted(lengths.values(), reverse=True)
    total = sum(values)
    half = total / 2
    running = 0
    n50 = 0
    for v in values:
        running += v
        if running >= half:
            n50 = v
            break
    return {"n_contigs": len(values), "n50": n50, "total_length": total}


def is_at_terminus(start: int, end: int, contig_length: int, window_bp: int) -> bool:
    """True if the (1-based, inclusive) span [start, end] lies within
    `window_bp` of either end of its contig -- the structural signature of
    a gene model truncated by a contig/scaffold boundary."""
    if contig_length <= 0:
        return False
    return start <= window_bp or (contig_length - end) <= window_bp


def find_private_family_strain(matrix: PresenceMatrix, singleton_families: set[str]) -> dict[str, str]:
    """{family: strain} for every singleton (strain_count == 1) family --
    the one strain actually carrying it. A family with zero or more than
    one carrying strain (shouldn't happen for a true singleton bin, but a
    stale/edited frequency_table could disagree with the matrix) is
    silently skipped rather than guessed at."""
    result: dict[str, str] = {}
    for fam in singleton_families:
        carriers = [s for s in matrix.strains if matrix.is_present(fam, s)]
        if len(carriers) == 1:
            result[fam] = carriers[0]
    return result


def _rank(values: list[float]) -> list[float]:
    """Average-rank transform (1-based), ties sharing the mean rank of
    their block -- the standard tie-handling for a Spearman correlation."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    n = len(values)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg_rank = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = avg_rank
        i = j + 1
    return ranks


def pearson(x: list[float], y: list[float]) -> float:
    n = len(x)
    if n < 2:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sx = sum((xi - mx) ** 2 for xi in x)
    sy = sum((yi - my) ** 2 for yi in y)
    if sx == 0 or sy == 0:
        return float("nan")
    sxy = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    return sxy / math.sqrt(sx * sy)


def spearman(x: list[float], y: list[float]) -> float:
    return pearson(_rank(x), _rank(y))


def partial_spearman(x: list[float], y: list[float], z: list[float]) -> float:
    """Partial Spearman correlation of x and y, controlling for z -- Pearson
    on ranks (equivalent to Spearman), with the standard first-order
    partial-correlation formula. NaN if either controlling correlation is
    +-1 (degenerate) or any input correlation is NaN (e.g. < 2 data points,
    or one variable is constant)."""
    rxy, rxz, ryz = spearman(x, y), spearman(x, z), spearman(y, z)
    if any(math.isnan(v) for v in (rxy, rxz, ryz)):
        return float("nan")
    denom = math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    if denom == 0:
        return float("nan")
    return (rxy - rxz * ryz) / denom


def compute_per_strain_content(
    matrix: PresenceMatrix, freq_bin: dict[str, str], strains: list[str],
) -> dict[str, dict]:
    """{strain: {n_families_present, accessory_present, core_missing}}.
    `accessory_present` counts families present in that strain whose
    frequency_table bin is anything but "core" (soft_core/shell/cloud/
    singleton) -- a family with no frequency_table entry at all (should not
    happen for a family the matrix itself carries, but tolerated) is
    excluded from both `accessory_present` and the core set rather than
    guessed at either way."""
    core_families = {f for f, b in freq_bin.items() if b == "core"}
    content: dict[str, dict] = {}
    for strain in strains:
        n_present = 0
        accessory_present = 0
        for fam in matrix.families:
            if matrix.is_present(fam, strain):
                n_present += 1
                b = freq_bin.get(fam)
                if b is not None and b != "core":
                    accessory_present += 1
        core_missing = sum(1 for fam in core_families if not matrix.is_present(fam, strain))
        content[strain] = {
            "n_families_present": n_present,
            "accessory_present": accessory_present,
            "core_missing": core_missing,
        }
    return content


def compute_terminus_enrichment(
    gene_positions_path: str | Path,
    contig_lengths_by_strain: dict[str, dict[str, int]],
    member_to_rep: dict[str, str],
    private_families_by_strain: dict[str, set[str]],
    id_sep: str,
    window_bp: int,
) -> dict[str, dict]:
    """{strain: {n_all_proteins, n_all_at_terminus, n_private_proteins,
    n_private_at_terminus}}, restricted to strains present in
    `contig_lengths_by_strain` (the ones this run has DNA + assembly stats
    for)."""
    counts = {
        s: {"n_all_proteins": 0, "n_all_at_terminus": 0,
            "n_private_proteins": 0, "n_private_at_terminus": 0}
        for s in contig_lengths_by_strain
    }
    with open_maybe_compressed(gene_positions_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            short = row["Short"]
            if short not in counts:
                continue
            contig_len = contig_lengths_by_strain[short].get(row["contig"])
            if contig_len is None:
                continue
            start, end = int(row["start"]), int(row["end"])
            at_terminus = is_at_terminus(start, end, contig_len, window_bp)
            counts[short]["n_all_proteins"] += 1
            if at_terminus:
                counts[short]["n_all_at_terminus"] += 1

            member_id = f"{short}{id_sep}{row['protein_id']}"
            family = member_to_rep.get(member_id)
            if family is not None and family in private_families_by_strain.get(short, ()):
                counts[short]["n_private_proteins"] += 1
                if at_terminus:
                    counts[short]["n_private_at_terminus"] += 1
    return counts


def _pct(numerator: int, denominator: int) -> str:
    return f"{100 * numerator / denominator:.2f}" if denominator else "-"


def build_per_strain_table(
    strains: list[str],
    assembly_stats: dict[str, dict],
    content: dict[str, dict],
    terminus_counts: dict[str, dict],
) -> list[dict]:
    rows = []
    for short in strains:
        stats = assembly_stats[short]
        c = content[short]
        t = terminus_counts.get(short, {
            "n_all_proteins": 0, "n_all_at_terminus": 0,
            "n_private_proteins": 0, "n_private_at_terminus": 0,
        })
        rows.append({
            "Short": short,
            "n_contigs": stats["n_contigs"],
            "n50": stats["n50"],
            "total_length": stats["total_length"],
            "n_families_present": c["n_families_present"],
            "accessory_present": c["accessory_present"],
            "core_missing": c["core_missing"],
            "n_private_proteins": t["n_private_proteins"],
            "pct_private_at_terminus": _pct(t["n_private_at_terminus"], t["n_private_proteins"]),
            "n_all_proteins": t["n_all_proteins"],
            "pct_all_at_terminus": _pct(t["n_all_at_terminus"], t["n_all_proteins"]),
        })
    return rows


PER_STRAIN_FIELDNAMES = [
    "Short", "n_contigs", "n50", "total_length", "n_families_present",
    "accessory_present", "core_missing", "n_private_proteins",
    "pct_private_at_terminus", "n_all_proteins", "pct_all_at_terminus",
]

CORRELATION_COMPARISONS = [
    ("accessory_present", "n_contigs"),
    ("accessory_present", "n50"),
    ("core_missing", "n_contigs"),
    ("core_missing", "n50"),
]


def build_correlation_table(per_strain_rows: list[dict]) -> list[dict]:
    """Spearman rho, raw and partial (controlling total_length), for each
    of CORRELATION_COMPARISONS. Rows with fewer than 3 strains report "-"
    for every rho (a 2-point Spearman/partial correlation is a degenerate
    fit, not a real estimate)."""
    n = len(per_strain_rows)
    total_length = [float(r["total_length"]) for r in per_strain_rows]
    out = []
    for y_name, x_name in CORRELATION_COMPARISONS:
        y = [float(r[y_name]) for r in per_strain_rows]
        x = [float(r[x_name]) for r in per_strain_rows]
        if n < 3:
            rho_raw, rho_partial = float("nan"), float("nan")
        else:
            rho_raw = spearman(y, x)
            rho_partial = partial_spearman(y, x, total_length)
        out.append({
            "comparison": f"{y_name}_vs_{x_name}",
            "y": y_name,
            "x": x_name,
            "rho_raw": f"{rho_raw:.4f}" if not math.isnan(rho_raw) else "-",
            "rho_partial_length": f"{rho_partial:.4f}" if not math.isnan(rho_partial) else "-",
            "n_strains": n,
        })
    return out


def render_report_markdown(
    correlation_rows: list[dict],
    per_strain_rows: list[dict],
    rho_warn_threshold: float,
    window_bp: int,
) -> tuple[str, bool]:
    """Returns (markdown_text, warning_triggered). Warns when |rho| (raw OR
    partial) for either accessory_present comparison exceeds
    `rho_warn_threshold` -- accessory-genome comparisons between strains are
    confounded by assembly quality on this dataset, matching issue #130's
    acceptance criterion."""
    lines = ["# Assembly Quality vs Pangenome Content QC", ""]
    lines += [
        "Fragmented assemblies can manufacture spurious \"accessory\" gene "
        "families: a gene broken across a contig boundary is predicted as "
        "two partial proteins, neither of which clusters with the intact "
        "family, so each becomes a fake strain-private family. This report "
        "checks whether that confound is present on this run.", "",
    ]

    lines += ["## Correlations: assembly quality vs pangenome content", ""]
    lines += [
        "| Comparison | rho (raw) | rho (partial, controlling total_length) | n strains |",
        "|---|---|---|---|",
    ]
    for row in correlation_rows:
        lines.append(
            f"| {row['y']} vs {row['x']} | {row['rho_raw']} | "
            f"{row['rho_partial_length']} | {row['n_strains']} |"
        )
    lines.append("")

    warning_triggered = False
    for row in correlation_rows:
        if row["y"] != "accessory_present":
            continue
        for key in ("rho_raw", "rho_partial_length"):
            if row[key] == "-":
                continue
            if abs(float(row[key])) > rho_warn_threshold:
                warning_triggered = True

    if warning_triggered:
        lines += [
            f"**WARNING: |rho| exceeds {rho_warn_threshold} for accessory-family "
            "content against assembly quality.** Between-strain accessory-genome "
            "comparisons (accessory island counts, gain/loss calls, Leiden "
            "trans-modules) are confounded on this dataset: \"strain X has more "
            "accessory genes\" may mean \"strain X has a worse assembly\", not a "
            "real biological difference. This is a diagnostic only -- no "
            "presence call or strain has been altered or excluded.", "",
        ]

    n_all = sum(int(r["n_all_proteins"]) for r in per_strain_rows)
    n_all_term = sum(
        round(float(r["pct_all_at_terminus"]) / 100 * int(r["n_all_proteins"]))
        if r["pct_all_at_terminus"] != "-" else 0
        for r in per_strain_rows
    )
    n_priv = sum(int(r["n_private_proteins"]) for r in per_strain_rows)
    n_priv_term = sum(
        round(float(r["pct_private_at_terminus"]) / 100 * int(r["n_private_proteins"]))
        if r["pct_private_at_terminus"] != "-" else 0
        for r in per_strain_rows
    )
    lines += ["## Contig-terminus enrichment (mechanism check)", ""]
    lines += [
        f"Proteins in strain-private (singleton) families sit within "
        f"{window_bp} bp of a contig end **{_pct(n_priv_term, n_priv)}%** of "
        f"the time (n={n_priv}), against **{_pct(n_all_term, n_all)}%** of "
        f"all {n_all} proteins genome-wide. See "
        "assembly_quality_vs_content.tsv for the per-strain breakdown.", "",
    ]

    return "\n".join(lines), warning_triggered


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, help="samplesheet CSV")
    ap.add_argument("--data_dir", required=True, help="directory holding dna/<Sample.DNA> FASTAs")
    ap.add_argument("--matrix", required=True, help="presence matrix TSV (rescued, if rescue is enabled)")
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--gene_positions", required=True, help="gene_positions.tsv[.gz|.zst]")
    ap.add_argument("--cluster_tsv", required=True, help="tier-1 mmseqs/diamond cluster TSV")
    ap.add_argument("--ingroup_label", default="IN",
                     help="samplesheet GROUP value identifying ingroup strains (default: IN)")
    ap.add_argument("--id_sep", default="|")
    ap.add_argument("--terminus_window_bp", type=int, default=1000,
                     help="a protein span within this many bp of a contig end counts as "
                          "'at a contig terminus' (default: 1000)")
    ap.add_argument("--rho_warn_threshold", type=float, default=0.3,
                     help="emit a WARNING (stderr + report.md) when |rho| for accessory-vs-"
                          "assembly-quality exceeds this (default: 0.3)")
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    samples = parse_config(args.config)
    ingroup = [s for s in samples if s.group == args.ingroup_label]
    dna_paths = {s.short: Path(args.data_dir) / "dna" / s.dna for s in ingroup if s.dna}

    matrix = PresenceMatrix.from_tsv(args.matrix)
    with open(args.frequency_table, newline="") as fh:
        freq_rows = list(csv.DictReader(fh, delimiter="\t"))
    freq_bin = {r["family"]: r["bin"] for r in freq_rows}
    singleton_families = {f for f, b in freq_bin.items() if b == "singleton"}
    member_to_rep = read_cluster_tsv(args.cluster_tsv)

    assembly_stats: dict[str, dict] = {}
    contig_lengths_by_strain: dict[str, dict[str, int]] = {}
    for s in ingroup:
        path = dna_paths.get(s.short)
        if path is None or not Path(path).exists():
            print(f"WARNING: no DNA FASTA for ingroup strain {s.short}; excluded from "
                  "assembly-quality QC", file=sys.stderr)
            continue
        lengths = parse_contig_lengths(path)
        contig_lengths_by_strain[s.short] = lengths
        assembly_stats[s.short] = assembly_stats_from_lengths(lengths)

    strains = [s.short for s in ingroup if s.short in assembly_stats and s.short in matrix.strains]
    missing_from_matrix = [
        s.short for s in ingroup if s.short in assembly_stats and s.short not in matrix.strains
    ]
    if missing_from_matrix:
        print(f"WARNING: {len(missing_from_matrix)} ingroup strains with assembly stats are "
              f"absent from the presence matrix and are excluded (first few: "
              f"{missing_from_matrix[:5]})", file=sys.stderr)

    content = compute_per_strain_content(matrix, freq_bin, strains)
    private_family_strain = find_private_family_strain(matrix, singleton_families)
    private_families_by_strain: dict[str, set[str]] = {}
    for fam, strain in private_family_strain.items():
        private_families_by_strain.setdefault(strain, set()).add(fam)

    terminus_counts = compute_terminus_enrichment(
        args.gene_positions, contig_lengths_by_strain, member_to_rep,
        private_families_by_strain, args.id_sep, args.terminus_window_bp,
    )

    per_strain_rows = build_per_strain_table(strains, assembly_stats, content, terminus_counts)
    write_tsv(out_dir / "assembly_quality_vs_content.tsv", PER_STRAIN_FIELDNAMES, per_strain_rows)

    correlation_rows = build_correlation_table(per_strain_rows)
    write_tsv(
        out_dir / "assembly_quality_correlations.tsv",
        ["comparison", "y", "x", "rho_raw", "rho_partial_length", "n_strains"],
        correlation_rows,
    )

    report_md, warning_triggered = render_report_markdown(
        correlation_rows, per_strain_rows, args.rho_warn_threshold, args.terminus_window_bp,
    )
    (out_dir / "assembly_quality_report.md").write_text(report_md)

    print(
        f"pangenome_assembly_quality_qc: wrote assembly_quality_vs_content.tsv "
        f"({len(per_strain_rows)} strains), assembly_quality_correlations.tsv, "
        f"assembly_quality_report.md to {out_dir}", file=sys.stderr,
    )
    if warning_triggered:
        print(
            f"WARNING: assembly-quality confound detected (|rho| > "
            f"{args.rho_warn_threshold} for accessory-family content vs assembly "
            "quality) -- see assembly_quality_report.md", file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
