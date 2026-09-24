#!/usr/bin/env python3
"""Surface pipeline diagnostics (issue #134) where a reader cannot skip them.

The precedent this closes: RESCUE_PASS's per-file funnel counts (rows
parsed, rows passed threshold, structural-filter rejections, cells applied)
existed only as stderr prints. A warning that only exists in a log does not
exist -- a 530-strain study ran to "SUCCESS" while every downstream table was
computed on unrescued data, and nobody noticed until a review went looking.

Advisory by default: a triggered diagnostic is reported, never fails the
run. `--pangenome_strict` promotes any triggered diagnostic to a run
failure -- see main()'s docstring for the one exception (the zero-hit guard,
issue #126, which already always fails inside pangenome_rescue_pass.py and
is not re-implemented here).

Diagnostics wired to REAL computed data in this pipeline today:
  - rescue_redundancy (issue #133): from pangenome_rescue_pass.py's funnel-
    stats sidecar (--rescue_funnel), written when RESCUE_PASS is given
    --funnel_tsv.
  - assembly_quality_confound (issue #130): from ASSEMBLY_QUALITY_QC's
    assembly_quality_correlations.tsv (--assembly_correlations). Same trip
    rule as that script's own report warning, so the two never disagree.

Diagnostics named in issue #134's table that this pipeline does not yet
compute the underlying statistic for (PCoA/split-half-ARI clade validity,
gain/loss event classification) are reported with status=not_computed rather than
fabricated -- see not_computed_diagnostic()'s docstring. Wiring a real one
in means adding an `evaluate_*` function here fed by that statistic's own
pipeline output, the same way evaluate_rescue_redundancy() is fed by
pangenome_rescue_pass.py's sidecar.

Usage:
  pangenome_diagnostics.py --rescue_funnel rescue_funnel.tsv \\
      --assembly_correlations assembly_quality_correlations.tsv \\
      --out_dir diagnostics/ [--pangenome_strict]

Writes, into --out_dir:
  diagnostics.tsv          machine-readable (one row per diagnostic)
  diagnostics_banner.md    Markdown block, prepend-ready for report.md
  diagnostics_banner.html  HTML snippet, insert-ready into a report page
"""
from __future__ import annotations

import argparse
import csv
import html as html_module
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Diagnostic:
    diagnostic_id: str
    status: str            # "triggered" | "ok" | "not_computed"
    would_fail_strict: bool
    detail: str
    prune_options: list[str]


def read_funnel_tsv(path: "str | Path") -> dict[str, int]:
    """Reads a pangenome_rescue_pass.py --funnel_tsv sidecar (metric\\tvalue
    rows) into {metric: int value}."""
    funnel: dict[str, int] = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)  # header
        for row in reader:
            if len(row) != 2:
                continue
            metric, value = row
            funnel[metric] = int(value)
    return funnel


def evaluate_rescue_redundancy(funnel: dict[str, int], threshold: float = 0.5) -> Diagnostic:
    """Issue #133/#134: trips when more than `threshold` (default 50%) of
    identity/coverage-qualifying ("rescuable") rescue-pass hits are rejected
    by the structural filter's overlap check -- the hit lands on a predicted
    gene already assigned to a DIFFERENT family in that strain, i.e. the same
    locus double-counted rather than a genuine annotation dropout.

    `rows_passed_threshold` == 0 (no rescuable cells at all, e.g. a study
    with zero ABSENT calls or rescue disabled upstream) is reported `ok`, not
    a fabricated 0/0 fraction."""
    passed = funnel.get("rows_passed_threshold", 0)
    overlap = funnel.get("rows_rejected_overlap", 0)
    if passed == 0:
        return Diagnostic(
            "rescue_redundancy", "ok", False,
            "0 rescuable (identity/coverage-passing) rescue-pass hits -- "
            "nothing to assess.",
            [],
        )
    fraction = overlap / passed
    detail = (
        f"{overlap}/{passed} ({fraction * 100:.1f}%) rescuable cells overlap "
        "a predicted gene already assigned to a different family in that "
        "strain -- the same locus double-counted, not a genuine annotation "
        "dropout (issue #133)."
    )
    if fraction > threshold:
        return Diagnostic(
            "rescue_redundancy", "triggered", True, detail,
            [
                "do not enable rescue (--pangenome_rescue_enable false); "
                "see clustering sweep #132 before re-enabling",
                "verify --pangenome_rescue_structural_filter is enabled "
                "(default true) -- it is the mechanism that removes these "
                "cells, not merely detects them",
            ],
        )
    return Diagnostic("rescue_redundancy", "ok", False, detail, [])


def read_correlations_tsv(path: "str | Path") -> list[dict[str, str]]:
    """Reads pangenome_assembly_quality_qc.py's
    assembly_quality_correlations.tsv into a list of row dicts (values kept
    as strings; an undefined rho is the literal "-")."""
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def evaluate_assembly_quality_confound(rows: list[dict[str, str]], threshold: float = 0.3) -> Diagnostic:
    """Issue #130/#134: trips when |rho| (raw OR partial, controlling
    total_length) exceeds `threshold` for any accessory_present comparison
    (vs n_contigs or n50). This is the same rule
    pangenome_assembly_quality_qc.py uses for its own report WARNING, and
    the default threshold matches that script's --rho_warn_threshold
    default (0.3). core_missing rows are ignored here, as they are there.

    Every accessory_present rho undefined ("-", fewer than 3 strains) is
    reported `ok` with that stated -- never a fabricated correlation."""
    worst_name, worst_value = None, None
    for row in rows:
        if row.get("y") != "accessory_present":
            continue
        for key in ("rho_raw", "rho_partial_length"):
            value = row.get(key, "-")
            if value in ("-", ""):
                continue
            rho = float(value)
            if worst_value is None or abs(rho) > abs(worst_value):
                label = "raw" if key == "rho_raw" else "partial, controlling total_length"
                worst_name, worst_value = f"accessory_present vs {row['x']} ({label})", rho
    if worst_value is None:
        return Diagnostic(
            "assembly_quality_confound", "ok", False,
            "every accessory_present rho is undefined (fewer than 3 strains) "
            "-- nothing to assess.",
            [],
        )
    n_strains = next((r.get("n_strains") for r in rows if r.get("y") == "accessory_present"), "?")
    detail = (
        f"max |rho| = {abs(worst_value):.2f} (rho = {worst_value:.4f}, "
        f"{worst_name}, n = {n_strains} strains); threshold {threshold}. "
        "See assembly_quality_report.md (issue #130)."
    )
    if abs(worst_value) > threshold:
        return Diagnostic(
            "assembly_quality_confound", "triggered", True,
            detail + " Between-strain accessory-content comparisons are "
            "confounded by assembly fragmentation on this dataset.",
            [
                "exclude or flag the most fragmented assemblies (highest "
                "n_contigs in assembly_quality_vs_content.tsv) and rerun",
                "report accessory-content comparisons with assembly quality "
                "as a covariate, not as raw per-strain counts",
            ],
        )
    return Diagnostic("assembly_quality_confound", "ok", False, detail, [])


def not_computed_diagnostic(diagnostic_id: str, tracking_issue: str, description: str) -> Diagnostic:
    """A diagnostic named in issue #134's table whose underlying statistic
    this pipeline does not yet compute anywhere (no N50/PCoA/gain-loss
    module exists to feed it). Reported honestly as `not_computed` --
    never a fabricated pass/fail -- so the reader knows the check was
    DECLARED, not silently skipped, and can find the tracking issue that
    will wire in the real computation."""
    return Diagnostic(
        diagnostic_id, "not_computed", False,
        f"trips when: {description}. Not yet computed by this pipeline -- "
        f"see {tracking_issue}.",
        [],
    )


def default_not_computed_diagnostics() -> list[Diagnostic]:
    """The issue #134 diagnostics with no computed input in this pipeline
    yet (see module docstring)."""
    return [
        not_computed_diagnostic(
            "clade_structure_validity", "#131",
            "PCoA1 variance share < ~50%, split-half ARI < 0.8, or "
            "between/within clade distance ratio < ~2",
        ),
        not_computed_diagnostic(
            "gain_loss_skew", "issue #134's calibration reference",
            "gain:loss ratio > 20:1 or ambiguous-polarity fraction < 2%",
        ),
    ]


def render_banner_markdown(diagnostics: list[Diagnostic]) -> str:
    """A '## Pipeline diagnostics' Markdown section listing every
    diagnostic's status -- rendered even when nothing triggered, so a clean
    run is visibly CHECKED, not silently unchecked."""
    lines = ["## Pipeline diagnostics", ""]
    for d in diagnostics:
        label = {"triggered": "TRIGGERED", "ok": "OK", "not_computed": "not computed"}[d.status]
        lines.append(f"- **{d.diagnostic_id}** [{label}]: {d.detail}")
        for opt in d.prune_options:
            lines.append(f"  - prune option: {opt}")
    lines.append("")
    return "\n".join(lines)


def render_banner_html(diagnostics: list[Diagnostic]) -> str:
    """Self-contained HTML snippet (no external CSS/JS dependency) with the
    same content as render_banner_markdown(), for embedding at the top of a
    report HTML page. Every diagnostic-supplied string is html.escape()'d --
    this pipeline does not control what a future diagnostic's detail text
    might quote from upstream data."""
    status_class = {"triggered": "diag-triggered", "ok": "diag-ok", "not_computed": "diag-not-computed"}
    status_label = {"triggered": "TRIGGERED", "ok": "OK", "not_computed": "not computed"}
    items = []
    for d in diagnostics:
        opts = "".join(
            f"<li>prune option: {html_module.escape(opt)}</li>" for opt in d.prune_options
        )
        opts_html = f"<ul>{opts}</ul>" if opts else ""
        items.append(
            f'<li class="{status_class[d.status]}">'
            f"<strong>{html_module.escape(d.diagnostic_id)}</strong> "
            f"[{status_label[d.status]}]: {html_module.escape(d.detail)}"
            f"{opts_html}</li>"
        )
    return (
        '<section id="pangenome-diagnostics" '
        'style="border:1px solid #999;padding:0.75em;margin-bottom:1em;">'
        "<h2>Pipeline diagnostics</h2><ul>" + "".join(items) + "</ul></section>"
    )


def write_diagnostics_tsv(diagnostics: list[Diagnostic], path: "str | Path") -> None:
    with open(path, "w", newline="") as out:
        out.write("diagnostic_id\tstatus\twould_fail_strict\tdetail\tprune_options\n")
        for d in diagnostics:
            prune = ";".join(d.prune_options)
            out.write(f"{d.diagnostic_id}\t{d.status}\t{d.would_fail_strict}\t{d.detail}\t{prune}\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--rescue_funnel", default=None,
        help="pangenome_rescue_pass.py --funnel_tsv sidecar. Omit when "
        "--pangenome_rescue_enable is false (rescue_redundancy is then "
        "reported not_computed rather than erroring).",
    )
    ap.add_argument(
        "--pangenome_strict", action="store_true", default=False,
        help="Promote any TRIGGERED diagnostic to a run failure (exit 1). "
        "Advisory (warn-only) by default. The zero-hit rescue guard (issue "
        "#126) always fails regardless of this flag -- that is enforced "
        "inside pangenome_rescue_pass.py itself, upstream of this script, "
        "and is not re-implemented here.",
    )
    ap.add_argument(
        "--rescue_redundancy_threshold", type=float, default=0.5,
        help="Fraction of rescuable cells rejected as overlapping a "
        "different family's gene above which rescue_redundancy triggers "
        "(default 0.5).",
    )
    ap.add_argument(
        "--assembly_correlations", default=None,
        help="pangenome_assembly_quality_qc.py's "
        "assembly_quality_correlations.tsv. Omit to report "
        "assembly_quality_confound as not_computed.",
    )
    ap.add_argument(
        "--assembly_rho_threshold", type=float, default=0.3,
        help="|rho| above which assembly_quality_confound triggers "
        "(default 0.3, same as pangenome_assembly_quality_qc.py's "
        "--rho_warn_threshold default).",
    )
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    diagnostics: list[Diagnostic] = []
    if args.rescue_funnel:
        funnel = read_funnel_tsv(args.rescue_funnel)
        diagnostics.append(evaluate_rescue_redundancy(funnel, threshold=args.rescue_redundancy_threshold))
    else:
        diagnostics.append(not_computed_diagnostic(
            "rescue_redundancy", "issue #133",
            ">50% of rescuable cells overlap another family's gene",
        ))
    if args.assembly_correlations:
        diagnostics.append(evaluate_assembly_quality_confound(
            read_correlations_tsv(args.assembly_correlations),
            threshold=args.assembly_rho_threshold,
        ))
    else:
        diagnostics.append(not_computed_diagnostic(
            "assembly_quality_confound", "#130",
            "abs(rho) between accessory-family count and assembly N50 > ~0.3",
        ))
    diagnostics += default_not_computed_diagnostics()

    write_diagnostics_tsv(diagnostics, out_dir / "diagnostics.tsv")
    (out_dir / "diagnostics_banner.md").write_text(render_banner_markdown(diagnostics))
    (out_dir / "diagnostics_banner.html").write_text(render_banner_html(diagnostics))

    for d in diagnostics:
        print(f"diagnostic {d.diagnostic_id}: {d.status} -- {d.detail}", file=sys.stderr)

    triggered = [d for d in diagnostics if d.status == "triggered"]
    if args.pangenome_strict and triggered:
        names = ", ".join(d.diagnostic_id for d in triggered)
        print(
            f"ERROR: --pangenome_strict is set and {len(triggered)} "
            f"diagnostic(s) triggered: {names}. See {out_dir / 'diagnostics.tsv'}.",
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
