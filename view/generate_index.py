#!/usr/bin/env python3
"""
Regenerate the view/ landing pages.

Two things are produced:

  1. A per-project ``report.html`` inside each sub-project folder — a run
     summary reconstructed **by parsing the report HTML files themselves**
     (novelties.html / core.html / losses.html). Each of those reports embeds
     its data as JSON in a ``<script id="payload">`` block, which carries the
     ingroup/outgroup proteomes and the run's presence-fraction parameters, so
     the summary needs neither the original config CSV nor a pipeline rerun —
     a folder of reports copied off the cluster is enough. Existing
     report.html files are left alone unless ``--force`` is given.

  2. The top-level ``view/index.html`` — a landing page linking to every
     sub-project's report.html and its individual reports.

Scans immediate subdirectories of this script's directory for *.html files.
Re-run any time a new sub-project report is added or an existing one is
regenerated.

Example:
  ./generate_index.py                 # write index.html + any missing report.html
  ./generate_index.py --force         # also regenerate every report.html
  ./generate_index.py --output index.html
"""
from __future__ import annotations

import argparse
import datetime
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'lib'))
from index_page import render_gallery_page, render_project_page  # noqa: E402

# The folder title links here.
REPORT_NAME = "report.html"

# Sub-links always shown (in this order) when the file is present. These are
# also the reports parsed to reconstruct each folder's run summary.
SUBPAGES = ["novelties.html", "core.html", "losses.html"]

# Human labels + the payload row-count meaning for each report type.
REPORT_LABELS = {
    "novelties.html": ("Novelty candidates", "lineage-specific genes"),
    "core.html": ("Core genes", "near-universally conserved"),
    "losses.html": ("Candidate gene losses", "conserved in outgroup, lost in ingroup"),
}

# Run parameters that may appear in the payloads, with display labels, in the
# order we want to show them. Different reports contribute different keys.
PARAM_LABELS = [
    ("ingroup_min_frac", "Ingroup min fraction"),
    ("outgroup_min_frac", "Outgroup min fraction"),
    ("loss_ingroup_max_frac", "Loss ingroup max fraction"),
    ("core_min_frac", "Core min fraction"),
]

_PAYLOAD_RE = re.compile(
    r'<script[^>]*id="payload"[^>]*>(.*?)</script>', re.S
)


def extract_payload(html_path: Path) -> dict | None:
    """Return the embedded JSON payload from a NovInvenio report, or None.

    The reports escape ``</`` as ``<\\/`` before embedding, so no literal
    ``</script>`` appears inside the JSON and the first one is the real closing
    tag. ``\\/`` is a valid JSON escape for ``/``, so json.loads decodes it
    without any un-escaping on our side.
    """
    try:
        text = html_path.read_text(encoding="utf-8")
    except OSError:
        return None
    m = _PAYLOAD_RE.search(text)
    if not m:
        return None
    body = m.group(1).strip()
    if not body or body == "/*__PAYLOAD__*/":
        return None
    try:
        return json.loads(body)
    except json.JSONDecodeError:
        return None


def collect_run_summary(project_dir: Path) -> dict | None:
    """Reconstruct a run summary for one folder from its report payloads.

    Returns None when no report in the folder carries a parseable payload.
    """
    proteomes: dict[str, dict] = {}   # short -> proteome dict (IN/OUT membership)
    params: dict[str, object] = {}
    counts: dict[str, int] = {}       # filename -> row count
    project = project_dir.name

    for name in SUBPAGES:
        payload = extract_payload(project_dir / name)
        if not payload:
            continue
        project = payload.get("project", project)
        for p in payload.get("proteomes", []):
            short = p.get("short")
            if short and short not in proteomes:
                proteomes[short] = p
        for key, _ in PARAM_LABELS:
            if key in payload and key not in params:
                params[key] = payload[key]
        counts[name] = len(payload.get("rows", []))

    if not proteomes and not counts:
        return None

    ingroup = [p for p in proteomes.values() if p.get("group") == "IN"]
    outgroup = [p for p in proteomes.values() if p.get("group") == "OUT"]
    ingroup.sort(key=lambda p: p.get("short", ""))
    outgroup.sort(key=lambda p: p.get("short", ""))
    return {
        "project": project,
        "ingroup": ingroup,
        "outgroup": outgroup,
        "params": params,
        "counts": counts,
    }


def _report_cards(project_dir: Path, counts: dict[str, int]) -> list[dict]:
    cards = []
    for name in SUBPAGES:
        if not (project_dir / name).exists():
            continue
        label, blurb = REPORT_LABELS.get(name, (name, ""))
        n = counts.get(name)
        cards.append({
            "href": name,
            "title": label,
            "desc": blurb,
            "meta": f"{n:,} rows" if n is not None else "",
        })
    return cards


def _param_tiles(params: dict[str, object]) -> list[tuple[str, str]]:
    tiles = []
    for key, label in PARAM_LABELS:
        if key not in params:
            continue
        value = params[key]
        tiles.append((label, f"{value:g}" if isinstance(value, (int, float)) else str(value)))
    return tiles


def render_report(project_dir: Path, summary: dict) -> str:
    return render_project_page(
        project=summary["project"],
        reports=_report_cards(project_dir, summary["counts"]),
        tiles=(
            [
                ("Ingroup proteomes", str(len(summary["ingroup"]))),
                ("Outgroup proteomes", str(len(summary["outgroup"]))),
            ]
            + _param_tiles(summary["params"])
        ),
        ingroup=summary["ingroup"],
        outgroup=summary["outgroup"],
        note="Reconstructed from the data embedded in the reports in this folder — "
             "no config CSV or pipeline rerun needed.",
        footer=f"Generated {datetime.datetime.now():%Y-%m-%d %H:%M} by generate_index.py",
    )


def write_reports(root: Path, force: bool) -> None:
    """Write a report.html summary into each project folder.

    Skips folders whose report.html already exists unless *force* is set — a
    report.html generated by the pipeline (bin/make_index_report.py) is not
    clobbered by default.
    """
    for d in find_project_dirs(root):
        out = d / REPORT_NAME
        if out.exists() and not force:
            print(f"Skip {out} (exists; use --force to overwrite)")
            continue
        summary = collect_run_summary(d)
        if not summary:
            print(f"Skip {d.name}/ (no parseable report payloads)")
            continue
        out.write_text(render_report(d, summary), encoding="utf-8")
        print(f"Wrote {out}")


# ---- top-level view/index.html -------------------------------------------

def find_project_dirs(root: Path) -> list[Path]:
    return sorted(
        p for p in root.iterdir()
        if p.is_dir() and any(p.glob("*.html"))
    )


def project_title(dirname: str) -> str:
    return dirname.replace("_", " ").replace("-", " ")


def page_title(filename: str) -> str:
    return filename[:-len(".html")].replace("_", " ").replace("-", " ")


def _project_card(d: Path) -> dict:
    """One gallery card, with the run's actual shape rather than just its name.

    The species and row counts are already sitting in the report payloads (see
    collect_run_summary) -- the old gallery parsed none of it and showed a bare
    directory name and an mtime, which told a first-time visitor nothing.
    """
    report = d / REPORT_NAME
    main = report if report.exists() else sorted(d.glob("*.html"))[0]
    summary = collect_run_summary(d)

    desc = ""
    meta_bits = []
    if summary:
        n_in, n_out = len(summary["ingroup"]), len(summary["outgroup"])
        if n_in or n_out:
            desc = f"{n_in} ingroup · {n_out} outgroup proteomes"
        counts = summary["counts"]
        for name in SUBPAGES:
            if name in counts:
                label = REPORT_LABELS.get(name, (name, ""))[0]
                meta_bits.append(f"{counts[name]:,} {label.lower()}")
    mtime = datetime.datetime.fromtimestamp(main.stat().st_mtime)
    meta_bits.append(f"updated {mtime:%Y-%m-%d}")

    return {
        "title": project_title(d.name),
        "href": f"{d.name}/{main.name}",
        "desc": desc,
        "meta": " · ".join(meta_bits),
        "subpages": [
            {"href": f"{d.name}/{name}", "title": page_title(name)}
            for name in SUBPAGES if (d / name).exists()
        ],
    }


def render(root: Path) -> str:
    return render_gallery_page(
        projects=[_project_card(d) for d in find_project_dirs(root)],
        footer=f"Regenerated {datetime.datetime.now():%Y-%m-%d %H:%M} by generate_index.py",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--output", default=Path(__file__).parent / "index.html", type=Path,
        help="Path to write the generated top-level index (default: index.html next to this script)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Overwrite existing per-folder report.html files (default: only create missing ones)",
    )
    args = parser.parse_args()

    root = Path(__file__).parent
    # Build per-folder report.html summaries first, so the index links to them.
    write_reports(root, args.force)
    args.output.write_text(render(root))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
