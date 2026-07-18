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
import html
import json
import re
from pathlib import Path

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


REPORT_PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>__PROJECT__ — NovInvenio run summary</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
  body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
    max-width: 60rem; margin: 2rem auto; padding: 0 1rem; color: #1a1a1a; background: #fff;
  }
  h1 { font-size: 1.5rem; margin-bottom: 0.25rem; }
  h2 { font-size: 1.1rem; margin-top: 1.75rem; }
  .sub { color: #666; font-size: 0.9rem; margin-top: 0; }
  a { color: #0645ad; text-decoration: none; }
  a:hover { text-decoration: underline; }
  .reports { list-style: none; padding: 0; }
  .reports li {
    border: 1px solid #ddd; border-radius: 6px; padding: 0.6rem 0.9rem; margin-bottom: 0.5rem;
  }
  .reports .meta { color: #666; font-size: 0.85rem; }
  table { border-collapse: collapse; width: 100%; font-size: 0.9rem; margin-top: 0.5rem; }
  th, td { text-align: left; padding: 0.35rem 0.6rem; border-bottom: 1px solid #eee; }
  th { border-bottom: 1px solid #bbb; }
  .mono { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
  .groups { display: grid; grid-template-columns: 1fr 1fr; gap: 1.25rem; }
  @media (max-width: 760px) { .groups { grid-template-columns: 1fr; } }
  .params { list-style: none; padding: 0; display: flex; flex-wrap: wrap; gap: 1.5rem; }
  .params .label { color: #666; font-size: 0.8rem; }
  .params .value { font-size: 1.1rem; font-weight: 600; }
  footer { margin-top: 2rem; color: #999; font-size: 0.8rem; }
  @media (prefers-color-scheme: dark) {
    body { background: #1a1a1a; color: #ddd; }
    .reports li { border-color: #444; }
    th { border-color: #555; } td { border-color: #333; }
    a { color: #6cb2ff; } .sub, .reports .meta, .params .label { color: #999; }
    footer { color: #777; }
  }
</style>
</head>
<body>
<h1>__PROJECT__</h1>
<p class="sub">NovInvenio run summary · reconstructed from the reports in this folder</p>

<h2>Reports</h2>
<ul class="reports">
__REPORTS__
</ul>

<h2>Run parameters</h2>
<ul class="params">
__PARAMS__
</ul>

<div class="groups">
  <div>
    <h2>Ingroup (__N_IN__)</h2>
    __INGROUP__
  </div>
  <div>
    <h2>Outgroup (__N_OUT__)</h2>
    __OUTGROUP__
  </div>
</div>

<footer>Generated __GENERATED__ by generate_index.py</footer>
</body>
</html>
"""


def _species_table(proteomes: list[dict]) -> str:
    if not proteomes:
        return '<p class="sub">None found in the reports.</p>'
    rows = []
    for p in proteomes:
        name = html.escape(p.get("species", ""))
        if p.get("strain"):
            name += " " + html.escape(p["strain"])
        rows.append(
            "<tr><td class=\"mono\">{short}</td><td>{name}</td><td>{taxon}</td></tr>".format(
                short=html.escape(p.get("short", "")),
                name=name,
                taxon=html.escape(p.get("taxon", "")),
            )
        )
    return (
        '<table><thead><tr><th>Short</th><th>Species</th><th>Taxon group</th></tr></thead>'
        "<tbody>" + "\n".join(rows) + "</tbody></table>"
    )


def _reports_list(project_dir: Path, counts: dict[str, int]) -> str:
    items = []
    for name in SUBPAGES:
        if not (project_dir / name).exists():
            continue
        label, blurb = REPORT_LABELS.get(name, (name, ""))
        n = counts.get(name)
        meta = blurb
        if n is not None:
            meta = f"{n:,} rows · {blurb}" if blurb else f"{n:,} rows"
        items.append(
            '  <li><a href="{href}">{label}</a><div class="meta">{meta}</div></li>'.format(
                href=html.escape(name), label=html.escape(label), meta=html.escape(meta),
            )
        )
    return "\n".join(items) or '  <li class="meta">No reports found in this folder.</li>'


def _params_list(params: dict[str, object]) -> str:
    items = []
    for key, label in PARAM_LABELS:
        if key not in params:
            continue
        value = params[key]
        text = f"{value:g}" if isinstance(value, (int, float)) else str(value)
        items.append(
            '  <li><div class="label">{label}</div><div class="value">{value}</div></li>'.format(
                label=html.escape(label), value=html.escape(text),
            )
        )
    return "\n".join(items) or '  <li class="label">No run parameters recorded in the reports.</li>'


def render_report(project_dir: Path, summary: dict) -> str:
    generated = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    return (
        REPORT_PAGE_TEMPLATE
        .replace("__PROJECT__", html.escape(summary["project"]))
        .replace("__REPORTS__", _reports_list(project_dir, summary["counts"]))
        .replace("__PARAMS__", _params_list(summary["params"]))
        .replace("__INGROUP__", _species_table(summary["ingroup"]))
        .replace("__OUTGROUP__", _species_table(summary["outgroup"]))
        .replace("__N_IN__", str(len(summary["ingroup"])))
        .replace("__N_OUT__", str(len(summary["outgroup"])))
        .replace("__GENERATED__", generated)
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

PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>NovInvenio Reports</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
    max-width: 48rem;
    margin: 2rem auto;
    padding: 0 1rem;
    color: #1a1a1a;
    background: #fff;
  }}
  h1 {{ font-size: 1.5rem; }}
  ul {{ list-style: none; padding: 0; }}
  li {{
    border: 1px solid #ddd;
    border-radius: 6px;
    padding: 0.75rem 1rem;
    margin-bottom: 0.5rem;
  }}
  a {{ font-size: 1.1rem; text-decoration: none; color: #0645ad; }}
  a:hover {{ text-decoration: underline; }}
  .meta {{ color: #666; font-size: 0.85rem; margin-top: 0.2rem; }}
  .pages {{ list-style: none; padding: 0; margin: 0.5rem 0 0; }}
  .pages li {{
    border: none;
    padding: 0.15rem 0;
    margin: 0;
  }}
  .pages a {{ font-size: 0.95rem; }}
  footer {{ margin-top: 2rem; color: #999; font-size: 0.8rem; }}
  @media (prefers-color-scheme: dark) {{
    body {{ background: #1a1a1a; color: #ddd; }}
    li {{ border-color: #444; }}
    a {{ color: #6cb2ff; }}
    .meta {{ color: #999; }}
    footer {{ color: #777; }}
  }}
</style>
</head>
<body>
<h1>NovInvenio Reports</h1>
<ul>
{items}
</ul>
<footer>Regenerated {generated} by generate_index.py</footer>
</body>
</html>
"""

ITEM_TEMPLATE = """  <li>
    <a href="{href}">{title}</a>
    <div class="meta">updated {mtime}</div>{pages}
  </li>"""

PAGE_ITEM_TEMPLATE = """      <li><a href="{href}">{title}</a></li>"""

PAGES_LIST_TEMPLATE = """
    <ul class="pages">
{pages}
    </ul>"""


def find_project_dirs(root: Path) -> list[Path]:
    return sorted(
        p for p in root.iterdir()
        if p.is_dir() and any(p.glob("*.html"))
    )


def project_title(dirname: str) -> str:
    return dirname.replace("_", " ").replace("-", " ")


def page_title(filename: str) -> str:
    return filename[:-len(".html")].replace("_", " ").replace("-", " ")


def render(root: Path) -> str:
    project_dirs = find_project_dirs(root)
    items = []
    for d in project_dirs:
        # The folder title links to report.html when present; otherwise fall
        # back to the first available HTML page so the title is always a link.
        report = d / REPORT_NAME
        main = report if report.exists() else sorted(d.glob("*.html"))[0]

        subpages = [name for name in SUBPAGES if (d / name).exists()]
        pages_html = ""
        if subpages:
            page_items = "\n".join(
                PAGE_ITEM_TEMPLATE.format(
                    href=f"{d.name}/{name}",
                    title=html.escape(page_title(name)),
                )
                for name in subpages
            )
            pages_html = PAGES_LIST_TEMPLATE.format(pages=page_items)

        mtime = datetime.datetime.fromtimestamp(main.stat().st_mtime)
        items.append(ITEM_TEMPLATE.format(
            href=f"{d.name}/{main.name}",
            title=html.escape(project_title(d.name)),
            mtime=mtime.strftime("%Y-%m-%d"),
            pages=pages_html,
        ))
    generated = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    return PAGE_TEMPLATE.format(items="\n".join(items), generated=generated)


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
