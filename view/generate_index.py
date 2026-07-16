#!/usr/bin/env python3
"""
Regenerate view/index.html: a landing page linking to each sub-project's
report.html.

Scans immediate subdirectories of this script's directory for a report.html
file and lists them as links, sorted alphabetically. Re-run this any time a
new sub-project report is added or an existing one is regenerated.

Example:
  ./generate_index.py
  ./generate_index.py --output index.html
"""
import argparse
import datetime
import html
from pathlib import Path

REPORT_NAME = "report.html"

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
    <div class="meta">updated {mtime}</div>
  </li>"""


def find_reports(root: Path) -> list[Path]:
    return sorted(
        p for p in root.iterdir()
        if p.is_dir() and (p / REPORT_NAME).is_file()
    )


def project_title(dirname: str) -> str:
    return dirname.replace("_", " ").replace("-", " ")


def render(root: Path) -> str:
    project_dirs = find_reports(root)
    items = []
    for d in project_dirs:
        report = d / REPORT_NAME
        mtime = datetime.datetime.fromtimestamp(report.stat().st_mtime)
        items.append(ITEM_TEMPLATE.format(
            href=f"{d.name}/{REPORT_NAME}",
            title=html.escape(project_title(d.name)),
            mtime=mtime.strftime("%Y-%m-%d"),
        ))
    generated = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    return PAGE_TEMPLATE.format(items="\n".join(items), generated=generated)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", default=Path(__file__).parent / "index.html", type=Path,
        help="Path to write the generated index (default: index.html next to this script)",
    )
    args = parser.parse_args()

    root = Path(__file__).parent
    args.output.write_text(render(root))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
