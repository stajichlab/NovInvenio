#!/usr/bin/env python3
"""
Render a self-contained landing page (report.html) that links to the three
NovInvenio result reports — novelties.html, core.html, losses.html — and
records what job produced them: the ingroup and outgroup proteomes, the search
tool, and the presence-fraction thresholds.

Designed to sit in view/<project>/ next to copies of the three reports, so a
whole result set can be opened from one file:// URL and shared as a folder.
No network access is required to open any of the pages.

Example:
  make_index_report.py \
      --config configs/pezio4_asco.csv \
      --project pezio4_asco \
      --run_tool phmmer \
      --output view/pezio4_asco/report.html
"""
import argparse
import html
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import parse_config  # noqa: E402
from report_common import BASE_PAGE_CSS, THEME_TOGGLE_JS, THEME_VARS_CSS  # noqa: E402
from report_data import INGROUP_ROLES, OUTGROUP_ROLES  # noqa: E402

# Each entry: (relative filename, title, one-line description).
REPORTS = [
    ('novelties.html', 'Novelty candidates',
     'Lineage-specific genes: present across the ingroup, absent from every outgroup proteome.'),
    ('core.html', 'Core genes',
     'Near-universally conserved genes shared by essentially every proteome, ingroup and outgroup.'),
    ('losses.html', 'Candidate gene losses',
     'Gene families conserved across the outgroups but (nearly) absent from the ingroup.'),
]


def _species_rows(samples, roles):
    cells = []
    for s in samples:
        if s.group not in roles:
            continue
        name = html.escape(s.species + (f' {s.strain}' if s.strain else ''))
        cells.append(
            f'<tr><td class="mono">{html.escape(s.short)}</td>'
            f'<td>{name}</td><td>{html.escape(s.taxon_group or "")}</td></tr>'
        )
    return '\n'.join(cells)


def _report_cards(available):
    cards = []
    for fname, title, desc in REPORTS:
        if fname not in available:
            continue
        cards.append(
            f'<a class="report-card" href="{html.escape(fname)}">'
            f'<span class="report-title">{html.escape(title)}</span>'
            f'<span class="report-desc">{html.escape(desc)}</span>'
            f'<span class="report-file mono">{html.escape(fname)}</span></a>'
        )
    return '\n'.join(cards)


def _param_tiles(args, n_in, n_out):
    tiles = [
        ('Ingroup proteomes', str(n_in)),
        ('Outgroup proteomes', str(n_out)),
    ]
    if args.run_tool:
        tiles.append(('Search tool', args.run_tool))
    if args.ingroup_min_frac is not None:
        tiles.append(('Ingroup min fraction', f'{args.ingroup_min_frac:g}'))
    if args.outgroup_min_frac is not None:
        tiles.append(('Outgroup min fraction', f'{args.outgroup_min_frac:g}'))
    if args.loss_ingroup_max_frac is not None:
        tiles.append(('Loss ingroup max fraction', f'{args.loss_ingroup_max_frac:g}'))
    if args.core_min_frac is not None:
        tiles.append(('Core min fraction', f'{args.core_min_frac:g}'))
    return '\n'.join(
        f'<div><div class="tile-label">{html.escape(label)}</div>'
        f'<div class="tile-value">{html.escape(value)}</div></div>'
        for label, value in tiles
    )


PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__PROJECT_TITLE__ — NovInvenio reports</title>
<style>
""" + THEME_VARS_CSS + BASE_PAGE_CSS + r"""
  .reports { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 14px; }
  a.report-card {
    display: flex; flex-direction: column; gap: 6px; text-decoration: none;
    background: var(--surface-1); border: 1px solid var(--border); border-radius: 10px; padding: 16px;
    color: var(--text-primary);
  }
  a.report-card:hover { border-color: var(--series-1); background: var(--hover-wash); }
  .report-title { font-size: 15px; font-weight: 600; }
  .report-desc { font-size: 12px; color: var(--text-secondary); }
  .report-file { font-size: 11px; color: var(--muted); }
  .mono { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
  .groups { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
  @media (max-width: 800px) { .groups { grid-template-columns: 1fr; } }
</style>
</head>
<body>
<div class="wrap">
  <header class="top">
    <div class="titles">
      <h1>__PROJECT_TITLE__ — NovInvenio reports</h1>
      <p class="sub">Lineage-specific gene analysis · one folder, three interactive reports</p>
    </div>
    <button id="theme-toggle" class="btn-ghost" type="button" aria-label="Toggle colour theme">Theme</button>
  </header>

  <section class="card">
    <h2 class="card-title">Reports</h2>
    <p class="card-note">Each page is self-contained — open it directly, no network needed.</p>
    <div class="reports">
__REPORT_CARDS__
    </div>
  </section>

  <section class="card">
    <h2 class="card-title">Run summary</h2>
    <p class="card-note">Parameters and proteomes used to generate these reports.</p>
    <div class="tiles">
__PARAM_TILES__
    </div>
  </section>

  <div class="groups">
    <section class="card">
      <h2 class="card-title">Ingroup (__N_IN__)</h2>
      <p class="card-note">Genes are called lineage-specific relative to this clade.</p>
      <table class="data">
        <thead><tr><th scope="col">Short</th><th scope="col">Species</th><th scope="col">Taxon group</th></tr></thead>
        <tbody>
__INGROUP_ROWS__
        </tbody>
      </table>
    </section>
    <section class="card">
      <h2 class="card-title">Outgroup (__N_OUT__)</h2>
      <p class="card-note">Background proteomes used to exclude broadly conserved genes.</p>
      <table class="data">
        <thead><tr><th scope="col">Short</th><th scope="col">Species</th><th scope="col">Taxon group</th></tr></thead>
        <tbody>
__OUTGROUP_ROWS__
        </tbody>
      </table>
    </section>
  </div>
</div>
<script>
(function () {
  "use strict";
""" + THEME_TOGGLE_JS + r"""
})();
</script>
</body>
</html>
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', required=True, help='Analysis config CSV')
    ap.add_argument('--project', default=None, help='Project name shown in the title')
    ap.add_argument('--output', required=True, help='Output report.html')
    ap.add_argument('--reports_dir', default=None,
                    help='Directory the linked reports live in; a report is linked only if '
                         'its file is present there (default: the --output parent directory)')
    ap.add_argument('--run_tool', default=None, help='Search tool used (phmmer/diamond/blast)')
    ap.add_argument('--ingroup_min_frac', type=float, default=None)
    ap.add_argument('--outgroup_min_frac', type=float, default=None)
    ap.add_argument('--loss_ingroup_max_frac', type=float, default=None)
    ap.add_argument('--core_min_frac', type=float, default=None)
    args = ap.parse_args()

    samples = parse_config(args.config)
    project = args.project or Path(args.config).stem

    out = Path(args.output)
    reports_dir = Path(args.reports_dir) if args.reports_dir else out.parent
    available = {fname for fname, _, _ in REPORTS if (reports_dir / fname).exists()}
    # If none are present yet (e.g. the collate step stages them alongside this
    # output later), link all three by their conventional names.
    if not available:
        available = {fname for fname, _, _ in REPORTS}

    n_in = sum(1 for s in samples if s.group in INGROUP_ROLES)
    n_out = sum(1 for s in samples if s.group in OUTGROUP_ROLES)

    doc = (PAGE
           .replace('__PROJECT_TITLE__', html.escape(project))
           .replace('__REPORT_CARDS__', _report_cards(available))
           .replace('__PARAM_TILES__', _param_tiles(args, n_in, n_out))
           .replace('__INGROUP_ROWS__', _species_rows(samples, INGROUP_ROLES))
           .replace('__OUTGROUP_ROWS__', _species_rows(samples, OUTGROUP_ROLES))
           .replace('__N_IN__', str(n_in))
           .replace('__N_OUT__', str(n_out)))

    if out.parent != Path(''):
        out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(doc, encoding='utf-8')
    print(f'Wrote {out}: index of {len(available)} reports '
          f'({n_in} ingroup, {n_out} outgroup proteomes)', file=sys.stderr)


if __name__ == '__main__':
    main()
