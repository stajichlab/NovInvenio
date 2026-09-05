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
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import parse_config  # noqa: E402
from index_page import render_project_page  # noqa: E402
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
    return [
        {
            'short': s.short,
            'species': s.species,
            'strain': s.strain,
            'taxon': s.taxon_group,
        }
        for s in samples if s.group in roles
    ]


def _param_tiles(args, n_in, n_out):
    tiles = [
        ('Ingroup proteomes', str(n_in)),
        ('Outgroup proteomes', str(n_out)),
    ]
    if args.run_tool:
        tiles.append(('Search tool', args.run_tool))
    for label, value in (
        ('Ingroup min fraction', args.ingroup_min_frac),
        ('Outgroup min fraction', args.outgroup_min_frac),
        ('Loss ingroup max fraction', args.loss_ingroup_max_frac),
        ('Core min fraction', args.core_min_frac),
    ):
        if value is not None:
            tiles.append((label, f'{value:g}'))
    return tiles


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

    doc = render_project_page(
        project=project,
        reports=[
            {'href': fname, 'title': title, 'desc': desc}
            for fname, title, desc in REPORTS if fname in available
        ],
        tiles=_param_tiles(args, n_in, n_out),
        ingroup=_species_rows(samples, INGROUP_ROLES),
        outgroup=_species_rows(samples, OUTGROUP_ROLES),
        note=f'Parameters and proteomes from {Path(args.config).name}.',
    )

    if out.parent != Path(''):
        out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(doc, encoding='utf-8')
    print(f'Wrote {out}: index of {len(available)} reports '
          f'({n_in} ingroup, {n_out} outgroup proteomes)', file=sys.stderr)


if __name__ == '__main__':
    main()
