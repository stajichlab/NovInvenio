"""
Landing pages for a NovInvenio result set: the per-project ``report.html`` and
the top-level ``view/index.html`` gallery.

Both are rendered here so there is exactly one design for each. Before this,
``bin/make_index_report.py`` (pipeline, driven by the config CSV) and
``view/generate_index.py`` (post-hoc, driven by the payloads embedded in the
report HTML) each rendered their *own* report.html with different markup and
different CSS -- so whichever ran last decided what the shared folder looked
like, and the top-level gallery matched neither. The two callers still differ
in where their facts come from; they no longer differ in what they draw.

Same constraints as the reports themselves: self-contained, opens from
file://, no network. Callers pass already-escaped-free plain text -- escaping
happens here.
"""
from __future__ import annotations

import html

from report_common import (
    BASE_PAGE_CSS,
    SKIN_BOOT_JS,
    SKIN_PICKER_HTML,
    SKIN_PICKER_JS,
    SKIN_VARS_CSS,
)

# What NovInvenio is, for a reader who arrived at the gallery from a link and
# has no idea what any of this is. The gallery is a public GitHub Pages front
# door, so it cannot assume the pipeline's own vocabulary.
BLURB = (
    'NovInvenio identifies lineage-specific genes — proteins conserved across a '
    'defined ingroup but absent from every outgroup proteome — by combining '
    'pairwise protein search, gene-family clustering, TBLASTN validation against '
    'outgroup genomes, and Pfam/SwissProt annotation. Each project below is one '
    'analysis run; open its summary for the species and thresholds used.'
)

_EXTRA_CSS = r"""
  .cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 14px; }
  a.card-link {
    display: flex; flex-direction: column; gap: 6px; text-decoration: none;
    background: var(--surface-1); border: 1px solid var(--border); border-radius: 10px;
    padding: 16px; color: var(--text-primary);
  }
  a.card-link:hover { border-color: var(--series-1); background: var(--hover-wash); }
  .card-link-title { font-size: 15px; font-weight: 600; }
  .card-link-desc { font-size: 12px; color: var(--text-secondary); }
  .card-link-meta { font-size: 11px; color: var(--muted); }
  .groups { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
  @media (max-width: 800px) { .groups { grid-template-columns: 1fr; } }
  .subpages { display: flex; flex-wrap: wrap; gap: 6px; margin-top: 4px; }
  .subpages a {
    font-size: 11px; padding: 3px 8px; border: 1px solid var(--border); border-radius: 999px;
    text-decoration: none; color: var(--text-secondary); background: var(--page);
  }
  .subpages a:hover { border-color: var(--series-1); color: var(--text-primary); }
  .intro { max-width: 70ch; color: var(--text-secondary); font-size: 13px; }
  footer { margin-top: 32px; color: var(--muted); font-size: 12px; }
"""


def _page(title: str, body: str) -> str:
    return (
        '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f'<title>{html.escape(title)}</title>\n'
        '<style>\n' + SKIN_VARS_CSS + BASE_PAGE_CSS + _EXTRA_CSS + '</style>\n'
        '<script>' + SKIN_BOOT_JS + '</script>\n'
        '</head>\n<body>\n<div class="wrap">\n'
        + body +
        '\n</div>\n<script>\n(function () {\n  "use strict";\n'
        + SKIN_PICKER_JS +
        '})();\n</script>\n</body>\n</html>\n'
    )


def _header(title: str, subtitle: str) -> str:
    return (
        '  <header class="top">\n'
        '    <div class="titles">\n'
        f'      <h1>{html.escape(title)}</h1>\n'
        f'      <p class="sub">{html.escape(subtitle)}</p>\n'
        '    </div>\n'
        + SKIN_PICKER_HTML + '\n'
        '  </header>\n'
    )


def _tiles(tiles: list[tuple[str, str]]) -> str:
    if not tiles:
        return '<p class="card-note">No run parameters recorded.</p>'
    return '\n'.join(
        f'      <div><div class="tile-label">{html.escape(label)}</div>'
        f'<div class="tile-value">{html.escape(str(value))}</div></div>'
        for label, value in tiles
    )


def _species_table(rows: list[dict]) -> str:
    if not rows:
        return '<p class="card-note">None recorded.</p>'
    cells = []
    for r in rows:
        name = html.escape(r.get('species', ''))
        if r.get('strain'):
            name += ' ' + html.escape(r['strain'])
        cells.append(
            f'<tr><td class="mono">{html.escape(r.get("short", ""))}</td>'
            f'<td>{name}</td><td>{html.escape(r.get("taxon", "") or "")}</td></tr>'
        )
    return (
        '<table class="data">\n'
        '<thead><tr><th scope="col">Short</th><th scope="col">Species</th>'
        '<th scope="col">Taxon group</th></tr></thead>\n<tbody>\n'
        + '\n'.join(cells) +
        '\n</tbody></table>'
    )


def _report_cards(reports: list[dict]) -> str:
    if not reports:
        return '<p class="card-note">No reports found in this folder.</p>'
    out = []
    for r in reports:
        meta = (
            f'<span class="card-link-meta">{html.escape(r["meta"])}</span>'
            if r.get('meta') else ''
        )
        out.append(
            f'      <a class="card-link" href="{html.escape(r["href"])}">'
            f'<span class="card-link-title">{html.escape(r["title"])}</span>'
            f'<span class="card-link-desc">{html.escape(r.get("desc", ""))}</span>'
            f'{meta}</a>'
        )
    return '\n'.join(out)


def render_project_page(
    *,
    project: str,
    reports: list[dict],
    tiles: list[tuple[str, str]],
    ingroup: list[dict],
    outgroup: list[dict],
    note: str = '',
    footer: str = '',
) -> str:
    """Render one project folder's ``report.html``.

    ``reports`` entries are ``{href, title, desc, meta}``; ``ingroup`` and
    ``outgroup`` entries are ``{short, species, strain, taxon}``. ``note``
    says where the summary came from (a config CSV, or the report payloads),
    since the two callers reconstruct it differently.
    """
    body = (
        _header(f'{project} — NovInvenio reports',
                'Lineage-specific gene analysis · one folder, three interactive reports')
        + '\n  <section class="card">\n'
          '    <h2 class="card-title">Reports</h2>\n'
          '    <p class="card-note">Each page is self-contained — open it directly, '
          'no network needed.</p>\n'
          '    <div class="cards">\n' + _report_cards(reports) + '\n    </div>\n'
          '  </section>\n'
          '\n  <section class="card">\n'
          '    <h2 class="card-title">Run summary</h2>\n'
        + (f'    <p class="card-note">{html.escape(note)}</p>\n' if note else '')
        + '    <div class="tiles">\n' + _tiles(tiles) + '\n    </div>\n'
          '  </section>\n'
          '\n  <div class="groups">\n'
          '    <section class="card">\n'
          f'      <h2 class="card-title">Ingroup ({len(ingroup)})</h2>\n'
          '      <p class="card-note">Genes are called lineage-specific relative to '
          'this clade.</p>\n'
        + _species_table(ingroup) + '\n    </section>\n'
          '    <section class="card">\n'
          f'      <h2 class="card-title">Outgroup ({len(outgroup)})</h2>\n'
          '      <p class="card-note">Background proteomes used to exclude broadly '
          'conserved genes.</p>\n'
        + _species_table(outgroup) + '\n    </section>\n'
          '  </div>\n'
        + (f'  <footer>{html.escape(footer)}</footer>\n' if footer else '')
    )
    return _page(f'{project} — NovInvenio reports', body)


def render_gallery_page(*, projects: list[dict], footer: str = '') -> str:
    """Render the top-level ``view/index.html``.

    ``projects`` entries are ``{title, href, meta, tiles, subpages}`` where
    ``subpages`` is a list of ``{href, title}``. This is the GitHub Pages front
    door, so it leads with what the project *is* rather than a bare list of
    directory names -- the counts shown per card come from the same report
    payloads the per-project summaries are built from.
    """
    cards = []
    for p in projects:
        subs = ''
        if p.get('subpages'):
            links = '\n'.join(
                f'        <a href="{html.escape(sp["href"])}">{html.escape(sp["title"])}</a>'
                for sp in p['subpages']
            )
            subs = f'\n      <div class="subpages">\n{links}\n      </div>'
        meta = (
            f'<span class="card-link-meta">{html.escape(p["meta"])}</span>'
            if p.get('meta') else ''
        )
        cards.append(
            f'      <div class="card-wrap">\n'
            f'      <a class="card-link" href="{html.escape(p["href"])}">'
            f'<span class="card-link-title">{html.escape(p["title"])}</span>'
            f'<span class="card-link-desc">{html.escape(p.get("desc", ""))}</span>'
            f'{meta}</a>{subs}\n      </div>'
        )

    body = (
        _header('NovInvenio reports', 'Lineage-specific gene discovery across fungal genomes')
        + '\n  <section class="card">\n'
          f'    <p class="intro">{html.escape(BLURB)}</p>\n'
          '  </section>\n'
          '\n  <section class="card">\n'
          '    <h2 class="card-title">Projects</h2>\n'
          '    <p class="card-note">One analysis run each. Counts come from the '
          'reports themselves.</p>\n'
          '    <div class="cards">\n'
        + ('\n'.join(cards) if cards else
           '      <p class="card-note">No project folders found.</p>')
        + '\n    </div>\n  </section>\n'
        + (f'  <footer>{html.escape(footer)}</footer>\n' if footer else '')
    )
    return _page('NovInvenio reports', body)
