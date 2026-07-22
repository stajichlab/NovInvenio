"""Publication-quality matplotlib figures for the PDF summary report (issue #19).

Renders a fixed figure set from the same payloads lib/report_data.py builds for the HTML
reports — no new computation. Colours reuse the reports' evidence language so the PDF reads
as the same visual system: series-1 blue = protein-search presence, series-2 green =
TBLASTN genome hit. The categorical palette (blue, orange, purple, green) is CVD- and
print-safe (validated with the dataviz palette checker).

Each fig_* function takes a payload dict (from report_data.build_payload /
build_losses_payload) and returns a matplotlib Figure, so they are unit-testable without
writing a PDF. No I/O at import time; the Agg backend is selected so it works headless.
"""
import matplotlib
matplotlib.use('Agg')  # headless — no display needed
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

# ---- palette (matches the HTML reports; CVD/print-safe categorical) -----------------
BLUE = '#2a78d6'    # series-1: protein-search presence
GREEN = '#008300'   # series-2: TBLASTN hit
ORANGE = '#e69f00'
PURPLE = '#cc79a7'
GRID = '#e1e0d9'
AXIS = '#c3c2b7'
INK = '#0b0b0b'
MUTED = '#898781'
ABSENT = '#ecebe4'  # light cell for "absent" in the heatmap
CATEGORICAL = [BLUE, ORANGE, PURPLE, GREEN]

_STYLE = {
    'figure.facecolor': 'white', 'axes.facecolor': 'white',
    'axes.edgecolor': AXIS, 'axes.linewidth': 0.8,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.titlesize': 12, 'axes.titleweight': 'bold', 'axes.titlecolor': INK,
    'axes.labelcolor': MUTED, 'axes.labelsize': 10,
    'text.color': INK, 'font.size': 10,
    'xtick.color': MUTED, 'ytick.color': MUTED,
    'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'figure.dpi': 150,
}


def _fields(payload):
    return {name: i for i, name in enumerate(payload['fields'])}


def _apply_style():
    plt.rcParams.update(_STYLE)


# --------------------------------------------------------------------------- figures
def fig_summary(payload, losses_payload=None):
    """Title/run-summary page: proteome counts, thresholds, headline candidate counts."""
    _apply_style()
    proteomes = payload['proteomes']
    n_in = sum(1 for p in proteomes if p['group'] == 'IN')
    n_out = len(proteomes) - n_in
    F = _fields(payload)
    n_rows = len(payload['rows'])
    n_nov = sum(1 for r in payload['rows'] if r[F['nov']])
    n_fam = len(payload.get('families', []))

    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.5, 0.90, payload.get('project', 'NovInvenio'), ha='center',
             fontsize=22, fontweight='bold', color=INK)
    fig.text(0.5, 0.865, 'Lineage-specific gene summary', ha='center',
             fontsize=13, color=MUTED)

    lines = [
        ('Ingroup proteomes', str(n_in)),
        ('Outgroup proteomes', str(n_out)),
        ('Proteins scored', f'{n_rows:,}'),
        ('Novelty candidates', f'{n_nov:,}'),
        ('Gene families', f'{n_fam:,}'),
        ('Ingroup presence threshold', f"≥ {payload.get('ingroup_min_frac', 0.75):.0%}"),
    ]
    if losses_payload is not None:
        lines.append(('Loss candidates', f"{len(losses_payload['rows']):,}"))
    y = 0.72
    for label, value in lines:
        fig.text(0.30, y, label, ha='left', fontsize=12, color=MUTED)
        fig.text(0.70, y, value, ha='right', fontsize=12, fontweight='bold', color=INK)
        y -= 0.045
    fig.text(0.5, 0.06, 'Colour key:  blue = protein-search presence   ·   '
             'green = TBLASTN genome hit', ha='center', fontsize=9, color=MUTED)
    return fig


def fig_novelty_per_species(payload):
    """Bar: novelty candidates originating in each ingroup species."""
    _apply_style()
    F = _fields(payload)
    proteomes = payload['proteomes']
    ingroup_idx = [i for i, p in enumerate(proteomes) if p['group'] == 'IN']
    counts = {i: 0 for i in ingroup_idx}
    for r in payload['rows']:
        if r[F['nov']] and r[F['src']] in counts:
            counts[r[F['src']]] += 1
    labels = [proteomes[i]['short'] for i in ingroup_idx]
    values = [counts[i] for i in ingroup_idx]
    order = sorted(range(len(values)), key=lambda k: values[k], reverse=True)
    labels = [labels[k] for k in order]
    values = [values[k] for k in order]

    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    bars = ax.bar(range(len(labels)), values, color=BLUE, width=0.72)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel('novelty candidates')
    ax.set_title('Novelty candidates per ingroup species')
    ax.grid(axis='x', visible=False)
    for b, v in zip(bars, values):
        if v:
            ax.text(b.get_x() + b.get_width() / 2, v, str(v), ha='center', va='bottom',
                    fontsize=7, color=MUTED)
    fig.tight_layout()
    return fig


def _top_novelty_rows(payload, top_n):
    F = _fields(payload)
    novel = [r for r in payload['rows'] if r[F['nov']]]
    novel.sort(key=lambda r: r[F['pres']].count('1'), reverse=True)
    return novel[:top_n], F


def fig_presence_heatmap(payload, top_n=40):
    """Static twin of the interactive heatmap: top-N novelty candidates × proteomes.

    Cells are present (blue) / absent (light); a TBLASTN genome hit is marked green.
    """
    _apply_style()
    rows, F = _top_novelty_rows(payload, top_n)
    shorts = [p['short'] for p in payload['proteomes']]
    tb_genomes = payload.get('tblastn_genomes', [])
    n_in = sum(1 for p in payload['proteomes'] if p['group'] == 'IN')

    if not rows:
        fig, ax = plt.subplots(figsize=(8.5, 4))
        ax.text(0.5, 0.5, 'No novelty candidates', ha='center', va='center', color=MUTED)
        ax.axis('off')
        return fig

    mat = [[1 if r[F['pres']][j] == '1' else 0 for j in range(len(shorts))] for r in rows]
    fig, ax = plt.subplots(figsize=(8.5, max(3.0, 0.22 * len(rows) + 1.5)))
    ax.imshow(mat, aspect='auto', cmap=ListedColormap([ABSENT, BLUE]),
              interpolation='nearest', vmin=0, vmax=1)
    # TBLASTN hits (green dots) over the outgroup columns they map to.
    for ri, r in enumerate(rows):
        tb = r[F['tb']]
        for gi, g in enumerate(tb_genomes):
            if gi < len(tb) and tb[gi] == '1' and g in shorts:
                ax.plot(shorts.index(g), ri, 'o', color=GREEN, markersize=4)
    ax.set_xticks(range(len(shorts)))
    ax.set_xticklabels(shorts, rotation=90, fontsize=6)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[F['id']] for r in rows], fontsize=6)
    ax.axvline(n_in - 0.5, color=AXIS, linewidth=1.2)  # ingroup | outgroup divider
    ax.set_title(f'Presence of the top {len(rows)} novelty candidates '
                 '(blue = present · green = TBLASTN hit)')
    ax.grid(False)
    fig.tight_layout()
    return fig


def fig_family_sizes(payload):
    """Histogram of gene-family sizes."""
    _apply_style()
    sizes = [fam.get('size', 0) for fam in payload.get('families', [])]
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    if sizes:
        ax.hist(sizes, bins=range(2, max(sizes) + 2), color=BLUE, edgecolor='white',
                linewidth=0.5, align='left')
    ax.set_xlabel('members per family')
    ax.set_ylabel('families')
    ax.set_title('Gene-family size distribution')
    ax.grid(axis='x', visible=False)
    fig.tight_layout()
    return fig


def fig_pfam_frequency(payload, top_n=15):
    """Horizontal bar: most frequent Pfam domains among novelty candidates."""
    _apply_style()
    F = _fields(payload)
    counts = {}
    for r in payload['rows']:
        if not r[F['nov']]:
            continue
        for dom in (r[F['pfam_n']] or '').split(','):
            dom = dom.strip()
            if dom:
                counts[dom] = counts.get(dom, 0) + 1
    top = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)[:top_n]
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    if top:
        labels = [k for k, _ in top][::-1]
        values = [v for _, v in top][::-1]
        ax.barh(range(len(labels)), values, color=ORANGE, height=0.72)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        for i, v in enumerate(values):
            ax.text(v, i, f' {v}', va='center', fontsize=7, color=MUTED)
    else:
        ax.text(0.5, 0.5, 'No Pfam annotation on novelty candidates', ha='center',
                va='center', color=MUTED, transform=ax.transAxes)
    ax.set_xlabel('novelty candidates with the domain')
    ax.set_title('Top Pfam domains among novelty candidates')
    ax.grid(axis='y', visible=False)
    fig.tight_layout()
    return fig


def fig_annotation_source(payload):
    """Bar: annotation source of the novelty candidates (model-org / Pfam / SwissProt / none)."""
    _apply_style()
    F = _fields(payload)
    fsources = payload.get('fsources', [])
    counts = {}
    for r in payload['rows']:
        if not r[F['nov']]:
            continue
        idx = r[F['fsrc']]
        key = fsources[idx] if 0 <= idx < len(fsources) else 'none'
        counts[key] = counts.get(key, 0) + 1
    items = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    if items:
        labels = [k for k, _ in items]
        values = [v for _, v in items]
        colors = [MUTED if lab.lower() == 'none' else CATEGORICAL[i % len(CATEGORICAL)]
                  for i, lab in enumerate(labels)]
        bars = ax.bar(range(len(labels)), values, color=colors, width=0.66)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=8)
        for b, v in zip(bars, values):
            ax.text(b.get_x() + b.get_width() / 2, v, str(v), ha='center', va='bottom',
                    fontsize=7, color=MUTED)
    ax.set_ylabel('novelty candidates')
    ax.set_title('Annotation source of novelty candidates')
    ax.grid(axis='x', visible=False)
    fig.tight_layout()
    return fig


def fig_loss_scatter(losses_payload):
    """Scatter of candidate losses: outgroup breadth vs ingroup species still retaining."""
    _apply_style()
    F = _fields(losses_payload)
    xs = [r[F['out_breadth']] for r in losses_payload['rows']]
    ys = [r[F['in_retained']] for r in losses_payload['rows']]
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    if xs:
        # clean losses (in_retained == 0) in green, the rest muted.
        clean = [(x, y) for x, y in zip(xs, ys) if y == 0]
        rest = [(x, y) for x, y in zip(xs, ys) if y > 0]
        if rest:
            ax.scatter([x for x, _ in rest], [y for _, y in rest], s=22, color=MUTED,
                       alpha=0.6, edgecolor='white', linewidth=0.4, label='retained in some ingroup')
        if clean:
            ax.scatter([x for x, _ in clean], [y for _, y in clean], s=26, color=GREEN,
                       edgecolor='white', linewidth=0.4, label='clean loss (0 ingroup)')
        ax.legend(frameon=False, fontsize=8)
    else:
        ax.text(0.5, 0.5, 'No loss candidates', ha='center', va='center', color=MUTED,
                transform=ax.transAxes)
    ax.set_xlabel('outgroup species carrying the family (breadth)')
    ax.set_ylabel('ingroup species still retaining a member')
    ax.set_title('Candidate gene losses')
    fig.tight_layout()
    return fig


def build_all(payload, losses_payload=None):
    """The ordered figure set for the PDF (novelty payload; losses optional)."""
    figs = [
        fig_summary(payload, losses_payload),
        fig_novelty_per_species(payload),
        fig_presence_heatmap(payload),
        fig_family_sizes(payload),
        fig_pfam_frequency(payload),
        fig_annotation_source(payload),
    ]
    if losses_payload is not None:
        figs.append(fig_loss_scatter(losses_payload))
    return figs
