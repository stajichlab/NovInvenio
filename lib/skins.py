"""
Colour-skin registry for every NovInvenio HTML report.

This is the single source of truth for the CSS custom properties the report
pages paint with. Before this module the same token block was written out four
times -- inline in lib/report_template.py, again as report_common.THEME_VARS_CSS,
and twice more as raw hex in view/generate_index.py -- so a fifth palette meant
editing four files and remembering lib/losses_report_template.py's hardcoded
``.badge.warn`` escape hatch. Everything now reads ``skins_css()``.

Selection model -- three states, matching the artifact/theme convention:

  * no ``data-skin`` attribute on <html>  -> follow the OS (paper / dark)
  * ``data-skin="<id>"``                  -> that skin, overriding the OS
  * ``@media print``                      -> always the paper tokens

``skin_boot_js()`` belongs in <head> so a stored choice is applied before first
paint; ``skin_picker_js()`` goes with the rest of the page script and calls the
page's ``window.onSkinChange`` hook (lib/report_template.py uses it to repaint
the canvas heatmap, which reads these same tokens through getComputedStyle).

Adding a skin
-------------
Add an entry to SKINS. Every skin must define *every* token in REQUIRED_TOKENS
-- tests/test_skins.py enforces that, and also enforces the contrast floors
below, so a new palette cannot quietly regress legibility.

Colour-role contract (see CLAUDE.md's report constraints): ``--series-1`` marks
presence from the protein search and ``--series-2`` marks a TBLASTN genome hit;
``--on-series`` is the label colour drawn on top of either fill (white over the
mid-tone blues and greens, near-black over the neon skins, where white would be
unreadable).
Hue carries evidence *type*; ingroup vs outgroup is carried by column position.
Both series must stay separable under red-green colour-vision deficiency, which
is why the neon skin pairs cyan with amber (separated on the blue-yellow axis,
which every common form of CVD preserves) rather than the obvious cyberpunk
cyan/magenta -- magenta desaturates toward blue-grey under deuteranopia and
converges with cyan.
"""
from __future__ import annotations

# WCAG 2.1 floors asserted by tests/test_skins.py.
MIN_TEXT_CONTRAST = 4.5      # body text on the page background (AA, normal text)
MIN_SERIES_CONTRAST = 3.0    # data marks / UI graphics on the page background (AA)
# Label text drawn *on top of* a --series-* fill (the presence chips). The floor
# is 3.0 rather than 4.5 because the shipped dark palette's white-on-#3987e5
# already sits at 3.64; raising it would mean recolouring a validated palette.
MIN_ON_SERIES_CONTRAST = 3.0

_SYSTEM_UI = (
    'system-ui, -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif'
)
# Self-contained pages can never fetch a webfont (CLAUDE.md: the reports open
# from file:// with no network), so every stack here is system fonts only.
# DejaVu Sans Mono covers the Linux hosts these reports get opened on.
_MONO = (
    'ui-monospace, SFMono-Regular, Menlo, "DejaVu Sans Mono", Consolas, monospace'
)

# Static scanlines for the neon skin. Deliberately not animated: animated
# texture over a 20k-row canvas heatmap is both a legibility and an
# accessibility problem. This is painted as a background layer on `body`
# (BASE_PAGE_CSS), so it sits in the page's negative space and can never land
# on a card, table, tooltip or detail panel -- an element's background always
# paints below its descendants. Printing pins the paper tokens, where this is
# `none`, so it never reaches paper either.
# The alpha can afford to be this visible precisely because the texture is
# confined to the page ground: it never sits behind text or heatmap cells, so
# it costs no legibility anywhere.
_SCANLINES = (
    'repeating-linear-gradient(to bottom, '
    'rgba(0, 229, 255, 0.06) 0px, rgba(0, 229, 255, 0.06) 1px, '
    'transparent 1px, transparent 3px)'
)

SKINS: dict[str, dict] = {
    # The two original palettes, carried over token-for-token so the refactor
    # that introduced this module changed no pixels. Both were already
    # validated for colour-blind separation and contrast -- re-validate before
    # editing (CLAUDE.md).
    'paper': {
        'label': 'Paper',
        'scheme': 'light',
        'note': 'The default. Neutral light palette for figures and screenshots.',
        'tokens': {
            '--page': '#f9f9f7',
            '--surface-1': '#fcfcfb',
            '--text-primary': '#0b0b0b',
            '--text-secondary': '#52514e',
            '--muted': '#898781',
            '--grid': '#e1e0d9',
            '--axis': '#c3c2b7',
            '--border': 'rgba(11, 11, 11, 0.10)',
            '--series-1': '#2a78d6',
            '--series-2': '#008300',
            '--on-series': '#ffffff',
            '--wash': 'rgba(42, 120, 214, 0.10)',
            '--hover-wash': 'rgba(11, 11, 11, 0.05)',
            '--warn': '#a15c00',
            '--shadow': '0 6px 24px rgba(11, 11, 11, 0.18)',
            '--font-ui': _SYSTEM_UI,
            '--font-mono': _MONO,
            '--glow': 'none',
            '--overlay': 'none',
        },
    },
    'dark': {
        'label': 'Dark',
        'scheme': 'dark',
        'note': 'The default when the OS is in dark mode.',
        'tokens': {
            '--page': '#0d0d0d',
            '--surface-1': '#1a1a19',
            '--text-primary': '#ffffff',
            '--text-secondary': '#c3c2b7',
            '--muted': '#898781',
            '--grid': '#2c2c2a',
            '--axis': '#383835',
            '--border': 'rgba(255, 255, 255, 0.10)',
            '--series-1': '#3987e5',
            '--series-2': '#008300',
            '--on-series': '#ffffff',
            '--wash': 'rgba(57, 135, 229, 0.16)',
            '--hover-wash': 'rgba(255, 255, 255, 0.06)',
            '--warn': '#e0a030',
            '--shadow': '0 6px 24px rgba(0, 0, 0, 0.55)',
            '--font-ui': _SYSTEM_UI,
            '--font-mono': _MONO,
            '--glow': 'none',
            '--overlay': 'none',
        },
    },
    # Near-black rather than #000: pure black behind saturated neon is the
    # single biggest source of eye strain in this style. Neon is reserved for
    # data marks and accents -- body text stays a desaturated blue-grey.
    'neuromancer': {
        'label': 'Neuromancer',
        'scheme': 'dark',
        'note': 'Terminal-phosphor neon. Cyan presence, amber genome hits.',
        'tokens': {
            '--page': '#0a0e14',
            '--surface-1': '#111820',
            '--text-primary': '#c6d3d6',
            '--text-secondary': '#8fa3a7',
            '--muted': '#6b8085',
            '--grid': '#1b2733',
            '--axis': '#2b3b48',
            '--border': 'rgba(0, 229, 255, 0.16)',
            '--series-1': '#00e5ff',
            '--series-2': '#ff9d00',
            '--on-series': '#0a0e14',
            '--wash': 'rgba(0, 229, 255, 0.14)',
            '--hover-wash': 'rgba(0, 229, 255, 0.07)',
            '--warn': '#ff5f56',
            '--shadow': '0 6px 24px rgba(0, 0, 0, 0.65)',
            '--font-ui': _MONO,
            '--font-mono': _MONO,
            '--glow': '0 0 6px rgba(0, 229, 255, 0.35)',
            '--overlay': _SCANLINES,
        },
    },
    'contrast': {
        'label': 'High contrast',
        'scheme': 'dark',
        'note': 'Maximum separation for projectors and low-vision reading.',
        'tokens': {
            '--page': '#000000',
            '--surface-1': '#0a0a0a',
            '--text-primary': '#ffffff',
            '--text-secondary': '#e6e6e6',
            '--muted': '#b8b8b8',
            '--grid': '#3a3a3a',
            '--axis': '#6a6a6a',
            '--border': 'rgba(255, 255, 255, 0.35)',
            '--series-1': '#4db8ff',
            '--series-2': '#ffd400',
            '--on-series': '#000000',
            '--wash': 'rgba(77, 184, 255, 0.22)',
            '--hover-wash': 'rgba(255, 255, 255, 0.12)',
            '--warn': '#ff6b6b',
            '--shadow': '0 0 0 1px rgba(255, 255, 255, 0.45), 0 6px 24px rgba(0, 0, 0, 0.8)',
            '--font-ui': _SYSTEM_UI,
            '--font-mono': _MONO,
            '--glow': 'none',
            '--overlay': 'none',
        },
    },
}

# The skin used when the OS reports light mode / no preference, and the one
# @media print forces so a printed report is legible whatever is on screen.
DEFAULT_SKIN = 'paper'
# The skin used when the OS reports dark mode and no explicit choice is stored.
DEFAULT_DARK_SKIN = 'dark'

# Every skin must define exactly this set -- a missing token would silently
# inherit from whichever skin was applied before it.
REQUIRED_TOKENS = frozenset(SKINS[DEFAULT_SKIN]['tokens'])

# localStorage key. Shared across novelties.html / core.html / losses.html /
# report.html so a choice carries between the reports of one project (same
# origin on GitHub Pages; browsers that give each file:// document an opaque
# origin simply won't persist it, which degrades to "follow system").
STORAGE_KEY = 'novinvenio.skin'


def _block(tokens: dict[str, str], scheme: str, indent: str = '      ') -> str:
    lines = [f'{indent}color-scheme: {scheme};']
    lines += [f'{indent}{name}: {value};' for name, value in tokens.items()]
    return '\n'.join(lines)


def skins_css() -> str:
    """Return the full ``:root`` custom-property CSS for every skin.

    Emits, in order: the default (paper) tokens on bare ``:root``; a dark
    override for OS dark mode that applies only when no explicit skin is set;
    one ``[data-skin="<id>"]`` block per skin so an explicit choice wins in
    both directions; and a print override pinning the paper tokens.
    """
    parts = [
        '  :root {',
        _block(SKINS[DEFAULT_SKIN]['tokens'], SKINS[DEFAULT_SKIN]['scheme'], '    '),
        '  }',
        '  @media (prefers-color-scheme: dark) {',
        '    :root:where(:not([data-skin])) {',
        _block(
            SKINS[DEFAULT_DARK_SKIN]['tokens'],
            SKINS[DEFAULT_DARK_SKIN]['scheme'],
        ),
        '    }',
        '  }',
    ]
    for skin_id, skin in SKINS.items():
        parts += [
            f'  :root[data-skin="{skin_id}"] {{',
            _block(skin['tokens'], skin['scheme'], '    '),
            '  }',
        ]
    # Printing a neon skin wastes toner and is unreadable; a printed report is
    # a figure, so it always prints as paper.
    parts += [
        '  @media print {',
        '    :root, :root[data-skin] {',
        _block(SKINS[DEFAULT_SKIN]['tokens'], SKINS[DEFAULT_SKIN]['scheme']),
        '    }',
        '  }',
    ]
    return '\n'.join(parts) + '\n'


def skin_picker_html(indent: str = '    ') -> str:
    """Return the skin ``<select>`` markup for a report's header."""
    opts = [f'{indent}  <option value="">Follow system</option>']
    opts += [
        f'{indent}  <option value="{sid}">{s["label"]}</option>'
        for sid, s in SKINS.items()
    ]
    return (
        f'{indent}<label class="sr-only" for="skin">Colour skin</label>\n'
        f'{indent}<select id="skin" title="Colour skin">\n'
        + '\n'.join(opts)
        + f'\n{indent}</select>'
    )


def skin_boot_js() -> str:
    """Return the <head> snippet that applies a stored skin before first paint.

    Runs synchronously in <head> precisely so the page never flashes the
    default palette before the stored choice lands. Kept independent of
    ``skin_picker_js()`` (which runs at the end of <body>, where the <select>
    exists) for that reason.
    """
    valid = ', '.join(f'"{s}"' for s in SKINS)
    return f"""
(function () {{
  try {{
    var v = localStorage.getItem("{STORAGE_KEY}");
    if (v && [{valid}].indexOf(v) >= 0) {{
      document.documentElement.setAttribute("data-skin", v);
    }}
  }} catch (e) {{ /* opaque file:// origin, or storage disabled */ }}
}})();
"""


def skin_picker_js() -> str:
    """Return the skin ``<select>`` wiring.

    Calls ``window.onSkinChange`` (if a page defines one) after any change to
    the effective palette, including an OS light/dark flip while "Follow
    system" is selected -- lib/report_template.py hooks that to repaint its
    canvas heatmap, which reads the tokens via getComputedStyle.
    """
    return f"""
  var skinSel = document.getElementById("skin");
  function notifySkinChange() {{
    if (typeof window.onSkinChange === "function") window.onSkinChange();
  }}
  function storeSkin(v) {{
    try {{
      if (v) localStorage.setItem("{STORAGE_KEY}", v);
      else localStorage.removeItem("{STORAGE_KEY}");
    }} catch (e) {{ /* storage unavailable -- selection still applies for this page */ }}
  }}
  skinSel.value = document.documentElement.getAttribute("data-skin") || "";
  skinSel.addEventListener("change", function () {{
    var v = skinSel.value;
    if (v) document.documentElement.setAttribute("data-skin", v);
    else document.documentElement.removeAttribute("data-skin");
    storeSkin(v);
    notifySkinChange();
  }});
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", notifySkinChange);
"""


# ---- contrast helpers (used by tests/test_skins.py) ------------------------

def _srgb_to_linear(c: float) -> float:
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def relative_luminance(hex_colour: str) -> float:
    """WCAG relative luminance of an opaque ``#rrggbb`` colour."""
    h = hex_colour.strip().lstrip('#')
    if len(h) == 3:
        h = ''.join(ch * 2 for ch in h)
    if len(h) != 6:
        raise ValueError(f'not an opaque hex colour: {hex_colour!r}')
    r, g, b = (int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))
    return (
        0.2126 * _srgb_to_linear(r)
        + 0.7152 * _srgb_to_linear(g)
        + 0.0722 * _srgb_to_linear(b)
    )


def contrast_ratio(fg: str, bg: str) -> float:
    """WCAG contrast ratio between two opaque ``#rrggbb`` colours."""
    a, b = relative_luminance(fg), relative_luminance(bg)
    lo, hi = sorted((a, b))
    return (hi + 0.05) / (lo + 0.05)
