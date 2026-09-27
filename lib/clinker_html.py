"""Slim clinker 0.0.32 HTML pages (spec section 8).

clinker embeds its data as `const data={...};` in one inline <script> and
copies both genes' DNA and protein sequence into every link. clinker's own
page scripts never read the `sequence` or `translation` fields (checked in
clinker 0.0.32's plot/clinker.js and clustermap.min.js), so removing them
only removes bytes: 4.13 MB -> 0.81 MB on the 2026-09-26 spike page and
5.35 MB -> 1.01 MB on a real 12-strain locus, with the same drawn figure
(tests/test_clinker_render.py). No I/O at import time.
"""
from __future__ import annotations

import json

MARKER = "const data="
STRIP = ("sequence", "translation")

# C2/C3 (spec section 8): verified 2026-09-27 against a real clinker 0.0.32
# page (spike run under $SCRATCH/clinker-venv) -- see
# tests/test_clinker_html.py's REAL_CLINKER_PAGE for the exact markup this
# matches. clinker's #div-floater options sidebar starts open and, in a
# panel-sized iframe, covers most of the figure (C2); its body sets
# `overflow: hidden`, which clips a figure taller than the iframe instead
# of scrolling (C3), and the browser's default body margin can clip the
# left-most locus labels in a narrow iframe.
_SIDEBAR_ACTIVE = 'class="collapsible active"'
_SIDEBAR_COLLAPSED = 'class="collapsible"'
_SUMMARY_OPEN = '<div id="div-summary">'
_SUMMARY_COLLAPSED = '<div id="div-summary" style="display:none">'
_HEAD_CLOSE = "</head>"
_UI_STYLE = "<style>body{overflow:auto;margin-left:16px}</style>"

# L (final review, verified fix): clinker's own default plot() config aligns
# every locus label to the left edge of its drawn genes; combined with a
# long cluster/locus name, the label can be drawn at a negative x and
# clipped off the left of the figure. Turning on clinker's own "hide locus
# coordinates" option (config + its sidebar checkbox) removes one source of
# label width, and the small window-load script below pans the whole
# drawing right if the label group's bounding box still starts left of the
# SVG's own left edge. Verified 2026-09-27 against a real clinker 0.0.32
# page (markup below matches tests/test_clinker_html.py's
# REAL_PAGE_WITH_LABEL_MARKUP) and, per the final-review brief, against a
# real generated page under headless Chromium (`--dump-dom`): cluster
# label boxes at x >= 0.
_ALIGN_LABELS_CONFIG = "alignLabels: true,"
_HIDE_COORDS_CONFIG = "alignLabels: true,\n        hideLocusCoordinates: true,"
_HIDE_COORDS_CHECKBOX = '<input id="input-cluster-hide-coords" type="checkbox">'
_HIDE_COORDS_CHECKBOX_CHECKED = '<input id="input-cluster-hide-coords" type="checkbox" checked>'
_BODY_CLOSE = "</body>"
_LABEL_FIT_SCRIPT = (
    '<script>(function(){function fit(){var svg=document.querySelector("svg.clusterMap"),'
    'g=svg&&svg.querySelector("g.clusterMapG");if(!g||!window.d3)return;var t=d3.zoomTransform(svg),'
    'b=g.getBBox(),left=t.x+b.x*t.k,pad=8;if(left>=pad)return;var nt=d3.zoomIdentity.translate('
    't.x+pad-left,t.y).scale(t.k);svg.__zoom=nt;g.setAttribute("transform",nt.toString());}'
    'window.addEventListener("load",function(){setTimeout(fit,600);});})();</script>'
)


def inject_ui_fixes(html: str) -> str:
    """clinker's page with the options sidebar starting collapsed (its own
    unmodified toggleActive() still opens it: it reads the inline
    `display:none` this sets, so the first click reopens it correctly),
    scrolling/a left margin restored, locus coordinates hidden by default
    (L), and a window-load script that pans the drawing right if a label
    is still clipped off the left edge (L). A no-op for any marker not
    found, so a clinker version whose markup differs never crashes the
    run. Idempotent: calling it again on its own output is a no-op."""
    out = html.replace(_SIDEBAR_ACTIVE, _SIDEBAR_COLLAPSED, 1)
    out = out.replace(_SUMMARY_OPEN, _SUMMARY_COLLAPSED, 1)
    # Fix round 1, item 2: idempotent -- a second call (e.g. a retried
    # pipeline step given an already-fixed page) must not add a second
    # <style> tag.
    if _HEAD_CLOSE in out and _UI_STYLE not in out:
        out = out.replace(_HEAD_CLOSE, _UI_STYLE + _HEAD_CLOSE, 1)
    if _HIDE_COORDS_CONFIG not in out and _ALIGN_LABELS_CONFIG in out:
        out = out.replace(_ALIGN_LABELS_CONFIG, _HIDE_COORDS_CONFIG, 1)
    if _HIDE_COORDS_CHECKBOX in out:
        out = out.replace(_HIDE_COORDS_CHECKBOX, _HIDE_COORDS_CHECKBOX_CHECKED, 1)
    if _LABEL_FIT_SCRIPT not in out and _BODY_CLOSE in out:
        out = out.replace(_BODY_CLOSE, _LABEL_FIT_SCRIPT + _BODY_CLOSE, 1)
    return out


def _strip_gene(gene: dict) -> None:
    for key in STRIP:
        gene.pop(key, None)


def slim_clinker_html(html: str) -> str:
    """The page with sequence/translation removed from every gene in
    data.clusters[].loci[].genes[] and data.links[].query/target, and the
    UI fixes above applied. `</` in the re-serialised data is escaped so a
    label cannot close the script. Raises ValueError when the page has no
    `const data=` object."""
    start = html.find(MARKER)
    if start < 0:
        raise ValueError("no 'const data=' block: not a clinker HTML page")
    begin = start + len(MARKER)
    data, end = json.JSONDecoder().raw_decode(html, begin)
    for cluster in data.get("clusters", []):
        for locus in cluster.get("loci", []):
            for gene in locus.get("genes", []):
                _strip_gene(gene)
    for link in data.get("links", []):
        for side in ("query", "target"):
            if isinstance(link.get(side), dict):
                _strip_gene(link[side])
    blob = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    return inject_ui_fixes(html[:begin] + blob + html[end:])
