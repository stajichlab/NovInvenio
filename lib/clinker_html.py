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


def _strip_gene(gene: dict) -> None:
    for key in STRIP:
        gene.pop(key, None)


def slim_clinker_html(html: str) -> str:
    """The page with sequence/translation removed from every gene in
    data.clusters[].loci[].genes[] and data.links[].query/target. `</` in
    the re-serialised data is escaped so a label cannot close the script.
    Raises ValueError when the page has no `const data=` object."""
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
    return html[:begin] + blob + html[end:]
