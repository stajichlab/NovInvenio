"""
Small reusable HTML/CSS/JS string fragments shared by the NovInvenio report
pages (lib/report_template.py, lib/core_report_template.py,
lib/losses_report_template.py, bin/make_index_report.py).

Colour tokens are *not* defined here -- they come from lib/skins.py, the single
registry every page paints from. This module holds the page chrome (CSS) and
the DOM/link helpers (JS) that sit on top of those tokens.

lib/report_template.py -- the canvas-heatmap novelty report -- historically
duplicated all of this inline. It now shares the skin CSS and the linkout
helpers (the two things that had actually diverged between the reports), while
keeping its own copy of the heatmap-specific chrome: that part is large,
thoroughly covered by tests/test_report_data.py's self-contained-HTML checks,
and has no counterpart in the single-table pages.

Every fragment here is plain text meant to be concatenated into a larger
HTML_TEMPLATE-style triple-quoted string (see lib/report_template.py's module
docstring for the __PROJECT_TITLE__ / /*__PAYLOAD__*/ substitution convention)
-- nothing here does its own token substitution.
"""

import base64
from pathlib import Path

from skins import skin_boot_js, skin_picker_html, skin_picker_js, skins_css

# The full skin registry as CSS. Named for what it is now; the old
# THEME_VARS_CSS spelling is gone along with the light/dark-only model.
SKIN_VARS_CSS = skins_css()

# Logo/favicon, embedded as base64 data URIs -- these pages must open from
# file:// with no network access, so a relative path to assets/logo/ (which
# works fine for NII's own docs/ site, a normal directory tree) won't do here;
# a data: URI is the only way to carry the image inside a single self-contained
# HTML file. Read from the actual asset files (not hardcoded as a giant string
# literal) so the source-of-truth stays the PNG/ICO in assets/logo/, not a
# second copy baked into this module.
_ASSETS_LOGO_DIR = Path(__file__).resolve().parent.parent / "assets" / "logo"


def _data_uri(filename: str, mime: str) -> str:
    data = (_ASSETS_LOGO_DIR / filename).read_bytes()
    return f"data:{mime};base64,{base64.b64encode(data).decode('ascii')}"


FAVICON_DATA_URI = _data_uri("NI_logo_favicon.ico", "image/x-icon")
LOGO_DATA_URI = _data_uri("NI_logo_card-96.png", "image/png")

# <link rel="icon"> tag, ready to drop into any page's <head>.
FAVICON_LINK_HTML = f'<link rel="icon" href="{FAVICON_DATA_URI}">'

# Small header logo, sized to sit next to the h1 inside header.top .titles'
# sibling position (see LOGO_CSS below and each template's <header class="top">).
LOGO_IMG_HTML = f'<img class="logo" src="{LOGO_DATA_URI}" alt="NovInvenio logo">'

LOGO_CSS = r"""
  header.top .logo { width: 40px; height: 40px; border-radius: 8px; flex: 0 0 auto; }
"""

# Top-of-page "up" navigation for the published two-tier gallery layout NII's
# bin/sync_reports.sh/generate_docs.py relocates these pages into --
# docs/<domain>/<set>/*.html sitting under docs/<domain>/index.html sitting
# under docs/index.html (see that repo's DESIGN.md Sec 3/8). Every report page
# a study produces (report.html, novelties.html, core.html, losses.html,
# alignment.html) lives in the same <set> directory, so the relative paths
# below are the same from any of them -- no per-page computation needed.
# Degrades harmlessly (a dead relative link, nothing else on the page depends
# on it) when a page is opened outside that layout, e.g. NovInvenio's own flat
# results/<project>/ copies, or this repo's own -profile test smoke output --
# same self-contained-but-context-optional philosophy as the rest of these
# pages (they already assume no server, no guaranteed sibling files beyond
# what each page embeds itself).
BREADCRUMB_NAV_CSS = r"""
  nav.breadcrumb { display: flex; flex-wrap: wrap; gap: 6px; margin-bottom: 12px; }
  nav.breadcrumb a {
    font-size: 11px; padding: 3px 9px; border: 1px solid var(--border); border-radius: 999px;
    text-decoration: none; color: var(--text-secondary); background: var(--surface-1);
  }
  nav.breadcrumb a:hover { border-color: var(--series-1); color: var(--text-primary); }
"""


def breadcrumb_nav_html(*, study: bool = True) -> str:
    """The nav bar itself: Study folder -> Group gallery -> All studies (site root).

    ``study=False`` omits the study-folder link -- pass it on report.html
    itself, where linking to "report.html" would just point at the current
    page.
    """
    links = []
    if study:
        links.append('<a href="report.html">&#128193; Study</a>')
    links.append('<a href="../index.html">&#128194; Group</a>')
    links.append('<a href="../../index.html">&#127968; All studies</a>')
    return '<nav class="breadcrumb">' + ''.join(links) + '</nav>'

# <head> snippet -- must run before first paint so a stored skin choice does
# not flash the default palette. Wrap in <script>...</script> at the call site.
SKIN_BOOT_JS = skin_boot_js()

# End-of-body wiring for the header's skin <select>.
SKIN_PICKER_JS = skin_picker_js()

# The header control itself.
SKIN_PICKER_HTML = skin_picker_html()

# Page chrome, filter bar, single data table and detail panel -- the shape
# every single-table report needs. No canvas/heatmap rules here; that stays
# specific to lib/report_template.py.
BASE_PAGE_CSS = r"""
  * { box-sizing: border-box; }
  body {
    margin: 0;
    /* Skin-supplied page texture (scanlines for the neon skin, `none`
       everywhere else) painted as a background layer on the page ground.
       An element's background always paints below all of its descendants, so
       the texture can never land on a card, a table, the detail panel or a
       tooltip -- it shows in the page's negative space only. An overlay
       element positioned above the content cannot make that guarantee without
       every content region opting out by z-index, which is how the first
       attempt at this put scanlines across the detail panel and the heatmap
       tooltip. Deliberately static, never animated. */
    background: var(--page) var(--overlay);
    color: var(--text-primary);
    font: 14px/1.5 var(--font-ui);
  }
  .wrap { max-width: min(95vw, 2100px); margin: 0 auto; padding: 24px 20px 64px; }
  header.top { display: flex; align-items: flex-start; gap: 16px; margin-bottom: 24px; }
  header.top .titles { flex: 1; min-width: 0; }
  h1 { margin: 0 0 4px; font-size: 22px; font-weight: 600; letter-spacing: -0.01em; text-shadow: var(--glow); }
  .sub { margin: 0; color: var(--text-secondary); font-size: 13px; }
  .mono { font-family: var(--font-mono); font-size: 11px; }
  button, select, input[type="search"] {
    font: inherit;
    color: var(--text-primary);
    background: var(--surface-1);
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 6px 10px;
  }
  button { cursor: pointer; }
  button:hover { background: var(--hover-wash); }
  button:focus-visible, select:focus-visible, input:focus-visible {
    outline: 2px solid var(--series-1);
    outline-offset: 2px;
  }
  .btn-ghost { background: transparent; }
  section.card {
    background: var(--surface-1);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 18px;
    margin-bottom: 16px;
  }
  .card-title { margin: 0 0 2px; font-size: 14px; font-weight: 600; }
  .card-note { margin: 0 0 16px; color: var(--text-secondary); font-size: 12px; }
  .tiles { display: flex; flex-wrap: wrap; gap: 28px; }
  .tile-label { color: var(--text-secondary); font-size: 12px; }
  .tile-value { font-size: 20px; font-weight: 600; }
  .filters { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; margin-bottom: 14px; }
  .filters input[type="search"] { min-width: 260px; }
  .spacer { flex: 1; }
  .count { color: var(--text-secondary); font-size: 13px; font-variant-numeric: tabular-nums; }
  .explorer { display: grid; grid-template-columns: 1fr 380px; gap: 16px; align-items: start; }
  @media (max-width: 1180px) { .explorer { grid-template-columns: 1fr; } }
  .tbl-scroll { overflow: auto; height: 560px; }
  table.data { border-collapse: collapse; width: 100%; font-size: 12px; }
  table.data th, table.data td {
    text-align: left; padding: 6px 10px; border-bottom: 1px solid var(--grid); white-space: nowrap;
  }
  table.data th {
    position: sticky; top: 0; background: var(--surface-1); font-weight: 600; z-index: 1;
    border-bottom: 1px solid var(--axis);
  }
  table.data th.sortable { cursor: pointer; user-select: none; }
  table.data th.sortable:hover { color: var(--series-1); }
  table.data th.sortable::after { content: ""; margin-left: 4px; color: var(--series-1); }
  table.data th.sortable.sort-active::after { content: "\25BE"; }
  table.data td.num { font-variant-numeric: tabular-nums; }
  /* Text columns grow with the viewport (floor keeps them readable, ceiling bounded). */
  table.data td.wrap-cell, table.data th.wrap-cell {
    white-space: normal; overflow-wrap: anywhere;
    min-width: 200px; max-width: clamp(260px, 26vw, 640px);
  }
  table.data td.cell, table.data th.cell { padding: 6px 5px; text-align: center; }
  table.data tbody tr:hover { background: var(--hover-wash); }
  table.data tbody tr.sel { background: var(--wash); }
  .empty { padding: 40px 8px; color: var(--text-secondary); font-size: 13px; }
  .detail { position: sticky; top: 16px; }
  .detail h3 { margin: 0 0 2px; font-size: 15px; font-weight: 600; word-break: break-all; }
  .detail .species { color: var(--text-secondary); font-size: 12px; font-style: italic; margin-bottom: 14px; }
  .field { margin-bottom: 12px; }
  .field-label { font-size: 11px; text-transform: uppercase; letter-spacing: 0.04em; color: var(--muted); margin-bottom: 3px; }
  .field-value { font-size: 13px; word-break: break-word; }
  .chips { display: flex; flex-wrap: wrap; gap: 5px; }
  .chip {
    display: inline-block; padding: 2px 7px; border: 1px solid var(--border); border-radius: 999px;
    font-size: 11px; text-decoration: none; color: var(--text-primary); background: var(--page);
  }
  .chip:hover { border-color: var(--series-1); }
  a.pfam-link { color: var(--series-1); text-decoration: none; }
  a.pfam-link:hover { text-decoration: underline; }
  .links { display: flex; flex-wrap: wrap; gap: 6px; }
  .links a, .links button {
    font-size: 12px; padding: 5px 9px; border: 1px solid var(--border); border-radius: 6px;
    text-decoration: none; color: var(--text-primary); background: var(--page);
  }
  .links a:hover, .links button:hover { border-color: var(--series-1); }
  /* The "this protein has no annotation" follow-up cluster -- the population
     this pipeline exists to surface, so it gets its own visible grouping
     rather than being one more link in the row. */
  .links-note {
    font-size: 11px; color: var(--muted); margin: 8px 0 4px;
  }
  .badge {
    display: inline-block; padding: 1px 6px; border-radius: 999px; font-size: 10.5px;
    border: 1px solid var(--border); color: var(--text-secondary);
  }
  .badge.warn { color: var(--warn); border-color: var(--warn); }
  .placeholder { color: var(--text-secondary); font-size: 13px; }
  .hidden { display: none !important; }
  .sr-only {
    position: absolute; width: 1px; height: 1px; padding: 0; margin: -1px;
    overflow: hidden; clip: rect(0 0 0 0); white-space: nowrap; border: 0;
  }
  @media print {
    .filters, .tabs, header.top select, header.top button { display: none !important; }
    .tbl-scroll { height: auto; overflow: visible; }
    .explorer { grid-template-columns: 1fr; }
  }
"""

# ---- JS fragments ---------------------------------------------------------

# Untrusted-string discipline: every report must build DOM nodes with this
# helper (or plain textContent), never innerHTML -- see CLAUDE.md's report
# constraints. Protein IDs and annotation text come from FASTA headers and
# SwissProt/Pfam text and are not sanitised upstream.
EL_HELPER_JS = r"""
  function el(tag, cls, text) {
    var n = document.createElement(tag);
    if (cls) n.className = cls;
    if (text !== undefined && text !== null) n.textContent = text; // labels are untrusted data
    return n;
  }
  function extLink(label, href, title) {
    var a = el("a", null, label);
    a.href = href;
    a.target = "_blank";
    a.rel = "noopener noreferrer";
    if (title) a.title = title;
    return a;
  }
  // Candidate protein IDs are only unique within their own source proteome
  // (MAG/prodigal locus tags in particular collide across MAGs -- e.g.
  // "k141_81591_30" says nothing about which of the study's MAGs it came
  // from), so any human-facing label needs the source proteome's Short
  // prefixed. Some pipelines already bake a "<Short>__" prefix into
  // protein_id itself (see NII's bin/build_study_config.py header-rewrite
  // convention) -- id.indexOf check makes this a no-op for those instead of
  // doubling up ("UHM102bin47:UHM102bin47__k141..."). Never used for the raw
  // id itself (FASTA export, UniProt-record parsing, sort/search keys) --
  // only for what a reader sees.
  function displayId(id, sp) {
    if (!sp) return id;
    return id.indexOf(sp.short) === 0 ? id : sp.short + ":" + id;
  }
"""

# Every report resolves a protein to the same external record, so this lives in
# one place. lib/report_template.py imports it too -- the inline copies it and
# the single-table templates each carried had already drifted apart (the core
# and losses pages silently dropped the Pfam E-values the novelty page shows).
LINKOUT_HELPERS_JS = r"""
  function geneIdFromProteinId(pid) {
    return pid.replace(/-[Tt][^-]*(-p\d+)?$|-p\d+$/, "");
  }
  function uniprotAcc(sprot) {
    if (!sprot) return "";
    var m = /(?:^|\|)([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:\||\s|$)/.exec(sprot);
    return m ? m[1] : "";
  }
  // Best SwissProt hit as a UniProt hyperlink -- falls back to plain text when
  // the accession can't be parsed out of the hit string.
  function uniprotLinkNode(sprot) {
    var acc = uniprotAcc(sprot);
    if (!acc) return document.createTextNode(sprot || "");
    var a = el("a", "pfam-link", sprot);
    a.href = "https://www.uniprot.org/uniprotkb/" + acc + "/entry";
    a.target = "_blank"; a.rel = "noopener noreferrer";
    a.title = "UniProt " + acc;
    return a;
  }
  // Resolve an interned-string index (payload['descriptions']/'go_sets'/'ipr_sets' --
  // -1 for "no value", see lib/report_data.py's _StringTable) back to its string.
  function fromTable(table, i) {
    return i >= 0 && table ? (table[i] || "") : "";
  }
  // Resolve a search-hit target protein ID (row[F.tgt], comma-split per proteome --
  // see lib/report_data.py's ROW_FIELDS 'tgt' docs) to a display label via
  // payload['protein_names'] (built from bin/extract_protein_descriptions.py's
  // output). Falls back to the bare ID when there's no matching entry (not every
  // target_id resolves -- non-UniProt-sourced proteomes, or a run with no
  // --descriptions given at all), and to "" when there's no target_id at all.
  function targetLabel(targetId) {
    if (!targetId) return "";
    var info = DATA.protein_names && DATA.protein_names[targetId];
    if (!info) return targetId;
    var parts = [info.gene_name, info.description].filter(Boolean);
    return parts.length ? parts.join(" ") : targetId;
  }
  // Format a numeric-string E-value to 4 significant figures for display (e.g.
  // "4.549999999999999e-230" -> "4.550e-230") -- the raw strings survive
  // float->str round-tripping through the TSV sidecars with full double
  // precision, which is noise no report reader needs. Returns the input
  // unchanged if it doesn't parse as a number (defensive; shouldn't happen).
  function fmtEvalue(ev) {
    if (!ev) return "";
    var n = Number(ev);
    if (!isFinite(n)) return ev;
    return n.toPrecision(4);
  }
  // Pfam names as a chip row (detail panel) -- accession + optional E-value in the
  // link title, same InterPro/Pfam target used by every report's "Pfam domains" field.
  function pfamChipsNode(pfamNames, pfamAccs, pfamEvs) {
    var chips = el("div", "chips");
    var names = pfamNames ? pfamNames.split(",") : [];
    var accs = pfamAccs ? pfamAccs.split(",") : [];
    var evs = pfamEvs ? pfamEvs.split(",") : [];
    names.forEach(function (n, i) {
      var acc = (accs[i] || "").split(".")[0];
      var ev = fmtEvalue(evs[i]);
      var node;
      if (/^PF\d+$/.test(acc)) {
        node = el("a", "chip", n + (ev ? " · " + ev : ""));
        node.href = "https://www.ebi.ac.uk/interpro/entry/pfam/" + acc + "/";
        node.target = "_blank"; node.rel = "noopener noreferrer";
        node.title = acc + (ev ? " — hmmsearch E-value " + ev : "");
      } else {
        node = el("span", "chip", n);
      }
      chips.appendChild(node);
    });
    return chips;
  }
  // GO terms as a chip row (detail panel) -- payload's 'go' field is "|"-separated
  // "GO:nnnnnnn:EVIDENCE" entries (NII bin/extract_dat_annotations.py, curated from
  // the source UniProt record's own DR GO cross-references). Links to QuickGO.
  function goChipsNode(goIds) {
    var chips = el("div", "chips");
    (goIds ? goIds.split("|") : []).filter(Boolean).forEach(function (entry) {
      var parts = entry.split(":");
      var id = parts[0] + ":" + parts[1];
      var evidence = parts[2] || "";
      var a = el("a", "chip", id + (evidence ? " · " + evidence : ""));
      a.href = "https://www.ebi.ac.uk/QuickGO/term/" + id;
      a.target = "_blank"; a.rel = "noopener noreferrer";
      a.title = id + (evidence ? " — evidence code " + evidence : "");
      chips.appendChild(a);
    });
    return chips;
  }
  // InterPro domains as a chip row (detail panel) -- payload's 'ipr' field is
  // "|"-separated "IPRnnnnnn" entries (same UniProt DR cross-reference source).
  function interproChipsNode(iprIds) {
    var chips = el("div", "chips");
    (iprIds ? iprIds.split("|") : []).filter(Boolean).forEach(function (id) {
      var a = el("a", "chip", id);
      a.href = "https://www.ebi.ac.uk/interpro/entry/InterPro/" + id + "/";
      a.target = "_blank"; a.rel = "noopener noreferrer";
      a.title = id;
      chips.appendChild(a);
    });
    return chips;
  }
  // EC numbers as a chip row (detail panel) -- payload's 'ec' field is
  // comma-separated (bin/extract_dat_annotations.py's DE-line EC=... values,
  // matching the Pfam fields' separator convention). Links to ExPASy ENZYME.
  function ecChipsNode(ecNumbers) {
    var chips = el("div", "chips");
    (ecNumbers ? ecNumbers.split(",") : []).filter(Boolean).forEach(function (ec) {
      var a = el("a", "chip", "EC " + ec);
      a.href = "https://enzyme.expasy.org/EC/" + ec;
      a.target = "_blank"; a.rel = "noopener noreferrer";
      a.title = "ExPASy ENZYME: EC " + ec;
      chips.appendChild(a);
    });
    return chips;
  }
  // AlphaFold predicted-structure link (detail panel) -- payload's 'af' field is
  // a single AlphaFold DB accession (bin/extract_dat_annotations.py's DR
  // AlphaFoldDB cross-reference), or '' when absent. AlphaFold DB covers nearly
  // all of UniProt, so this is usually available and is a real structure link,
  // not the generic "no Pfam/SwissProt -- try a remote-homology search" fallback.
  function alphafoldLinkNode(afId) {
    var a = el("a", "chip", "View predicted structure");
    a.href = "https://alphafold.ebi.ac.uk/entry/" + afId;
    a.target = "_blank"; a.rel = "noopener noreferrer";
    a.title = "AlphaFold DB: " + afId;
    return a;
  }
  // The detail panel's protein-ID heading as a UniProt hotlink, when the ID is a
  // UniProt FASTA header token ("sp|ACC|NAME" or "tr|ACC|NAME" -- see NII's
  // lib/uniprot_ids.bare_accession, the same format this parses). Returns null
  // for a non-UniProt protein_id (e.g. a BFD/funannotate-style gene ID), so the
  // caller can fall back to plain text -- never assume every study is UniProt-sourced.
  function uniprotRecordLinkNode(proteinId) {
    var m = /^(?:sp|tr)\|([^|]+)\|/.exec(proteinId);
    if (!m) return null;
    var a = el("a", null, proteinId);
    a.href = "https://www.uniprot.org/uniprotkb/" + m[1] + "/entry";
    a.target = "_blank"; a.rel = "noopener noreferrer";
    a.title = "View " + m[1] + " on UniProt";
    return a;
  }
  // Compact comma-separated Pfam links for a table cell (same accession rule as
  // pfamChipsNode, without the chip styling -- a table row is dense already).
  function pfamLinksInline(pfamNames, pfamAccs) {
    var frag = document.createDocumentFragment();
    var names = pfamNames ? pfamNames.split(",") : [];
    var accs = pfamAccs ? pfamAccs.split(",") : [];
    names.forEach(function (n, i) {
      if (i > 0) frag.appendChild(document.createTextNode(", "));
      var acc = (accs[i] || "").split(".")[0];
      if (/^PF\d+$/.test(acc)) {
        var a = el("a", "pfam-link", n);
        a.href = "https://www.ebi.ac.uk/interpro/entry/pfam/" + acc + "/";
        a.target = "_blank"; a.rel = "noopener noreferrer";
        a.title = acc;
        frag.appendChild(a);
      } else {
        frag.appendChild(document.createTextNode(n));
      }
    });
    return frag;
  }

  // ---- genome-database gene records --------------------------------------
  // Driven by the config CSV's optional SourceDB column, carried through to
  // each payload proteome entry as `source_db` (see lib/config_parser.py).
  // Before this, a gene-record link was gated on the annotation source being
  // a ModelOrg_* entry, so only the handful of model organisms got one even
  // though most of the sample pool lives in FungiDB or MycoCosm.
  //
  // Accepted forms: "fungidb", "mycocosm:<portal>", "ensemblfungi:<species>",
  // "veupathdb:<project>", "ncbipep", or any URL template containing "{gene}".
  // A config CSV travels between users and projects, so a SourceDB value is
  // only semi-trusted: anything that ends up in an href must be checked for
  // its scheme, or "javascript:...{gene}" becomes a clickable link in the
  // report. The keyed forms below interpolate into a fixed https:// prefix,
  // so only the free-form template needs the check.
  function isSafeHttpUrl(url) {
    return /^https?:\/\//i.test(url);
  }
  function genomeDbLink(sourceDb, proteinId) {
    if (!sourceDb) return null;
    var gene = geneIdFromProteinId(proteinId);
    if (sourceDb.indexOf("{gene}") >= 0) {
      if (!isSafeHttpUrl(sourceDb)) return null;
      return extLink("Gene record", sourceDb.replace("{gene}", encodeURIComponent(gene)),
                     "Source database record for " + gene);
    }
    var sep = sourceDb.indexOf(":");
    var kind = (sep >= 0 ? sourceDb.slice(0, sep) : sourceDb).toLowerCase();
    var arg = sep >= 0 ? sourceDb.slice(sep + 1) : "";
    if (kind === "fungidb") {
      return extLink("FungiDB gene",
        "https://fungidb.org/fungidb/app/record/gene/" + encodeURIComponent(gene),
        "FungiDB gene record for " + gene);
    }
    if (kind === "mycocosm" && arg) {
      return extLink("JGI MycoCosm",
        "https://mycocosm.jgi.doe.gov/cgi-bin/dispGeneModel?db=" +
          encodeURIComponent(arg) + "&id=" + encodeURIComponent(proteinId),
        "JGI MycoCosm gene model in portal " + arg);
    }
    if (kind === "ensemblfungi" && arg) {
      return extLink("Ensembl Fungi",
        "https://fungi.ensembl.org/" + encodeURIComponent(arg) +
          "/Gene/Summary?g=" + encodeURIComponent(gene),
        "Ensembl Fungi gene record for " + gene);
    }
    if (kind === "veupathdb" && arg) {
      return extLink("VEuPathDB gene",
        "https://" + encodeURIComponent(arg) + ".org/" + encodeURIComponent(arg) +
          "/app/record/gene/" + encodeURIComponent(gene),
        "VEuPathDB gene record for " + gene);
    }
    if (kind === "ncbipep") {
      // Uses proteinId as-is, NOT geneIdFromProteinId(proteinId) -- NCBI Protein
      // pages are keyed by the literal accession (e.g. "NP_001020957.2"), and
      // stripping a FungiDB-style transcript/protein suffix here would be wrong
      // for species with no such suffix convention (a no-op for RefSeq
      // accessions today, but the two ID spaces are not interchangeable).
      return extLink("NCBI Protein",
        "https://www.ncbi.nlm.nih.gov/protein/" + encodeURIComponent(proteinId),
        "NCBI Protein record for " + proteinId);
    }
    return null;
  }

  function taxonomyLink(proteome) {
    if (!proteome) return null;
    if (proteome.taxid) {
      return extLink("NCBI Taxonomy",
        "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=" +
          encodeURIComponent(proteome.taxid),
        "NCBI Taxonomy taxid " + proteome.taxid);
    }
    if (proteome.species) {
      return extLink("NCBI Taxonomy",
        "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name=" +
          encodeURIComponent(proteome.species),
        "NCBI Taxonomy search for " + proteome.species);
    }
    return null;
  }

  // ---- sequence-driven links ---------------------------------------------
  function fastaText(pid, seq) {
    return ">" + pid + "\n" + (seq.match(/.{1,60}/g) || []).join("\n") + "\n";
  }
  function copyFastaButton(pid, seq) {
    var cp = el("button", null, "Copy FASTA");
    cp.type = "button";
    cp.addEventListener("click", function () {
      var fa = fastaText(pid, seq);
      // navigator.clipboard is undefined on some file:// origins; fall back to
      // selecting a scratch textarea so the button is never a dead control.
      function done() {
        cp.textContent = "Copied";
        setTimeout(function () { cp.textContent = "Copy FASTA"; }, 1400);
      }
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(fa).then(done, function () { legacyCopy(fa, done); });
      } else {
        legacyCopy(fa, done);
      }
    });
    return cp;
  }
  function legacyCopy(text, done) {
    var ta = document.createElement("textarea");
    ta.value = text;
    ta.setAttribute("readonly", "");
    ta.style.position = "fixed";
    ta.style.opacity = "0";
    document.body.appendChild(ta);
    ta.select();
    try { document.execCommand("copy"); done(); } catch (e) { /* leave label alone */ }
    ta.remove();
  }
  // NCBI BLAST accepts the query in the URL, but a long protein blows past the
  // ~8 KB practical URL ceiling and the request silently fails in some
  // browsers and proxies. Above the threshold, submit the same query as a POST
  // form instead -- still no external assets, still works from file://.
  var BLAST_GET_MAX_AA = 1200;
  function ncbiBlastpNode(pid, seq) {
    var base = "https://blast.ncbi.nlm.nih.gov/Blast.cgi";
    if (seq.length <= BLAST_GET_MAX_AA) {
      return extLink("BLASTP at NCBI",
        base + "?PAGE=Proteins&PROGRAM=blastp&BLAST_PROGRAMS=blastp" +
        "&DATABASE=nr&CMD=Web&QUERY=" + encodeURIComponent(seq),
        "blastp against NCBI nr");
    }
    var b = el("button", null, "BLASTP at NCBI");
    b.type = "button";
    b.title = "Sequence is " + seq.length + " aa — submitted as a POST form";
    b.addEventListener("click", function () {
      var f = document.createElement("form");
      f.method = "POST";
      f.action = base;
      f.target = "_blank";
      f.rel = "noopener noreferrer";
      var fields = { CMD: "Put", PROGRAM: "blastp", DATABASE: "nr", QUERY: fastaText(pid, seq) };
      Object.keys(fields).forEach(function (k) {
        var i = document.createElement("input");
        i.type = "hidden"; i.name = k; i.value = fields[k];
        f.appendChild(i);
      });
      document.body.appendChild(f);
      f.submit();
      f.remove();
    });
    return b;
  }
"""

# Every report resolves a protein to the same external record via this one
# builder, so this lives alongside LINKOUT_HELPERS_JS (imported/concatenated
# together at each template's <script> assembly site) rather than inside it,
# purely so the DR-cross-reference registry has its own named block.
EXTERNAL_LINKS_JS = r"""
  // ---- the whole "External resources" block ------------------------------
  // One builder for all three reports so they can no longer drift apart.
  //   o.id        protein ID
  //   o.gene      gene_name ('' if none)
  //   o.geneUrl   Model_Org_Gene_URL -- modelorgs.yaml's gene_url_template
  //               resolved by bin/annotate_presence_matrix.py against the
  //               gene_names_csv LOOKUP KEY (e.g. a UniProt accession), not
  //               against o.gene itself (see ModelOrgAnnotator.gene_url()'s
  //               docstring for why those can differ), or '' if unset/unresolved
  //   o.sprot     Best_Swissprot ('' if none)
  //   o.pfam      Pfam_Names ('' if none)
  //   o.fsrcName  annotation source label ('' if none)
  //   o.seq       protein sequence ('' when the payload carries no sequences)
  //   o.proteome  payload proteomes[] entry for the row's source species
  //   o.xrefs     UniProt DR-line cross-references, "DB:id|DB:id" ('' if none)
  // ---- generic UniProt DR-cross-reference linkouts ------------------------
  // o.xrefs is "DB:id|DB:id" (NII bin/extract_dat_annotations.py's DR-line
  // extraction, lib/report_data.py's 'xrefs' field) -- split each entry on
  // the FIRST colon only (never a fixed-arity split like goChipsNode's
  // go_ids parsing): KEGG ("ncr:NCU10683") and VEuPathDB ("FungiDB:NCU10683")
  // ids carry their own colons and must not be truncated.
  var VEUPATHDB_HOSTS = {
    FungiDB: "fungidb.org/fungidb",
    PlasmoDB: "plasmodb.org/plasmo",
    ToxoDB: "toxodb.org/toxo"
  };
  var XREF_LINK_TEMPLATES = {
    VEuPathDB: {
      render: function (id) {
        var i = id.indexOf(":");
        var project = i === -1 ? "VEuPathDB" : id.slice(0, i);
        var geneId = i === -1 ? id : id.slice(i + 1);
        var host = VEUPATHDB_HOSTS[project] || ("veupathdb.org/" + project.toLowerCase());
        return {
          label: project + " " + geneId,
          href: "https://" + host + "/app/record/gene/" + encodeURIComponent(geneId),
          title: project + " gene record for " + geneId
        };
      }
    },
    GeneID: {
      render: function (id) {
        return { label: "NCBI Gene " + id,
                 href: "https://www.ncbi.nlm.nih.gov/gene/" + encodeURIComponent(id),
                 title: "NCBI Gene " + id };
      }
    },
    RefSeq: {
      render: function (id) {
        return { label: "RefSeq " + id,
                 href: "https://www.ncbi.nlm.nih.gov/protein/" + encodeURIComponent(id),
                 title: "RefSeq protein " + id };
      }
    },
    KEGG: {
      render: function (id) {
        return { label: "KEGG " + id,
                 href: "https://www.genome.jp/dbget-bin/www_bget?" + encodeURIComponent(id),
                 title: "KEGG entry " + id };
      }
    },
    EnsemblFungi: {
      render: function (id) {
        return { label: "Ensembl Fungi " + id,
                 href: "https://fungi.ensembl.org/id/" + encodeURIComponent(id),
                 title: "Ensembl Fungi gene " + id };
      }
    },
    EnsemblBacteria: {
      render: function (id) {
        return { label: "Ensembl Bacteria " + id,
                 href: "https://bacteria.ensembl.org/id/" + encodeURIComponent(id),
                 title: "Ensembl Bacteria gene " + id };
      }
    }
  };
  // Returns {links: [<a> node, ...], hasVEuPathDB: bool} -- the caller (see
  // externalLinksNode below) needs to know whether a VEuPathDB xref fired so
  // it can suppress genomeDbLink's separate, UniProt-accession-keyed FungiDB
  // link for the same row (that one resolves against o.id, which for a
  // UniProt-sourced protein is the wrong ID space -- see the dedup comment in
  // externalLinksNode below).
  function xrefLinkNodes(xrefsStr) {
    var links = [];
    var hasVEuPathDB = false;
    (xrefsStr ? xrefsStr.split("|") : []).filter(Boolean).forEach(function (entry) {
      var i = entry.indexOf(":");
      if (i === -1) return;
      var db = entry.slice(0, i);
      var id = entry.slice(i + 1);
      var tmpl = XREF_LINK_TEMPLATES[db];
      if (!tmpl || !id) return;
      var r = tmpl.render(id);
      links.push(extLink(r.label, r.href, r.title));
      if (db === "VEuPathDB") hasVEuPathDB = true;
    });
    return { links: links, hasVEuPathDB: hasVEuPathDB };
  }
  function externalLinksNode(o) {
    var box = document.createDocumentFragment();
    var links = el("div", "links");
    var acc = uniprotAcc(o.sprot);
    if (acc) {
      links.appendChild(extLink("UniProt " + acc,
        "https://www.uniprot.org/uniprotkb/" + acc + "/entry"));
      links.appendChild(extLink("AlphaFold",
        "https://alphafold.ebi.ac.uk/entry/" + acc,
        "Predicted structure for " + acc));
    }
    var db = genomeDbLink(o.proteome && o.proteome.source_db, o.id);
    var xrefResult = xrefLinkNodes(o.xrefs);
    // FungiDB dedup: genomeDbLink's "fungidb" branch resolves against
    // geneIdFromProteinId(o.id), which for a UniProt-sourced protein is the
    // UniProt accession, not a real FungiDB gene ID -- that link is already
    // broken for these rows. When the xrefs-derived VEuPathDB link is
    // present, it is the correct one; suppress genomeDbLink's instead of
    // showing both (one dead) side by side.
    var suppressGenomeDbFungiDb = xrefResult.hasVEuPathDB && db &&
      /fungidb\.org/.test(db.href || "");
    if (db && !suppressGenomeDbFungiDb) links.appendChild(db);
    xrefResult.links.forEach(function (a) { links.appendChild(a); });
    // Model-org gene lookup (modelorgs.yaml) is a *different* external record
    // than the candidate's own genome-db entry above: it's a lookup into
    // whatever database the model organism's gene_names_csv itself came from
    // (UniProt, FungiDB, ...), and o.geneUrl is already resolved against the
    // exact lookup key (see the o.geneUrl field note above) -- never on the
    // candidate's own protein_id, which a MAG's prodigal-called ID has no
    // relationship to at all. Only rendered when the modelorgs.yaml entry
    // that fired sets gene_url_template; there is no safe default to guess,
    // since the same fsrcName ("ModelOrg_<short>") fires for both FungiDB-
    // and UniProt-backed configs.
    if (o.gene && o.geneUrl && isSafeHttpUrl(o.geneUrl)) {
      links.appendChild(extLink("Model organism gene: " + o.gene, o.geneUrl,
        "Reference gene record for " + o.gene));
    }
    var tax = taxonomyLink(o.proteome);
    if (tax) links.appendChild(tax);
    if (!acc && !db) {
      var term = o.gene || geneIdFromProteinId(o.id);
      links.appendChild(extLink("Search NCBI Protein",
        "https://www.ncbi.nlm.nih.gov/protein/?term=" + encodeURIComponent(term),
        "No SwissProt hit or source-database link — worth a manual check"));
    }
    if (o.seq) {
      // Copy first: every tool below it is a paste target, so the old order
      // (paste-me links above the button that fills the clipboard) was
      // backwards.
      links.appendChild(copyFastaButton(o.id, o.seq));
      links.appendChild(ncbiBlastpNode(o.id, o.seq));
      links.appendChild(extLink("UniProt BLAST",
        "https://www.uniprot.org/blast?query=" + encodeURIComponent(o.seq),
        "blastp against UniProtKB"));
    }
    box.appendChild(links);

    // A candidate with neither a Pfam domain nor a SwissProt hit is the whole
    // point of this pipeline, and it is exactly the row for which an ID-based
    // NCBI search returns nothing. Give it the remote-homology and structure
    // tools that are the real next step.
    if (o.seq && !acc && !o.pfam) {
      box.appendChild(el("p", "links-note",
        "No Pfam domain and no SwissProt hit — remote-homology and structure searches:"));
      var more = el("div", "links");
      more.appendChild(extLink("HHpred",
        "https://toolkit.tuebingen.mpg.de/tools/hhpred",
        "MPI Bioinformatics Toolkit — paste the copied FASTA"));
      more.appendChild(extLink("Foldseek",
        "https://search.foldseek.com/search",
        "Structure search — paste the copied FASTA"));
      more.appendChild(extLink("InterProScan",
        "https://www.ebi.ac.uk/interpro/search/sequence/",
        "InterProScan — paste the copied FASTA"));
      box.appendChild(more);
    }
    return box;
  }
"""

DOWNLOAD_JS = r"""
  function download(name, text, type) {
    var blob = new Blob([text], { type: type });
    var url = URL.createObjectURL(blob);
    var a = document.createElement("a");
    a.href = url; a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 0);
  }
"""

# TBLASTN alignment popup -- docs/-only (see CLAUDE.md's report constraints):
# the results/ copy stays fully embedded/file://-safe and never includes this
# fragment. Callers gate inclusion behind their own `online` flag (see #75);
# nothing here decides that for itself. bin/build_alignment_shards.py (#72)
# is the producer: one gzip JSON shard per query genome under
# alignments/<genome>.json.gz (novelty) or loss_alignments/<genome>.json.gz
# (loss), keyed protein_id -> [hit, ...], schema_version-stamped in a sibling
# manifest.json.
ALIGNMENT_POPUP_CSS = r"""
  dialog.alignment {
    max-width: min(92vw, 900px);
    width: 100%;
    max-height: 85vh;
    padding: 0;
    border: 1px solid var(--border);
    border-radius: 10px;
    background: var(--surface-1);
    color: var(--text-primary);
    box-shadow: var(--shadow);
  }
  dialog.alignment::backdrop { background: rgba(0, 0, 0, 0.45); }
  dialog.alignment .align-head {
    display: flex; align-items: center; gap: 10px;
    padding: 14px 18px; border-bottom: 1px solid var(--border);
  }
  dialog.alignment .align-head h3 { margin: 0; font-size: 14px; font-weight: 600; flex: 1; min-width: 0; }
  dialog.alignment .align-head h3 .align-uniprot-link {
    font-weight: 400; font-size: 12px; color: var(--text-secondary); text-decoration: none;
  }
  dialog.alignment .align-head h3 .align-uniprot-link:hover { text-decoration: underline; }
  dialog.alignment .align-stats {
    padding: 10px 18px; font-size: 12px; color: var(--text-secondary);
    border-bottom: 1px solid var(--border);
  }
  dialog.alignment .align-body { padding: 12px 18px; overflow: auto; max-height: 60vh; }
  dialog.alignment pre.align-block {
    font-family: var(--font-mono); font-size: 11px; line-height: 1.5;
    white-space: pre; margin: 0 0 10px;
  }
  dialog.alignment .align-empty { padding: 24px 4px; color: var(--text-secondary); font-size: 13px; }
"""

# Bare <dialog>; content is populated entirely via textContent at render time
# (protein IDs, sequences and hit metadata are untrusted strings -- FASTA
# headers and BLAST output -- see CLAUDE.md's report constraints).
ALIGNMENT_POPUP_HTML = r"""
  <dialog class="alignment" id="alignment-dialog">
    <div class="align-head">
      <h3 id="alignment-title"></h3>
      <button type="button" class="btn-ghost" id="alignment-close" aria-label="Close">&times;</button>
    </div>
    <div class="align-stats" id="alignment-stats"></div>
    <div class="align-body" id="alignment-body"></div>
  </dialog>
"""

# Public surface: window.NIAlignments.open(baseUrl, genome, proteinId, hitIndex).
# baseUrl is the caller's own relative path to its shard directory
# ("alignments/" or "loss_alignments/"), so this fragment stays agnostic to
# which direction (novelty/loss) or which page is calling it.
ALIGNMENT_POPUP_JS = r"""
  (function () {
    var shardCache = new Map(); // url -> Promise<Object>, cached by URL so two
                                 // rapid clicks against the same genome never
                                 // fire a second fetch (see #74 design notes).

    function isGzip(bytes) {
      return bytes.length >= 2 && bytes[0] === 0x1f && bytes[1] === 0x8b;
    }

    async function fetchDataShard(url) {
      if (shardCache.has(url)) return shardCache.get(url);
      var promise = (async function () {
        var resp = await fetch(url);
        if (!resp.ok) throw new Error("fetch " + url + ": " + resp.status);
        var bytes = new Uint8Array(await resp.arrayBuffer());
        var text;
        if (isGzip(bytes)) {
          if (typeof DecompressionStream === "undefined") {
            throw new Error("DECOMPRESSION_UNSUPPORTED");
          }
          var stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
          text = await new Response(stream).text();
        } else {
          // Some hosts transparently decode gzip in transit (Content-Encoding) --
          // if the bytes we got back aren't gzip-magic, treat them as already-
          // decoded plain text rather than double-decompressing and throwing.
          text = new TextDecoder("utf-8").decode(bytes);
        }
        return JSON.parse(text);
      })();
      shardCache.set(url, promise);
      return promise;
    }

    // Positional match/mismatch line. tblastn is gapped by default, so qseq/
    // sseq are always equal length (padded with '-') -- never assume ungapped.
    function midline(qseq, sseq) {
      var out = "";
      for (var i = 0; i < qseq.length; i++) {
        var a = qseq[i], b = sseq[i];
        out += (a === "-" || b === "-") ? " " : (a === b ? "|" : " ");
      }
      return out;
    }

    // Query coordinates advance left-to-right per block (proteins have no
    // strand, so they always increase). Subject coordinates are genomic
    // nucleotide positions -- each residue of sseq corresponds to one codon
    // (3 nt) -- and move in the direction given by sframe: increasing for a
    // plus-strand hit (sframe > 0), decreasing for a minus-strand hit
    // (sframe < 0). Gap columns ('-') consume a query/subject residue on one
    // side only, so per-block spans are computed from each side's own
    // non-gap character count, not the raw slice width.
    function chunkAlignment(qseq, sseq, qstart, sstart, sframe, width) {
      width = width || 60;
      var mid = midline(qseq, sseq);
      var strand = sframe < 0 ? -1 : 1;
      var blocks = [];
      var qi = qstart;
      var si = sstart;
      for (var off = 0; off < qseq.length; off += width) {
        var qc = qseq.slice(off, off + width);
        var sc = sseq.slice(off, off + width);
        var mc = mid.slice(off, off + width);
        var qConsumed = qc.replace(/-/g, "").length;
        var sConsumed = sc.replace(/-/g, "").length;
        var sSpanNt = sConsumed * 3;
        var qBlockStart = qi;
        var qBlockEnd = qi + qConsumed - 1;
        var sBlockStart = si;
        var sBlockEnd = strand > 0 ? si + sSpanNt - 1 : si - sSpanNt + 1;
        blocks.push({
          qstart: qBlockStart, qend: qBlockEnd, qseq: qc, mid: mc, sseq: sc,
          sstart: sBlockStart, send: sBlockEnd,
        });
        qi += qConsumed;
        si = strand > 0 ? sBlockEnd + 1 : sBlockEnd - 1;
      }
      return blocks;
    }

    // UniProt FASTA-header-style IDs look like "sp|P12345|NAME_ORG" or
    // "tr|P12345|NAME_ORG" -- pull the accession out so the title can link
    // straight to the UniProt entry (that's the fastest way to get the
    // protein's real name today, short of embedding gene_name/description
    // into the shard itself). Not every proteinId matches (e.g. test
    // fixtures, or a future non-UniProt-sourced protein set), so this is a
    // best-effort addition, not a required field.
    function uniprotAccession(proteinId) {
      var m = /^(?:sp|tr)\|([^|]+)\|/.exec(proteinId);
      return m ? m[1] : null;
    }

    // queryDesc: the shard entry's own {gene_name?, description?} (schema v2) --
    // the query protein's real name, not the target's (see lib/build_alignment_
    // shards.py's module docstring). Undefined/empty fields are simply omitted,
    // so a shard built before --descriptions was wired in (or a query with no
    // matching FASTA header) still renders exactly as before.
    function renderHit(dialogEls, proteinId, hit, queryDesc) {
      dialogEls.title.textContent = "";
      dialogEls.title.appendChild(document.createTextNode(proteinId));
      var acc = uniprotAccession(proteinId);
      if (acc) {
        dialogEls.title.appendChild(document.createTextNode(" "));
        var link = document.createElement("a");
        link.href = "https://www.uniprot.org/uniprotkb/" + encodeURIComponent(acc);
        link.target = "_blank";
        link.rel = "noopener noreferrer";
        link.className = "align-uniprot-link";
        link.textContent = "(UniProt ↗)";
        dialogEls.title.appendChild(link);
      }
      dialogEls.title.appendChild(document.createTextNode(" vs " + hit.sseqid + " (" + hit.genome + ")"));
      var alignedNote = hit.aligned_as ? (" (shown via cluster representative " + hit.aligned_as + ")") : "";
      var queryNote = "";
      if (queryDesc && (queryDesc.gene_name || queryDesc.description)) {
        queryNote = "Query= " +
          [queryDesc.gene_name, queryDesc.description].filter(Boolean).join(" ") + "  ";
      }
      dialogEls.stats.textContent =
        queryNote +
        "Subject= " + hit.sseqid + "  evalue=" + hit.evalue + "  bitscore=" + hit.bitscore +
        "  pident=" + hit.pident.toFixed(1) + "%" +
        "  length=" + hit.length + "  range=" + hit.sstart + "-" + hit.send +
        " (frame " + hit.sframe + ")" + alignedNote;
      dialogEls.body.textContent = "";
      var blocks = chunkAlignment(hit.qseq, hit.sseq, hit.qstart, hit.sstart, hit.sframe);
      // Pad every block's position numbers to the same width (the widest
      // number seen across the whole HSP), so Query/Sbjct columns line up
      // block to block the way NCBI BLAST's own text output does.
      var numWidth = Math.max(
        String(hit.qstart).length, String(hit.qend).length,
        String(hit.sstart).length, String(hit.send).length
      );
      blocks.forEach(function (b) {
        var pre = document.createElement("pre");
        pre.className = "align-block";
        var qLabel = "Query  " + String(b.qstart).padStart(numWidth, " ") + "  ";
        var sLabel = "Sbjct  " + String(b.sstart).padStart(numWidth, " ") + "  ";
        var labelWidth = Math.max(qLabel.length, sLabel.length);
        qLabel = qLabel.padEnd(labelWidth, " ");
        sLabel = sLabel.padEnd(labelWidth, " ");
        pre.textContent =
          qLabel + b.qseq + "  " + b.qend + "\n" +
          " ".repeat(labelWidth) + b.mid + "\n" +
          sLabel + b.sseq + "  " + b.send;
        dialogEls.body.appendChild(pre);
      });
    }

    function renderEmpty(dialogEls, message) {
      dialogEls.title.textContent = "";
      dialogEls.stats.textContent = "";
      dialogEls.body.textContent = "";
      var p = document.createElement("p");
      p.className = "align-empty";
      p.textContent = message;
      dialogEls.body.appendChild(p);
    }

    function dialogEls() {
      return {
        dialog: document.getElementById("alignment-dialog"),
        title: document.getElementById("alignment-title"),
        stats: document.getElementById("alignment-stats"),
        body: document.getElementById("alignment-body"),
      };
    }

    async function open(baseUrl, genome, proteinId, hitIndex) {
      var els = dialogEls();
      if (!els.dialog) return; // page didn't include ALIGNMENT_POPUP_HTML
      els.dialog.showModal();
      renderEmpty(els, "Loading alignment...");
      try {
        var shard = await fetchDataShard(baseUrl + genome + ".json.gz");
        // Schema v2: shard[proteinId] = {hits: [...], gene_name?, description?}.
        // Schema v1 (pre-existing shards not yet rebuilt by a real pipeline run since
        // the v2 upgrade -- the report HTML embedding this script can be regenerated
        // standalone by bin/sync_reports.sh well before the shard files themselves
        // are, since those only come from an actual BUILD_ALIGNMENT_SHARDS run):
        // shard[proteinId] is the bare hits array itself, with no query name/
        // description available. Detect which shape this is rather than assuming
        // v2 and silently reporting "no archived alignment" for real v1 data.
        var raw = shard[proteinId];
        var hits, queryDesc;
        if (Array.isArray(raw)) {
          hits = raw;
          queryDesc = null;
        } else {
          hits = raw && raw.hits;
          queryDesc = raw;
        }
        if (!hits || !hits.length) {
          renderEmpty(els, "No archived alignment for " + proteinId + " vs " + genome + ".");
          return;
        }
        renderHit(els, proteinId, hits[hitIndex || 0], queryDesc);
      } catch (err) {
        if (err && err.message === "DECOMPRESSION_UNSUPPORTED") {
          renderEmpty(els, "Alignment viewer requires a modern browser (Chrome/Firefox/Safari, last ~2 years).");
        } else {
          renderEmpty(els, "Could not load alignment data (" + (err && err.message ? err.message : err) + ").");
        }
      }
    }

    // Shared click-handler entry point (issue #86): a plain click opens the
    // in-page dialog; Ctrl/Cmd-click opens the standalone alignment.html
    // page in a new tab instead -- a real bookmarkable URL, same convention
    // as any ordinary link's modifier-click behaviour. alignment.html lives
    // alongside the calling report (docs/<project>/alignment.html), so a
    // bare relative href resolves correctly from either novelties.html or
    // losses.html.
    function openOrNewTab(baseUrl, genome, proteinId, ev, hitIndex) {
      if (ev && (ev.ctrlKey || ev.metaKey)) {
        var qs = "dir=" + encodeURIComponent(baseUrl) +
          "&genome=" + encodeURIComponent(genome) +
          "&protein=" + encodeURIComponent(proteinId) +
          (hitIndex ? "&hit=" + encodeURIComponent(hitIndex) : "");
        window.open("alignment.html?" + qs, "_blank", "noopener");
        return;
      }
      open(baseUrl, genome, proteinId, hitIndex);
    }

    function wireClose() {
      var els = dialogEls();
      if (!els.dialog) return;
      var closeBtn = document.getElementById("alignment-close");
      if (closeBtn) closeBtn.addEventListener("click", function () { els.dialog.close(); });
      els.dialog.addEventListener("click", function (ev) {
        if (ev.target === els.dialog) els.dialog.close(); // click on backdrop
      });
    }
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", wireClose);
    } else {
      wireClose();
    }

    window.NIAlignments = { open: open, openOrNewTab: openOrNewTab, fetchDataShard: fetchDataShard };
  })();
"""

# Small, synchronous popup for a single presence-chip's search-hit evidence
# (e-value, target protein) -- everything it shows is already embedded in the
# report payload, so unlike ALIGNMENT_POPUP_* above it needs no fetch/loading
# state and works in both the online (docs/) and offline (results/) report
# copies alike. Percent identity is deliberately not shown yet -- the pairwise
# search pathway (lib/hits.py's Hit/parse_tabular) doesn't capture it.
HIT_INFO_POPUP_CSS = r"""
  dialog.hitinfo {
    max-width: min(92vw, 440px);
    width: 100%;
    padding: 0;
    border: 1px solid var(--border);
    border-radius: 10px;
    background: var(--surface-1);
    color: var(--text-primary);
    box-shadow: var(--shadow);
  }
  dialog.hitinfo::backdrop { background: rgba(0, 0, 0, 0.45); }
  dialog.hitinfo .hi-head {
    display: flex; align-items: center; gap: 10px;
    padding: 14px 18px; border-bottom: 1px solid var(--border);
  }
  dialog.hitinfo .hi-head h3 { margin: 0; font-size: 14px; font-weight: 600; flex: 1; min-width: 0; }
  dialog.hitinfo .hi-body { padding: 14px 18px; }
  dialog.hitinfo .hi-row { margin-bottom: 10px; }
  dialog.hitinfo .hi-row:last-child { margin-bottom: 0; }
  dialog.hitinfo .hi-label {
    font-size: 11px; text-transform: uppercase; letter-spacing: 0.04em;
    color: var(--text-secondary); margin-bottom: 3px;
  }
  dialog.hitinfo .hi-value { font-size: 13px; word-break: break-word; }
"""

HIT_INFO_POPUP_HTML = r"""
  <dialog class="hitinfo" id="hitinfo-dialog">
    <div class="hi-head">
      <h3 id="hitinfo-title"></h3>
      <button type="button" class="btn-ghost" id="hitinfo-close" aria-label="Close">&times;</button>
    </div>
    <div class="hi-body" id="hitinfo-body"></div>
  </dialog>
"""

# Public surface: window.NIHitInfo.open(info) -- info: {proteomeLabel, present,
# contextNote, evalue, targetId, targetLabel}. Every field is optional except
# proteomeLabel/present; the popup only renders rows it actually has data for.
HIT_INFO_POPUP_JS = r"""
  (function () {
    // Self-contained, deliberately not sharing LINKOUT_HELPERS_JS's uniprotAcc/
    // uniprotLinkNode -- this fragment can be spliced into a page at a point
    // outside whatever scope those happen to be declared in (report_template.py
    // wraps its whole main script in one big IIFE; this fragment is appended
    // after it closes), so it carries its own tiny copy instead of risking a
    // ReferenceError depending on insertion order.
    function uniprotAcc(id) {
      var m = /(?:^|\|)([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:\||\s|$)/.exec(id || "");
      return m ? m[1] : "";
    }

    function dialogEls() {
      return {
        dialog: document.getElementById("hitinfo-dialog"),
        title: document.getElementById("hitinfo-title"),
        body: document.getElementById("hitinfo-body"),
      };
    }

    function row(body, label, valueNode) {
      var r = document.createElement("div");
      r.className = "hi-row";
      var l = document.createElement("div");
      l.className = "hi-label";
      l.textContent = label;
      var v = document.createElement("div");
      v.className = "hi-value";
      if (typeof valueNode === "string") v.textContent = valueNode;
      else v.appendChild(valueNode);
      r.appendChild(l);
      r.appendChild(v);
      body.appendChild(r);
    }

    function open(info) {
      var els = dialogEls();
      if (!els.dialog) return; // page didn't include HIT_INFO_POPUP_HTML
      els.title.textContent = info.proteomeLabel;
      els.body.textContent = "";
      row(els.body, "Status", info.present ? "Present" + (info.contextNote || "") : "Absent" + (info.contextNote || ""));
      if (info.evalue) row(els.body, "Hit e-value (" + info.searchLabel + " search)", info.evalue);
      if (info.targetId) {
        var acc = uniprotAcc(info.targetId);
        var label = (info.targetLabel && info.targetLabel !== info.targetId)
          ? info.targetLabel + " (" + info.targetId + ")" : info.targetId;
        var node;
        if (acc) {
          node = document.createElement("a");
          node.href = "https://www.uniprot.org/uniprotkb/" + acc + "/entry";
          node.target = "_blank";
          node.rel = "noopener noreferrer";
          node.textContent = label;
        } else {
          node = label;
        }
        row(els.body, "Target protein", node);
      }
      els.dialog.showModal();
    }

    function wireClose() {
      var els = dialogEls();
      if (!els.dialog) return;
      var closeBtn = document.getElementById("hitinfo-close");
      if (closeBtn) closeBtn.addEventListener("click", function () { els.dialog.close(); });
      els.dialog.addEventListener("click", function (ev) {
        if (ev.target === els.dialog) els.dialog.close(); // click on backdrop
      });
    }
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", wireClose);
    } else {
      wireClose();
    }

    window.NIHitInfo = { open: open };
  })();
"""
