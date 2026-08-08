"""
Small reusable HTML/CSS/JS string fragments for the lighter, single-table
report pages (lib/core_report_template.py, and the future LOSSES report).

lib/report_template.py — the canvas-heatmap novelty report — predates this
module and is not wired to it: it is large, thoroughly covered by
tests/test_report_data.py's self-contained-HTML checks, and restructuring it
to pull from here is not worth the risk for the handful of lines it would
save. New single-table templates should build on these fragments instead of
re-copying them inline, so the two report families stay visually and
behaviourally consistent without duplicating boilerplate.

Every fragment here is plain text meant to be concatenated into a larger
HTML_TEMPLATE-style triple-quoted string (see lib/report_template.py's
module docstring for the __PROJECT_TITLE__ / /*__PAYLOAD__*/ substitution
convention) — nothing here does its own token substitution.
"""

# Colour tokens shared with lib/report_template.py so every report reads as
# one visual system in light and dark mode. series-2 (TBLASTN green) is
# intentionally omitted — these single-table reports have nothing to encode
# a second evidence colour for.
THEME_VARS_CSS = r"""
  :root {
    color-scheme: light;
    --page: #f9f9f7;
    --surface-1: #fcfcfb;
    --text-primary: #0b0b0b;
    --text-secondary: #52514e;
    --muted: #898781;
    --grid: #e1e0d9;
    --axis: #c3c2b7;
    --border: rgba(11, 11, 11, 0.10);
    --series-1: #2a78d6;
    --wash: rgba(42, 120, 214, 0.10);
    --hover-wash: rgba(11, 11, 11, 0.05);
  }
  @media (prefers-color-scheme: dark) {
    :root:where(:not([data-theme="light"])) {
      color-scheme: dark;
      --page: #0d0d0d;
      --surface-1: #1a1a19;
      --text-primary: #ffffff;
      --text-secondary: #c3c2b7;
      --muted: #898781;
      --grid: #2c2c2a;
      --axis: #383835;
      --border: rgba(255, 255, 255, 0.10);
      --series-1: #3987e5;
      --wash: rgba(57, 135, 229, 0.16);
      --hover-wash: rgba(255, 255, 255, 0.06);
    }
  }
  :root[data-theme="dark"] {
    color-scheme: dark;
    --page: #0d0d0d;
    --surface-1: #1a1a19;
    --text-primary: #ffffff;
    --text-secondary: #c3c2b7;
    --muted: #898781;
    --grid: #2c2c2a;
    --axis: #383835;
    --border: rgba(255, 255, 255, 0.10);
    --series-1: #3987e5;
    --wash: rgba(57, 135, 229, 0.16);
    --hover-wash: rgba(255, 255, 255, 0.06);
  }
"""

# Page chrome, filter bar, single data table and detail panel — the shape
# every single-table report needs. No canvas/heatmap rules here; that stays
# specific to lib/report_template.py.
BASE_PAGE_CSS = r"""
  * { box-sizing: border-box; }
  body {
    margin: 0;
    background: var(--page);
    color: var(--text-primary);
    font: 14px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif;
  }
  .wrap { max-width: min(95vw, 2100px); margin: 0 auto; padding: 24px 20px 64px; }
  header.top { display: flex; align-items: flex-start; gap: 16px; margin-bottom: 24px; }
  header.top .titles { flex: 1; min-width: 0; }
  h1 { margin: 0 0 4px; font-size: 22px; font-weight: 600; letter-spacing: -0.01em; }
  .sub { margin: 0; color: var(--text-secondary); font-size: 13px; }
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
  .placeholder { color: var(--text-secondary); font-size: 13px; }
  .hidden { display: none !important; }
  .sr-only {
    position: absolute; width: 1px; height: 1px; padding: 0; margin: -1px;
    overflow: hidden; clip: rect(0 0 0 0); white-space: nowrap; border: 0;
  }
"""

# ---- JS fragments ---------------------------------------------------------

# Untrusted-string discipline: every report must build DOM nodes with this
# helper (or plain textContent), never innerHTML — see CLAUDE.md's report
# constraints. Protein IDs and annotation text come from FASTA headers and
# SwissProt/Pfam text and are not sanitised upstream.
EL_HELPER_JS = r"""
  function el(tag, cls, text) {
    var n = document.createElement(tag);
    if (cls) n.className = cls;
    if (text !== undefined && text !== null) n.textContent = text; // labels are untrusted data
    return n;
  }
"""

# Shared with lib/report_template.py's detail panel (kept independent there,
# see this module's docstring) — same gene-id stripping and SwissProt
# accession regex, so the two reports resolve the same protein to the same
# FungiDB/UniProt link.
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
  // Pfam names as a chip row (detail panel) -- accession + optional E-value in the
  // link title, same InterPro/Pfam target used by every report's "Pfam domains" field.
  function pfamChipsNode(pfamNames, pfamAccs, pfamEvs) {
    var chips = el("div", "chips");
    var names = pfamNames ? pfamNames.split(",") : [];
    var accs = pfamAccs ? pfamAccs.split(",") : [];
    var evs = pfamEvs ? pfamEvs.split(",") : [];
    names.forEach(function (n, i) {
      var acc = (accs[i] || "").split(".")[0];
      var ev = evs[i];
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
"""

THEME_TOGGLE_JS = r"""
  var themeBtn = document.getElementById("theme-toggle");
  themeBtn.addEventListener("click", function () {
    var cur = document.documentElement.getAttribute("data-theme");
    if (!cur) {
      cur = window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
    }
    document.documentElement.setAttribute("data-theme", cur === "dark" ? "light" : "dark");
  });
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
