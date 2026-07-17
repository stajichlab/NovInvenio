"""
Self-contained HTML/CSS/JS template for the NovInvenio interactive report.

bin/make_report.py substitutes two tokens:
  __PROJECT_TITLE__   plain-text project name (HTML-escaped by the caller)
  /*__PAYLOAD__*/     the JSON payload from lib/report_data.build_payload()

The page has no external dependencies — it must open from a file:// URL on a
laptop with no network, which is how these reports get shared off the cluster.

Colour roles come from the validated data-viz palette: series-1 (blue) marks
protein-search presence, series-2 (green) marks TBLASTN genome hits.  Column
position — not hue — carries ingroup vs outgroup, so no colour does double duty.
"""

HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__PROJECT_TITLE__ — NovInvenio candidates</title>
<style>
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
    --series-2: #008300;
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
      --series-2: #008300;
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
    --series-2: #008300;
    --wash: rgba(57, 135, 229, 0.16);
    --hover-wash: rgba(255, 255, 255, 0.06);
  }

  * { box-sizing: border-box; }
  body {
    margin: 0;
    background: var(--page);
    color: var(--text-primary);
    font: 14px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif;
  }
  .wrap { max-width: 1560px; margin: 0 auto; padding: 24px 20px 64px; }

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
  button:focus-visible, select:focus-visible, input:focus-visible, .grid-scroll:focus-visible {
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
  .card-title {
    margin: 0 0 2px;
    font-size: 14px;
    font-weight: 600;
  }
  .card-note { margin: 0 0 16px; color: var(--text-secondary); font-size: 12px; }

  /* ---- run summary ---- */
  .summary { display: grid; grid-template-columns: minmax(220px, 300px) 1fr; gap: 32px; }
  @media (max-width: 900px) { .summary { grid-template-columns: 1fr; } }

  .hero-label { color: var(--text-secondary); font-size: 13px; margin-bottom: 2px; }
  .hero-value { font-size: 52px; line-height: 1.05; font-weight: 600; letter-spacing: -0.02em; }
  .hero-note { color: var(--text-secondary); font-size: 12px; margin-top: 6px; }

  .tiles { display: flex; flex-wrap: wrap; gap: 28px; margin-top: 22px; }
  .tile-label { color: var(--text-secondary); font-size: 12px; }
  .tile-value { font-size: 20px; font-weight: 600; }

  /* ---- per-species bars ---- */
  .bars { display: grid; grid-template-columns: auto 1fr; gap: 6px 12px; align-items: center; }
  .bar-name { font-size: 12px; color: var(--text-secondary); white-space: nowrap; text-align: right; }
  .bar-track { display: flex; align-items: center; gap: 8px; }
  .bar-fill {
    height: 16px;
    background: var(--series-1);
    border-radius: 0 4px 4px 0;
    min-width: 2px;
  }
  .bar-value { font-size: 12px; font-variant-numeric: tabular-nums; color: var(--text-secondary); }

  /* ---- filters ---- */
  .filters { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; margin-bottom: 14px; }
  .filters label.check { display: inline-flex; align-items: center; gap: 6px; font-size: 13px; cursor: pointer; }
  .filters input[type="search"] { min-width: 260px; }
  .spacer { flex: 1; }
  .count { color: var(--text-secondary); font-size: 13px; font-variant-numeric: tabular-nums; }

  /* ---- tabs ---- */
  .tabs { display: flex; gap: 4px; margin-bottom: 12px; }
  .tab {
    border: 1px solid transparent;
    background: transparent;
    border-radius: 6px 6px 0 0;
    padding: 6px 14px;
    color: var(--text-secondary);
  }
  .tab[aria-selected="true"] {
    color: var(--text-primary);
    background: var(--surface-1);
    border-color: var(--border);
    border-bottom-color: var(--surface-1);
    font-weight: 600;
  }

  /* ---- explorer ---- */
  .explorer { display: grid; grid-template-columns: 1fr 380px; gap: 16px; align-items: start; }
  @media (max-width: 1180px) { .explorer { grid-template-columns: 1fr; } }

  .legend { display: flex; flex-wrap: wrap; gap: 18px; margin-bottom: 12px; }
  .legend-item { display: inline-flex; align-items: center; gap: 7px; font-size: 12px; color: var(--text-secondary); }
  .swatch { width: 12px; height: 12px; border-radius: 2px; flex: none; }
  .swatch.pres { background: var(--series-1); }
  .swatch.tb { background: var(--series-2); }
  .swatch.absent { background: var(--grid); }

  .grid-hscroll { overflow-x: auto; }
  .grid-scroll { overflow-y: auto; height: 560px; position: relative; }
  .grid-spacer { position: relative; width: 100%; }
  canvas.grid { display: block; position: absolute; left: 0; top: 0; }
  canvas.head { display: block; }
  .empty {
    padding: 40px 8px;
    color: var(--text-secondary);
    font-size: 13px;
  }

  /* ---- tooltip ---- */
  .tip {
    position: fixed;
    z-index: 20;
    max-width: 340px;
    background: var(--surface-1);
    border: 1px solid var(--border);
    border-radius: 8px;
    box-shadow: 0 6px 24px rgba(0, 0, 0, 0.18);
    padding: 10px 12px;
    font-size: 12px;
    pointer-events: none;
    display: none;
  }
  .tip-value { font-size: 13px; font-weight: 600; margin-bottom: 2px; }
  .tip-row { color: var(--text-secondary); display: flex; gap: 8px; align-items: baseline; }
  .tip-key { display: inline-block; width: 10px; height: 2px; border-radius: 1px; flex: none; transform: translateY(-3px); }
  .tip-id { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 11px; color: var(--text-secondary); margin-bottom: 6px; word-break: break-all; }

  /* ---- detail panel ---- */
  .detail { position: sticky; top: 16px; }
  .detail h3 { margin: 0 0 2px; font-size: 15px; font-weight: 600; word-break: break-all; }
  .detail .species { color: var(--text-secondary); font-size: 12px; font-style: italic; margin-bottom: 14px; }
  .field { margin-bottom: 12px; }
  .field-label { font-size: 11px; text-transform: uppercase; letter-spacing: 0.04em; color: var(--muted); margin-bottom: 3px; }
  .field-value { font-size: 13px; word-break: break-word; }
  .mono { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 11px; }
  .chips { display: flex; flex-wrap: wrap; gap: 5px; }
  .chip {
    display: inline-block;
    padding: 2px 7px;
    border: 1px solid var(--border);
    border-radius: 999px;
    font-size: 11px;
    text-decoration: none;
    color: var(--text-primary);
    background: var(--page);
  }
  .chip:hover { border-color: var(--series-1); }
  .links { display: flex; flex-wrap: wrap; gap: 6px; }
  .links a, .links button {
    font-size: 12px;
    padding: 5px 9px;
    border: 1px solid var(--border);
    border-radius: 6px;
    text-decoration: none;
    color: var(--text-primary);
    background: var(--page);
  }
  .links a:hover, .links button:hover { border-color: var(--series-1); }
  .seq {
    font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
    font-size: 10.5px;
    line-height: 1.45;
    word-break: break-all;
    max-height: 130px;
    overflow-y: auto;
    background: var(--page);
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 8px;
    color: var(--text-secondary);
  }
  .presence-mini { display: flex; flex-wrap: wrap; gap: 4px; }
  .pm {
    font-size: 11px;
    padding: 2px 6px;
    border-radius: 4px;
    border: 1px solid var(--border);
    font-variant-numeric: tabular-nums;
  }
  .pm.on-pres { background: var(--series-1); color: #fff; border-color: transparent; }
  .pm.on-tb { background: var(--series-2); color: #fff; border-color: transparent; }
  .pm.off { color: var(--muted); }
  .placeholder { color: var(--text-secondary); font-size: 13px; }

  /* ---- table view ---- */
  .tbl-scroll { overflow: auto; height: 560px; }
  table.data { border-collapse: collapse; width: 100%; font-size: 12px; }
  table.data th, table.data td {
    text-align: left;
    padding: 6px 10px;
    border-bottom: 1px solid var(--grid);
    white-space: nowrap;
  }
  table.data th {
    position: sticky;
    top: 0;
    background: var(--surface-1);
    font-weight: 600;
    z-index: 1;
    border-bottom: 1px solid var(--axis);
  }
  table.data td.num { font-variant-numeric: tabular-nums; }
  table.data td.wrap-cell { white-space: normal; max-width: 320px; }
  table.data tbody tr:hover { background: var(--hover-wash); }
  .hidden { display: none !important; }
  .sr-only {
    position: absolute; width: 1px; height: 1px; padding: 0; margin: -1px;
    overflow: hidden; clip: rect(0 0 0 0); white-space: nowrap; border: 0;
  }
</style>
</head>
<body>
<div class="wrap">
  <header class="top">
    <div class="titles">
      <h1 id="title"></h1>
      <p class="sub" id="subtitle"></p>
    </div>
    <button id="theme-toggle" class="btn-ghost" type="button" aria-label="Toggle colour theme">Theme</button>
  </header>

  <!-- Run summary describes the whole run; the filter row below scopes only the explorer. -->
  <section class="card">
    <h2 class="card-title">Run summary</h2>
    <p class="card-note">Counts cover every protein in the presence matrix, independent of the filters below.</p>
    <div class="summary">
      <div>
        <div class="hero-label">Novelty candidates</div>
        <div class="hero-value" id="hero"></div>
        <div class="hero-note" id="hero-note"></div>
        <div class="tiles">
          <div>
            <div class="tile-label">Proteins in matrix</div>
            <div class="tile-value" id="t-total"></div>
          </div>
          <div>
            <div class="tile-label">Ingroup</div>
            <div class="tile-value" id="t-in"></div>
          </div>
          <div>
            <div class="tile-label">Outgroup</div>
            <div class="tile-value" id="t-out"></div>
          </div>
          <div>
            <div class="tile-label">Functionally annotated</div>
            <div class="tile-value" id="t-annot"></div>
          </div>
        </div>
      </div>
      <div>
        <h3 class="card-title">Novelty candidates per ingroup species</h3>
        <p class="card-note">Proteins originating in each ingroup species that pass the novelty filter.</p>
        <div class="bars" id="species-bars"></div>
      </div>
    </div>
  </section>

  <div class="filters" role="group" aria-label="Filter candidates">
    <input type="search" id="f-search" placeholder="Search ID, gene, product, Pfam…" aria-label="Search proteins">
    <select id="f-src" aria-label="Source proteome"></select>
    <select id="f-fsrc" aria-label="Annotation source"></select>
    <select id="f-family" aria-label="Gene family"></select>
    <label class="check"><input type="checkbox" id="f-nov" checked> Novelty candidates only</label>
    <label class="check"><input type="checkbox" id="f-pfam"> Has Pfam</label>
    <label class="check"><input type="checkbox" id="f-notb"> No TBLASTN hit</label>
    <select id="f-sort" aria-label="Sort by">
      <option value="ingroup">Sort: ingroup breadth</option>
      <option value="id">Sort: protein ID</option>
      <option value="src">Sort: source proteome</option>
      <option value="outgroup">Sort: fewest outgroup hits</option>
      <option value="tb">Sort: fewest TBLASTN hits</option>
      <option value="pfam">Sort: annotated first</option>
    </select>
    <button id="f-reset" type="button">Reset</button>
    <div class="spacer"></div>
    <span class="count" id="count"></span>
    <button id="dl-tsv" type="button">Download TSV</button>
    <button id="dl-fa" type="button">Download FASTA</button>
  </div>

  <div class="tabs" role="tablist">
    <button class="tab" id="tab-heatmap" role="tab" aria-selected="true" aria-controls="panel-heatmap">Heatmap</button>
    <button class="tab" id="tab-table" role="tab" aria-selected="false" aria-controls="panel-table">Table</button>
  </div>

  <div class="explorer">
    <section class="card" style="padding:14px">
      <div id="panel-heatmap" role="tabpanel" aria-labelledby="tab-heatmap">
        <div class="legend">
          <span class="legend-item"><span class="swatch pres"></span>Present (protein search)</span>
          <span class="legend-item"><span class="swatch tb"></span>TBLASTN hit (outgroup genome)</span>
          <span class="legend-item"><span class="swatch absent"></span>Absent</span>
        </div>
        <div class="grid-hscroll">
          <canvas class="head" id="head"></canvas>
          <div class="grid-scroll" id="scroll" tabindex="0" role="application"
               aria-label="Presence heatmap. Use arrow keys to move between proteins.">
            <div class="grid-spacer" id="spacer"><canvas class="grid" id="grid"></canvas></div>
          </div>
        </div>
        <p class="empty hidden" id="grid-empty">No proteins match the current filters.</p>
      </div>

      <div id="panel-table" role="tabpanel" aria-labelledby="tab-table" class="hidden">
        <div class="tbl-scroll" id="tbl-scroll">
          <table class="data">
            <caption class="sr-only">Presence and annotation for every protein matching the current filters</caption>
            <thead id="tbl-head"></thead>
            <tbody id="tbl-body"></tbody>
          </table>
        </div>
        <p class="empty hidden" id="tbl-empty">No proteins match the current filters.</p>
      </div>
    </section>

    <aside class="card detail" id="detail" aria-live="polite">
      <p class="placeholder">Select a protein in the heatmap or table to see its annotation and database links.</p>
    </aside>
  </div>
</div>

<div class="tip" id="tip" role="tooltip"></div>

<script type="application/json" id="payload">/*__PAYLOAD__*/</script>
<script>
(function () {
  "use strict";

  var DATA = JSON.parse(document.getElementById("payload").textContent);
  var F = {};
  DATA.fields.forEach(function (name, i) { F[name] = i; });

  var ROWS = DATA.rows;
  var PROTEOMES = DATA.proteomes;
  var TB_GENOMES = DATA.tblastn_genomes;
  var FAMILIES = DATA.families || [];
  var N_IN = PROTEOMES.filter(function (p) { return p.group === "IN"; }).length;

  // Gene families group candidates recovered independently in several ingroup
  // species (via mmseqs clustering of the candidate set) — see lib/clusters.py.
  // Membership is carried per-row via F.fam; family-level stats (size, species)
  // come precomputed from the payload rather than scanning ROWS at render time.
  function familyLabel(fam) {
    return fam.rep + " (" + fam.size + " in " + fam.species.length +
      (fam.species.length === 1 ? " species)" : " species)");
  }

  // ---- column model -------------------------------------------------------
  // Column position carries ingroup/outgroup/TBLASTN; hue only carries the kind
  // of evidence (protein search vs TBLASTN).
  var COLS = [];
  PROTEOMES.forEach(function (p, i) {
    COLS.push({ kind: "pres", idx: i, label: p.short, group: p.group, proteome: p });
  });
  TB_GENOMES.forEach(function (g, i) {
    var p = PROTEOMES.find(function (x) { return x.short === g; });
    COLS.push({ kind: "tb", idx: i, label: g, group: "TB", proteome: p });
  });

  // Band labels sit over a block only as wide as its columns, so each carries a
  // short form; drawHead picks the longest that fits rather than clipping.
  var GROUP_LABEL = {
    IN: ["Ingroup proteomes", "Ingroup"],
    OUT: ["Outgroup proteomes", "Outgroup"],
    TB: ["Outgroup genomes (TBLASTN)", "Genomes (TBLASTN)", "TBLASTN"]
  };

  // Contiguous runs of the same group, for the header bands and separators.
  var BLOCKS = [];
  COLS.forEach(function (c, i) {
    var last = BLOCKS[BLOCKS.length - 1];
    if (last && last.group === c.group) { last.end = i; }
    else { BLOCKS.push({ group: c.group, start: i, end: i }); }
  });

  // ---- derived per-row stats ---------------------------------------------
  var nRows = ROWS.length;
  var inN = new Int16Array(nRows);
  var outN = new Int16Array(nRows);
  var tbN = new Int16Array(nRows);
  var hasPfam = new Uint8Array(nRows);

  for (var r = 0; r < nRows; r++) {
    var pres = ROWS[r][F.pres];
    var ci = 0, cn = 0;
    for (ci = 0; ci < N_IN; ci++) { if (pres.charCodeAt(ci) === 49) cn++; }
    inN[r] = cn;
    cn = 0;
    for (ci = N_IN; ci < pres.length; ci++) { if (pres.charCodeAt(ci) === 49) cn++; }
    outN[r] = cn;
    var tb = ROWS[r][F.tb];
    cn = 0;
    for (ci = 0; ci < tb.length; ci++) { if (tb.charCodeAt(ci) === 49) cn++; }
    tbN[r] = cn;
    hasPfam[r] = ROWS[r][F.pfam_n] ? 1 : 0;
  }

  // Lowercased haystack per row, built once — search runs on every keystroke.
  // Includes the row's family representative ID, so searching for any one
  // member's family rep surfaces every other ingroup species' copy too.
  var HAY = new Array(nRows);
  for (var h = 0; h < nRows; h++) {
    var row = ROWS[h];
    var famRep = row[F.fam] >= 0 ? FAMILIES[row[F.fam]].rep : "";
    HAY[h] = (row[F.id] + " " + row[F.gene] + " " + row[F.prod] + " " +
              row[F.pfam_n] + " " + row[F.sprot] + " " + famRep).toLowerCase();
  }

  // ---- state --------------------------------------------------------------
  var state = {
    search: "",
    src: "",
    fsrc: "",
    family: -1,
    novOnly: true,
    pfamOnly: false,
    noTb: false,
    sort: "ingroup",
    view: "heatmap",
    selected: -1,
    hoverRow: -1,
    hoverCol: -1
  };
  var view = [];

  // ---- geometry -----------------------------------------------------------
  var ROW_H = 20;
  var CELL_W = 46;
  var GUTTER = 260;
  var BLOCK_GAP = 10; // extra air between column blocks
  var LABEL_ANGLE = Math.PI / 3;
  var LABEL_CAP = 130;  // longest column label we grow the header for
  // Both recomputed by measureHead() — Short IDs vary in length between configs
  // (Ncra vs AOBS17_30), and a fixed header clips the long ones into ambiguity.
  var HEAD_H = 78;
  var LABEL_MAX_W = 50;

  function colX(i) {
    var x = GUTTER;
    for (var b = 0; b < BLOCKS.length; b++) {
      var blk = BLOCKS[b];
      if (i >= blk.start && i <= blk.end) return x + (i - blk.start) * CELL_W;
      x += (blk.end - blk.start + 1) * CELL_W + BLOCK_GAP;
    }
    return x;
  }
  var TOTAL_W = colX(COLS.length - 1) + CELL_W + 8;

  function css(name) {
    return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  }
  function palette() {
    return {
      surface: css("--surface-1"),
      grid: css("--grid"),
      axis: css("--axis"),
      primary: css("--text-primary"),
      secondary: css("--text-secondary"),
      muted: css("--muted"),
      s1: css("--series-1"),
      s2: css("--series-2"),
      wash: css("--wash"),
      hover: css("--hover-wash")
    };
  }

  // ---- filtering & sorting -----------------------------------------------
  function applyFilters() {
    var q = state.search.trim().toLowerCase();
    var terms = q ? q.split(/\s+/) : [];
    var srcIdx = state.src ? PROTEOMES.findIndex(function (p) { return p.short === state.src; }) : -1;
    var fsrcIdx = state.fsrc ? DATA.fsources.indexOf(state.fsrc) : -2;

    view = [];
    for (var i = 0; i < nRows; i++) {
      var row = ROWS[i];
      // A family filter is a request to see every member regardless of that
      // member's own novelty status, so it overrides (not stacks with) novOnly.
      if (state.family >= 0) {
        if (row[F.fam] !== state.family) continue;
      } else if (state.novOnly && !row[F.nov]) continue;
      if (srcIdx >= 0 && row[F.src] !== srcIdx) continue;
      if (state.fsrc && row[F.fsrc] !== fsrcIdx) continue;
      if (state.pfamOnly && !hasPfam[i]) continue;
      if (state.noTb && tbN[i] > 0) continue;
      if (terms.length) {
        var hay = HAY[i], ok = true;
        for (var t = 0; t < terms.length; t++) {
          if (hay.indexOf(terms[t]) === -1) { ok = false; break; }
        }
        if (!ok) continue;
      }
      view.push(i);
    }
    sortView();
  }

  function sortView() {
    var s = state.sort;
    var cmpId = function (a, b) { return ROWS[a][F.id] < ROWS[b][F.id] ? -1 : ROWS[a][F.id] > ROWS[b][F.id] ? 1 : 0; };
    var cmp;
    if (s === "id") cmp = cmpId;
    else if (s === "src") cmp = function (a, b) { return (ROWS[a][F.src] - ROWS[b][F.src]) || cmpId(a, b); };
    else if (s === "outgroup") cmp = function (a, b) { return (outN[a] - outN[b]) || (inN[b] - inN[a]) || cmpId(a, b); };
    else if (s === "tb") cmp = function (a, b) { return (tbN[a] - tbN[b]) || (inN[b] - inN[a]) || cmpId(a, b); };
    else if (s === "pfam") cmp = function (a, b) { return (hasPfam[b] - hasPfam[a]) || (inN[b] - inN[a]) || cmpId(a, b); };
    else cmp = function (a, b) { return (inN[b] - inN[a]) || (outN[a] - outN[b]) || cmpId(a, b); };
    view.sort(cmp);
  }

  // ---- heatmap ------------------------------------------------------------
  var scrollEl = document.getElementById("scroll");
  var spacerEl = document.getElementById("spacer");
  var gridCanvas = document.getElementById("grid");
  var headCanvas = document.getElementById("head");
  var gctx = gridCanvas.getContext("2d");
  var hctx = headCanvas.getContext("2d");

  function sizeCanvas(canvas, ctx, w, h) {
    var dpr = window.devicePixelRatio || 1;
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
    canvas.style.width = w + "px";
    canvas.style.height = h + "px";
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  function ellipsize(ctx, text, maxW) {
    if (ctx.measureText(text).width <= maxW) return text;
    var lo = 0, hi = text.length;
    while (lo < hi) {
      var mid = (lo + hi) >> 1;
      if (ctx.measureText(text.slice(0, mid) + "…").width <= maxW) lo = mid + 1;
      else hi = mid;
    }
    return text.slice(0, Math.max(0, lo - 1)) + "…";
  }

  var COL_LABEL_FONT = "600 11px system-ui, -apple-system, 'Segoe UI', sans-serif";

  // Grow the header to fit the longest rotated label, so column identity is never
  // lost to an ellipsis. Must run before sizeCanvas, which resets canvas state.
  function measureHead() {
    hctx.setTransform(1, 0, 0, 1, 0, 0);
    hctx.font = COL_LABEL_FONT;
    var maxW = 0;
    COLS.forEach(function (c) { maxW = Math.max(maxW, hctx.measureText(c.label).width); });
    LABEL_MAX_W = Math.min(maxW, LABEL_CAP);
    HEAD_H = Math.max(78, Math.round(34 + LABEL_MAX_W * Math.sin(LABEL_ANGLE)));
  }

  function drawHead() {
    var P = palette();
    measureHead();
    sizeCanvas(headCanvas, hctx, TOTAL_W, HEAD_H);
    hctx.clearRect(0, 0, TOTAL_W, HEAD_H);
    hctx.fillStyle = P.surface;
    hctx.fillRect(0, 0, TOTAL_W, HEAD_H);

    // group band labels
    hctx.font = "600 11px system-ui, -apple-system, 'Segoe UI', sans-serif";
    hctx.textBaseline = "middle";
    BLOCKS.forEach(function (blk) {
      var x0 = colX(blk.start);
      var x1 = colX(blk.end) + CELL_W;
      hctx.fillStyle = P.muted;
      var avail = x1 - x0 - 4;
      var options = GROUP_LABEL[blk.group] || [blk.group];
      var label = options.find(function (o) { return hctx.measureText(o).width <= avail; });
      if (label === undefined) label = ellipsize(hctx, options[options.length - 1], avail);
      hctx.textAlign = "left";
      hctx.fillText(label, x0, 12);
      hctx.fillStyle = P.axis;
      hctx.fillRect(x0, 22, x1 - x0, 1);
    });

    // rotated per-column labels
    hctx.textAlign = "left";
    hctx.textBaseline = "middle";
    COLS.forEach(function (c, i) {
      var x = colX(i) + CELL_W / 2;
      hctx.save();
      hctx.translate(x, HEAD_H - 8);
      hctx.rotate(-LABEL_ANGLE);
      hctx.fillStyle = P.primary;
      hctx.font = COL_LABEL_FONT;
      hctx.fillText(ellipsize(hctx, c.label, LABEL_MAX_W), 0, 0);
      hctx.restore();
    });

    hctx.textAlign = "left";
    hctx.fillStyle = P.muted;
    hctx.font = "11px system-ui, -apple-system, 'Segoe UI', sans-serif";
    hctx.fillText("Protein", 8, HEAD_H - 14);
  }

  function drawGrid() {
    var P = palette();
    var vh = scrollEl.clientHeight;
    var top = scrollEl.scrollTop;
    sizeCanvas(gridCanvas, gctx, TOTAL_W, vh);
    gridCanvas.style.top = top + "px";

    gctx.clearRect(0, 0, TOTAL_W, vh);
    gctx.fillStyle = P.surface;
    gctx.fillRect(0, 0, TOTAL_W, vh);

    var first = Math.max(0, Math.floor(top / ROW_H));
    var last = Math.min(view.length - 1, Math.ceil((top + vh) / ROW_H));

    for (var vi = first; vi <= last; vi++) {
      var ri = view[vi];
      if (ri === undefined) continue;
      var row = ROWS[ri];
      var y = vi * ROW_H - top;

      var isSel = ri === state.selected;
      var isHover = ri === state.hoverRow;
      if (isSel) {
        gctx.fillStyle = P.wash;
        gctx.fillRect(0, y, TOTAL_W, ROW_H);
        gctx.fillStyle = P.s1;
        gctx.fillRect(0, y + 1, 3, ROW_H - 2);
      } else if (isHover) {
        gctx.fillStyle = P.hover;
        gctx.fillRect(0, y, TOTAL_W, ROW_H);
      }

      // gutter: protein ID, then gene/product in muted ink if there is room
      gctx.textBaseline = "middle";
      gctx.textAlign = "left";
      gctx.font = "11px ui-monospace, SFMono-Regular, Menlo, monospace";
      gctx.fillStyle = P.primary;
      var idText = ellipsize(gctx, row[F.id], GUTTER - 100);
      gctx.fillText(idText, 8, y + ROW_H / 2);

      var note = row[F.gene] || row[F.prod] || "";
      if (note) {
        gctx.font = "10px system-ui, -apple-system, 'Segoe UI', sans-serif";
        gctx.fillStyle = P.muted;
        gctx.textAlign = "right";
        gctx.fillText(ellipsize(gctx, note, 84), GUTTER - 10, y + ROW_H / 2);
        gctx.textAlign = "left";
      }

      // cells — 1px inset each side gives the 2px surface gap between fills
      for (var ci = 0; ci < COLS.length; ci++) {
        var c = COLS[ci];
        var on = c.kind === "pres"
          ? row[F.pres].charCodeAt(c.idx) === 49
          : row[F.tb].charCodeAt(c.idx) === 49;
        gctx.fillStyle = on ? (c.kind === "pres" ? P.s1 : P.s2) : P.grid;
        gctx.fillRect(colX(ci) + 1, y + 1, CELL_W - 2, ROW_H - 2);
      }

      if (isHover && state.hoverCol >= 0) {
        gctx.strokeStyle = P.primary;
        gctx.lineWidth = 1.5;
        gctx.strokeRect(colX(state.hoverCol) + 1.25, y + 1.25, CELL_W - 2.5, ROW_H - 2.5);
      }
    }
  }

  function rowAt(clientY) {
    var rect = scrollEl.getBoundingClientRect();
    var vi = Math.floor((clientY - rect.top + scrollEl.scrollTop) / ROW_H);
    return (vi >= 0 && vi < view.length) ? vi : -1;
  }
  function colAt(clientX) {
    var rect = gridCanvas.getBoundingClientRect();
    var x = clientX - rect.left;
    for (var i = 0; i < COLS.length; i++) {
      var cx = colX(i);
      if (x >= cx && x < cx + CELL_W) return i;
    }
    return -1;
  }

  // ---- tooltip ------------------------------------------------------------
  var tipEl = document.getElementById("tip");

  function el(tag, cls, text) {
    var n = document.createElement(tag);
    if (cls) n.className = cls;
    if (text !== undefined && text !== null) n.textContent = text; // labels are untrusted data
    return n;
  }

  function showTip(ri, ci, x, y) {
    var row = ROWS[ri];
    tipEl.textContent = "";

    var col = ci >= 0 ? COLS[ci] : null;
    if (col) {
      var on = col.kind === "pres"
        ? row[F.pres].charCodeAt(col.idx) === 49
        : row[F.tb].charCodeAt(col.idx) === 49;
      var verb = col.kind === "pres"
        ? (on ? "Present" : "Absent")
        : (on ? "TBLASTN hit" : "No TBLASTN hit");
      tipEl.appendChild(el("div", "tip-value", verb));       // value leads
      var sp = col.proteome;
      var name = sp ? (sp.species + (sp.strain ? " " + sp.strain : "")) : col.label;
      var r1 = el("div", "tip-row");
      var key = el("span", "tip-key");
      key.style.background = on ? (col.kind === "pres" ? "var(--series-1)" : "var(--series-2)") : "var(--grid)";
      r1.appendChild(key);
      r1.appendChild(el("span", null, name));
      tipEl.appendChild(r1);
      tipEl.appendChild(el("div", "tip-row", col.kind === "pres"
        ? (col.group === "IN" ? "Ingroup proteome" : "Outgroup proteome")
        : "Outgroup genome"));
    }

    tipEl.appendChild(el("div", "tip-id", row[F.id]));
    if (row[F.gene]) tipEl.appendChild(el("div", "tip-row", "Gene: " + row[F.gene]));
    if (row[F.prod]) tipEl.appendChild(el("div", "tip-row", row[F.prod]));
    tipEl.appendChild(el("div", "tip-row",
      "Ingroup " + inN[ri] + "/" + N_IN + " · outgroup " + outN[ri] + "/" + (PROTEOMES.length - N_IN) +
      (TB_GENOMES.length ? " · TBLASTN " + tbN[ri] + "/" + TB_GENOMES.length : "")));
    if (row[F.nov]) tipEl.appendChild(el("div", "tip-row", "Novelty candidate"));

    tipEl.style.display = "block";
    var w = tipEl.offsetWidth, hgt = tipEl.offsetHeight;
    var left = x + 14, topPos = y + 14;
    if (left + w > window.innerWidth - 8) left = x - w - 14;
    if (topPos + hgt > window.innerHeight - 8) topPos = y - hgt - 14;
    tipEl.style.left = Math.max(8, left) + "px";
    tipEl.style.top = Math.max(8, topPos) + "px";
  }
  function hideTip() { tipEl.style.display = "none"; }

  // ---- detail panel -------------------------------------------------------
  var detailEl = document.getElementById("detail");

  function geneIdFromProteinId(pid) {
    return pid.replace(/-[Tt][^-]*(-p\d+)?$|-p\d+$/, "");
  }
  function uniprotAcc(sprot) {
    if (!sprot) return "";
    var m = /(?:^|\|)([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:\||\s|$)/.exec(sprot);
    return m ? m[1] : "";
  }
  function field(label, valueNode) {
    var f = el("div", "field");
    f.appendChild(el("div", "field-label", label));
    if (typeof valueNode === "string") f.appendChild(el("div", "field-value", valueNode));
    else f.appendChild(valueNode);
    return f;
  }

  function renderDetail() {
    detailEl.textContent = "";
    var ri = state.selected;
    if (ri < 0) {
      detailEl.appendChild(el("p", "placeholder",
        "Select a protein in the heatmap or table to see its annotation and database links."));
      return;
    }
    var row = ROWS[ri];
    var sp = row[F.src] >= 0 ? PROTEOMES[row[F.src]] : null;

    detailEl.appendChild(el("h3", null, row[F.id]));
    if (sp) {
      detailEl.appendChild(el("div", "species",
        sp.species + (sp.strain ? " " + sp.strain : "") + " · " + sp.short +
        " · " + (sp.group === "IN" ? "ingroup" : "outgroup")));
    }

    if (row[F.nov]) {
      detailEl.appendChild(field("Status", "Novelty candidate — present in ≥ " +
        Math.round(DATA.ingroup_min_frac * 100) + "% of the ingroup and absent from every outgroup proteome."));
    }

    if (row[F.fam] >= 0) {
      var fam = FAMILIES[row[F.fam]];
      var famBox = el("div");
      famBox.appendChild(el("div", "field-value",
        fam.size + " members across " + fam.species.length + " species (" + fam.species.join(", ") + ")"));
      var famLinks = el("div", "links");
      var famBtn = el("button", null, "Show family members (" + fam.size + ")");
      famBtn.type = "button";
      famBtn.addEventListener("click", function () { setFamilyFilter(row[F.fam]); setView("table"); });
      famLinks.appendChild(famBtn);
      famBox.appendChild(famLinks);
      detailEl.appendChild(field("Gene family — independently recovered in multiple ingroup species", famBox));
    }

    // presence chips: colour + the short ID text, so identity is never colour-alone
    var mini = el("div", "presence-mini");
    PROTEOMES.forEach(function (p, i) {
      var on = row[F.pres].charCodeAt(i) === 49;
      var chip = el("span", "pm " + (on ? "on-pres" : "off"), p.short);
      chip.title = p.species + (p.strain ? " " + p.strain : "") + " — " + (on ? "present" : "absent");
      mini.appendChild(chip);
    });
    detailEl.appendChild(field("Presence (protein search) · ingroup " + inN[ri] + "/" + N_IN +
      ", outgroup " + outN[ri] + "/" + (PROTEOMES.length - N_IN), mini));

    if (TB_GENOMES.length) {
      var tbm = el("div", "presence-mini");
      TB_GENOMES.forEach(function (g, i) {
        var on = row[F.tb].charCodeAt(i) === 49;
        tbm.appendChild(el("span", "pm " + (on ? "on-tb" : "off"), g));
      });
      detailEl.appendChild(field("TBLASTN vs outgroup genomes · " + tbN[ri] + "/" + TB_GENOMES.length + " hit", tbm));
    }

    if (row[F.gene]) detailEl.appendChild(field("Gene name", row[F.gene]));
    if (row[F.prod]) detailEl.appendChild(field("Product", row[F.prod]));
    if (row[F.fsrc] >= 0) detailEl.appendChild(field("Annotation source", DATA.fsources[row[F.fsrc]]));
    if (row[F.sprot]) detailEl.appendChild(field("Best SwissProt hit", row[F.sprot]));

    if (row[F.pfam_n]) {
      var names = row[F.pfam_n].split(",");
      var accs = row[F.pfam_a] ? row[F.pfam_a].split(",") : [];
      var evs = row[F.pfam_e] ? row[F.pfam_e].split(",") : [];
      var chips = el("div", "chips");
      names.forEach(function (n, i) {
        var acc = (accs[i] || "").split(".")[0];
        var ev = evs[i];
        var node;
        if (/^PF\d+$/.test(acc)) {
          node = el("a", "chip", n + (ev ? " · " + ev : ""));
          node.href = "https://www.ebi.ac.uk/interpro/entry/pfam/" + acc + "/";
          node.target = "_blank";
          node.rel = "noopener noreferrer";
          node.title = acc + (ev ? " — hmmsearch E-value " + ev : "");
        } else {
          node = el("span", "chip", n);
        }
        chips.appendChild(node);
      });
      detailEl.appendChild(field("Pfam domains (" + names.length + ")", chips));
    }

    // external links
    var links = el("div", "links");
    var acc = uniprotAcc(row[F.sprot]);
    if (acc) {
      var up = el("a", null, "UniProt " + acc);
      up.href = "https://www.uniprot.org/uniprotkb/" + acc + "/entry";
      up.target = "_blank"; up.rel = "noopener noreferrer";
      links.appendChild(up);
      var af = el("a", null, "AlphaFold");
      af.href = "https://alphafold.ebi.ac.uk/entry/" + acc;
      af.target = "_blank"; af.rel = "noopener noreferrer";
      links.appendChild(af);
    }
    var fsrcName = row[F.fsrc] >= 0 ? DATA.fsources[row[F.fsrc]] : "";
    if (fsrcName.indexOf("ModelOrg_") === 0) {
      var fd = el("a", null, "FungiDB gene");
      fd.href = "https://fungidb.org/fungidb/app/record/gene/" + encodeURIComponent(geneIdFromProteinId(row[F.id]));
      fd.target = "_blank"; fd.rel = "noopener noreferrer";
      links.appendChild(fd);
    }
    var seq = row[F.seq];
    if (seq) {
      var bl = el("a", null, "BLASTP at NCBI");
      bl.href = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&PROGRAM=blastp" +
                "&BLAST_PROGRAMS=blastp&DATABASE=nr&CMD=Web&QUERY=" + encodeURIComponent(seq);
      bl.target = "_blank"; bl.rel = "noopener noreferrer";
      links.appendChild(bl);

      var ipr = el("a", null, "InterProScan");
      ipr.href = "https://www.ebi.ac.uk/interpro/search/sequence/";
      ipr.target = "_blank"; ipr.rel = "noopener noreferrer";
      ipr.title = "Opens InterProScan — paste the sequence copied below";
      links.appendChild(ipr);

      var cp = el("button", null, "Copy FASTA");
      cp.type = "button";
      cp.addEventListener("click", function () {
        var fa = ">" + row[F.id] + "\n" + (seq.match(/.{1,60}/g) || []).join("\n") + "\n";
        navigator.clipboard.writeText(fa).then(function () {
          cp.textContent = "Copied";
          setTimeout(function () { cp.textContent = "Copy FASTA"; }, 1400);
        });
      });
      links.appendChild(cp);
    }
    if (links.childNodes.length) detailEl.appendChild(field("External resources", links));

    if (seq) {
      detailEl.appendChild(field("Protein sequence (" + seq.length + " aa)", el("div", "seq", seq)));
    }
  }

  function select(ri) {
    state.selected = ri;
    renderDetail();
    drawGrid();
    if (state.view === "table") markTableSelection();
  }

  // ---- table view ---------------------------------------------------------
  var TBL_COLS = [
    { label: "Protein ID", get: function (r) { return ROWS[r][F.id]; }, cls: "mono" },
    { label: "Source", get: function (r) { return ROWS[r][F.src] >= 0 ? PROTEOMES[ROWS[r][F.src]].short : ""; } },
    { label: "Novelty", get: function (r) { return ROWS[r][F.nov] ? "yes" : "no"; } },
    {
      label: "Gene family", cls: "wrap-cell",
      get: function (r) {
        var fi = ROWS[r][F.fam];
        return fi >= 0 ? FAMILIES[fi].rep + " (" + FAMILIES[fi].size + ")" : "";
      }
    },
    { label: "Ingroup", get: function (r) { return inN[r] + "/" + N_IN; }, cls: "num" },
    { label: "Outgroup", get: function (r) { return outN[r] + "/" + (PROTEOMES.length - N_IN); }, cls: "num" },
    { label: "TBLASTN", get: function (r) { return TB_GENOMES.length ? tbN[r] + "/" + TB_GENOMES.length : "—"; }, cls: "num" },
    { label: "Gene", get: function (r) { return ROWS[r][F.gene]; } },
    { label: "Product", get: function (r) { return ROWS[r][F.prod]; }, cls: "wrap-cell" },
    { label: "Source of annotation", get: function (r) { return ROWS[r][F.fsrc] >= 0 ? DATA.fsources[ROWS[r][F.fsrc]] : ""; } },
    { label: "Pfam domains", get: function (r) { return ROWS[r][F.pfam_n]; }, cls: "wrap-cell" }
  ];
  // Per-proteome presence columns keep the table a true twin of the heatmap.
  PROTEOMES.forEach(function (p, i) {
    TBL_COLS.push({
      label: p.short + (p.group === "IN" ? " (in)" : " (out)"),
      get: function (r) { return ROWS[r][F.pres].charCodeAt(i) === 49 ? "1" : "0"; },
      cls: "num"
    });
  });
  TB_GENOMES.forEach(function (g, i) {
    TBL_COLS.push({
      label: g + " (tblastn)",
      get: function (r) { return ROWS[r][F.tb].charCodeAt(i) === 49 ? "1" : "0"; },
      cls: "num"
    });
  });

  var TBL_PAGE = 300;
  var tblShown = TBL_PAGE;
  var tblBody = document.getElementById("tbl-body");
  var tblScroll = document.getElementById("tbl-scroll");

  function renderTableHead() {
    var thead = document.getElementById("tbl-head");
    thead.textContent = "";
    var tr = document.createElement("tr");
    TBL_COLS.forEach(function (c) {
      var th = el("th", null, c.label);
      th.scope = "col";
      tr.appendChild(th);
    });
    thead.appendChild(tr);
  }

  function renderTable(reset) {
    if (reset) { tblShown = TBL_PAGE; tblBody.textContent = ""; tblScroll.scrollTop = 0; }
    var start = tblBody.childNodes.length;
    var end = Math.min(tblShown, view.length);
    var frag = document.createDocumentFragment();
    for (var i = start; i < end; i++) {
      var ri = view[i];
      var tr = document.createElement("tr");
      tr.dataset.ri = ri;
      if (ri === state.selected) tr.style.background = "var(--wash)";
      TBL_COLS.forEach(function (c) {
        var td = el("td", c.cls || null, c.get(ri) || "");
        tr.appendChild(td);
      });
      frag.appendChild(tr);
    }
    tblBody.appendChild(frag);
    document.getElementById("tbl-empty").classList.toggle("hidden", view.length > 0);
  }

  function markTableSelection() {
    Array.prototype.forEach.call(tblBody.childNodes, function (tr) {
      tr.style.background = Number(tr.dataset.ri) === state.selected ? "var(--wash)" : "";
    });
  }

  tblScroll.addEventListener("scroll", function () {
    if (state.view !== "table") return;
    if (tblScroll.scrollTop + tblScroll.clientHeight > tblScroll.scrollHeight - 200 && tblShown < view.length) {
      tblShown = Math.min(tblShown + TBL_PAGE, view.length);
      renderTable(false);
    }
  });
  tblBody.addEventListener("click", function (e) {
    var tr = e.target.closest("tr");
    if (tr && tr.dataset.ri !== undefined) select(Number(tr.dataset.ri));
  });

  // ---- downloads ----------------------------------------------------------
  function download(name, text, type) {
    var blob = new Blob([text], { type: type });
    var url = URL.createObjectURL(blob);
    var a = document.createElement("a");
    a.href = url; a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 0);
  }
  document.getElementById("dl-tsv").addEventListener("click", function () {
    var lines = [TBL_COLS.map(function (c) { return c.label; }).join("\t")];
    view.forEach(function (ri) {
      lines.push(TBL_COLS.map(function (c) { return String(c.get(ri) || "").replace(/[\t\n\r]/g, " "); }).join("\t"));
    });
    download(DATA.project + ".filtered.tsv", lines.join("\n") + "\n", "text/tab-separated-values");
  });
  document.getElementById("dl-fa").addEventListener("click", function () {
    var out = [], missing = 0;
    view.forEach(function (ri) {
      var s = ROWS[ri][F.seq];
      if (!s) { missing++; return; }
      out.push(">" + ROWS[ri][F.id] + "\n" + (s.match(/.{1,60}/g) || []).join("\n"));
    });
    if (!out.length) {
      alert("No sequences are embedded for the current selection. Re-run make_report.py " +
            "with --sequences novelties (or all) and --candidates_fa pointing at candidates.fa.");
      return;
    }
    if (missing) {
      // Not fatal: by default only novelty rows carry sequences.
      console.warn(missing + " filtered proteins have no embedded sequence and were skipped.");
    }
    download(DATA.project + ".filtered.faa", out.join("\n") + "\n", "text/plain");
  });

  // ---- summary ------------------------------------------------------------
  function renderSummary() {
    var novTotal = 0, annot = 0;
    var perSpecies = {};
    PROTEOMES.forEach(function (p) { if (p.group === "IN") perSpecies[p.short] = 0; });
    for (var i = 0; i < nRows; i++) {
      if (ROWS[i][F.nov]) {
        novTotal++;
        var s = ROWS[i][F.src] >= 0 ? PROTEOMES[ROWS[i][F.src]].short : null;
        if (s && perSpecies[s] !== undefined) perSpecies[s]++;
      }
      if (ROWS[i][F.fsrc] >= 0) annot++;
    }

    document.getElementById("hero").textContent = novTotal.toLocaleString();
    document.getElementById("hero-note").textContent =
      "Present in ≥ " + Math.round(DATA.ingroup_min_frac * 100) +
      "% of the ingroup and absent from every outgroup proteome.";
    document.getElementById("t-total").textContent = nRows.toLocaleString();
    document.getElementById("t-in").textContent = N_IN + (N_IN === 1 ? " species" : " species");
    document.getElementById("t-out").textContent = (PROTEOMES.length - N_IN) + " species";
    document.getElementById("t-annot").textContent =
      nRows ? Math.round((annot / nRows) * 100) + "%" : "—";

    var bars = document.getElementById("species-bars");
    bars.textContent = "";
    var maxV = Math.max.apply(null, Object.keys(perSpecies).map(function (k) { return perSpecies[k]; }).concat([1]));
    PROTEOMES.filter(function (p) { return p.group === "IN"; }).forEach(function (p) {
      var v = perSpecies[p.short] || 0;
      var name = el("div", "bar-name", p.species + (p.strain ? " " + p.strain : ""));
      name.title = p.short;
      bars.appendChild(name);
      var track = el("div", "bar-track");
      var fill = el("div", "bar-fill");
      fill.style.width = Math.max(2, Math.round((v / maxV) * 220)) + "px";
      track.appendChild(fill);
      track.appendChild(el("span", "bar-value", v.toLocaleString())); // value at the bar tip
      bars.appendChild(track);
    });
  }

  // ---- wiring -------------------------------------------------------------
  function refresh(resetScroll) {
    applyFilters();
    if (state.selected >= 0 && view.indexOf(state.selected) === -1) {
      state.selected = -1;
      renderDetail();
    }
    document.getElementById("count").textContent =
      "Showing " + view.length.toLocaleString() + " of " + nRows.toLocaleString() + " proteins";
    spacerEl.style.height = Math.max(1, view.length * ROW_H) + "px";
    if (resetScroll) scrollEl.scrollTop = 0;
    document.getElementById("grid-empty").classList.toggle("hidden", view.length > 0);
    drawGrid();
    renderTable(true);
  }

  function populateSelects() {
    var src = document.getElementById("f-src");
    src.appendChild(new Option("All source proteomes", ""));
    PROTEOMES.forEach(function (p) {
      src.appendChild(new Option(p.short + " — " + p.species + (p.group === "IN" ? " (in)" : " (out)"), p.short));
    });
    var fs = document.getElementById("f-fsrc");
    fs.appendChild(new Option("Any annotation source", ""));
    DATA.fsources.slice().sort().forEach(function (f) { fs.appendChild(new Option(f, f)); });

    var fam = document.getElementById("f-family");
    fam.appendChild(new Option(
      FAMILIES.length ? "All families (" + FAMILIES.length + ")" : "No multi-species families", ""
    ));
    FAMILIES.forEach(function (f, i) { fam.appendChild(new Option(familyLabel(f), String(i))); });
  }

  function setFamilyFilter(idx) {
    state.family = idx;
    document.getElementById("f-family").value = idx >= 0 ? String(idx) : "";
    refresh(true);
  }

  var searchTimer = null;
  document.getElementById("f-search").addEventListener("input", function (e) {
    var v = e.target.value;
    clearTimeout(searchTimer);
    searchTimer = setTimeout(function () { state.search = v; refresh(true); }, 140);
  });
  document.getElementById("f-src").addEventListener("change", function (e) { state.src = e.target.value; refresh(true); });
  document.getElementById("f-fsrc").addEventListener("change", function (e) { state.fsrc = e.target.value; refresh(true); });
  document.getElementById("f-family").addEventListener("change", function (e) {
    state.family = e.target.value === "" ? -1 : Number(e.target.value);
    refresh(true);
  });
  document.getElementById("f-nov").addEventListener("change", function (e) { state.novOnly = e.target.checked; refresh(true); });
  document.getElementById("f-pfam").addEventListener("change", function (e) { state.pfamOnly = e.target.checked; refresh(true); });
  document.getElementById("f-notb").addEventListener("change", function (e) { state.noTb = e.target.checked; refresh(true); });
  document.getElementById("f-sort").addEventListener("change", function (e) { state.sort = e.target.value; refresh(true); });
  document.getElementById("f-reset").addEventListener("click", function () {
    state.search = ""; state.src = ""; state.fsrc = ""; state.family = -1;
    state.novOnly = true; state.pfamOnly = false; state.noTb = false; state.sort = "ingroup";
    document.getElementById("f-search").value = "";
    document.getElementById("f-src").value = "";
    document.getElementById("f-fsrc").value = "";
    document.getElementById("f-family").value = "";
    document.getElementById("f-nov").checked = true;
    document.getElementById("f-pfam").checked = false;
    document.getElementById("f-notb").checked = false;
    document.getElementById("f-sort").value = "ingroup";
    refresh(true);
  });

  function setView(v) {
    state.view = v;
    document.getElementById("tab-heatmap").setAttribute("aria-selected", String(v === "heatmap"));
    document.getElementById("tab-table").setAttribute("aria-selected", String(v === "table"));
    document.getElementById("panel-heatmap").classList.toggle("hidden", v !== "heatmap");
    document.getElementById("panel-table").classList.toggle("hidden", v !== "table");
    if (v === "heatmap") drawGrid();
  }
  document.getElementById("tab-heatmap").addEventListener("click", function () { setView("heatmap"); });
  document.getElementById("tab-table").addEventListener("click", function () { setView("table"); });

  scrollEl.addEventListener("scroll", function () { drawGrid(); hideTip(); });
  scrollEl.addEventListener("pointermove", function (e) {
    var vi = rowAt(e.clientY);
    var ci = colAt(e.clientX);
    var ri = vi >= 0 ? view[vi] : -1;
    if (ri !== state.hoverRow || ci !== state.hoverCol) {
      state.hoverRow = ri; state.hoverCol = ci;
      drawGrid();
    }
    if (ri >= 0) showTip(ri, ci, e.clientX, e.clientY); else hideTip();
  });
  scrollEl.addEventListener("pointerleave", function () {
    state.hoverRow = -1; state.hoverCol = -1;
    hideTip(); drawGrid();
  });
  scrollEl.addEventListener("click", function (e) {
    var vi = rowAt(e.clientY);
    if (vi >= 0) select(view[vi]);
  });
  // Keyboard parity: arrows move the selection and surface the same detail as hover.
  scrollEl.addEventListener("keydown", function (e) {
    if (!view.length) return;
    var cur = view.indexOf(state.selected);
    var next = cur;
    if (e.key === "ArrowDown") next = Math.min(view.length - 1, cur < 0 ? 0 : cur + 1);
    else if (e.key === "ArrowUp") next = Math.max(0, cur < 0 ? 0 : cur - 1);
    else if (e.key === "Home") next = 0;
    else if (e.key === "End") next = view.length - 1;
    else if (e.key === "PageDown") next = Math.min(view.length - 1, (cur < 0 ? 0 : cur) + 20);
    else if (e.key === "PageUp") next = Math.max(0, (cur < 0 ? 0 : cur) - 20);
    else return;
    e.preventDefault();
    select(view[next]);
    var y = next * ROW_H;
    if (y < scrollEl.scrollTop) scrollEl.scrollTop = y;
    else if (y + ROW_H > scrollEl.scrollTop + scrollEl.clientHeight) {
      scrollEl.scrollTop = y - scrollEl.clientHeight + ROW_H;
    }
  });

  var themeBtn = document.getElementById("theme-toggle");
  themeBtn.addEventListener("click", function () {
    var cur = document.documentElement.getAttribute("data-theme");
    if (!cur) {
      cur = window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
    }
    document.documentElement.setAttribute("data-theme", cur === "dark" ? "light" : "dark");
    drawHead(); drawGrid();
  });
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", function () {
    drawHead(); drawGrid();
  });

  var resizeTimer = null;
  window.addEventListener("resize", function () {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(function () { drawHead(); drawGrid(); }, 100);
  });

  // ---- init ---------------------------------------------------------------
  document.getElementById("title").textContent = DATA.project + " — novelty candidates";
  var inNames = PROTEOMES.filter(function (p) { return p.group === "IN"; })
    .map(function (p) { return p.species; }).join(", ");
  var outNames = PROTEOMES.filter(function (p) { return p.group === "OUT"; })
    .map(function (p) { return p.species; }).join(", ");
  document.getElementById("subtitle").textContent =
    "Ingroup: " + inNames + "  ·  Outgroup: " + outNames;
  document.title = DATA.project + " — NovInvenio candidates";

  populateSelects();
  renderSummary();
  renderTableHead();
  drawHead();
  refresh(true);
  renderDetail();
})();
</script>
</body>
</html>
"""
