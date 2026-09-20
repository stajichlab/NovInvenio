"""
Self-contained HTML/CSS/JS template for the NovInvenio island synteny view.

bin/pangenome_island_synteny.py substitutes two tokens:
  __PROJECT_TITLE__   plain-text project name (HTML-escaped by the caller)
  /*__PAYLOAD__*/     the JSON payload from lib/island_synteny.build_payload()

Same file:// / no-network / textContent-only constraints as
lib/report_template.py (see that module's docstring and CLAUDE.md's
"Constraints to preserve when editing the page"). Built on
lib/core_report_template.py's skeleton -- the smallest of the existing
single-table pages -- but this page draws a canvas grid rather than a table,
so it also borrows the palette()/getComputedStyle repaint pattern from
lib/report_template.py's canvas heatmap.

The view: for ONE selected accessory island, rows are strains collapsed into
distinct haplotypes (identical presence/absence rows merged, carrying a
"x N" count badge) and columns are the island's member gene families IN
LOCUS ORDER. A filled cell means that family is present in that haplotype.
Because haplotype rows are sorted on their 0/1 pattern by default, a
deletion breakpoint shows up as a shared vertical edge across rows rather
than scattered gaps.

Design note -- why the glyph strip is one colour per island, not per family:
lib.island_synteny.build_payload() records `dominant_class` and `domains` at
the ISLAND level (the union of Pfam domain names across every member
family), not per family -- there is no column-by-column class breakdown in
the payload to draw from. So every glyph-strip tick for a given island is
coloured by that island's single `dominant_class`, and hovering a tick shows
the family ID for that column plus the island's full domain string (matching
what the payload actually carries, not an invented per-family split).

Design note -- glyph/class colour tokens: lib/skins.py's REQUIRED_TOKENS has
exactly two data-series colours (--series-1, --series-2), both already
spoken for elsewhere (protein-search presence / TBLASTN hit). Rather than
add a five-way categorical palette to every skin for a single new page, the
five fixed lib.pfam_classes.CLASS_ORDER keys are mapped to existing tokens
(see CLASS_TOKENS below) -- still entirely skin-driven, no new hex anywhere.
"""

from report_common import (
    BASE_PAGE_CSS,
    BREADCRUMB_NAV_CSS,
    EL_HELPER_JS,
    FAVICON_LINK_HTML,
    FOOTER_CSS,
    FOOTER_HTML,
    LOGO_CSS,
    LOGO_IMG_HTML,
    SKIN_BOOT_JS,
    SKIN_PICKER_HTML,
    SKIN_PICKER_JS,
    SKIN_VARS_CSS,
    breadcrumb_nav_html,
)

ISLAND_SYNTENY_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__PROJECT_TITLE__ — NovInvenio island synteny</title>
""" + FAVICON_LINK_HTML + r"""
<style>
""" + SKIN_VARS_CSS + BASE_PAGE_CSS + LOGO_CSS + BREADCRUMB_NAV_CSS + FOOTER_CSS + r"""
  /* This page selects skins via the `data-skin` attribute (lib/skins.py),
     never `data-theme` -- noted here, not left for a reader to guess, since
     `data-theme` is the attribute name a couple of other in-house HTML
     conventions use for the same idea. */
  :root { --bg: var(--page); } /* alias kept for tooling that greps for --bg */

  .isv-explorer { display: grid; grid-template-columns: 300px 1fr; gap: 16px; align-items: start; }
  @media (max-width: 1080px) { .isv-explorer { grid-template-columns: 1fr; } }

  .isv-sidebar { position: sticky; top: 16px; max-height: calc(100vh - 32px); display: flex; flex-direction: column; }
  .isv-sidebar-controls { display: flex; flex-direction: column; gap: 8px; margin-bottom: 10px; }
  .isv-sidebar-controls input, .isv-sidebar-controls select { width: 100%; }
  .isv-list { overflow-y: auto; flex: 1; border: 1px solid var(--border); border-radius: 8px; }
  .isv-item {
    display: block; width: 100%; text-align: left; border: none; border-radius: 0;
    border-bottom: 1px solid var(--border); background: var(--surface-1); padding: 10px 12px;
  }
  .isv-item:last-child { border-bottom: none; }
  .isv-item:hover { background: var(--hover-wash); }
  .isv-item.sel { background: var(--wash); box-shadow: inset 3px 0 0 var(--series-1); }
  .isv-item-id { font-family: var(--font-mono); font-size: 12px; font-weight: 600; word-break: break-all; }
  .isv-item-stats { color: var(--text-secondary); font-size: 11.5px; margin-top: 3px; }
  .isv-empty-list { padding: 16px 12px; color: var(--text-secondary); font-size: 12px; }

  .isv-chip {
    display: inline-block; padding: 1px 7px; border-radius: 999px; font-size: 10.5px;
    border: 1px solid var(--border); margin-top: 5px;
  }

  .isv-main-head { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; margin-bottom: 12px; }
  .isv-main-head h2 { margin: 0; font-size: 15px; font-weight: 600; font-family: var(--font-mono); }
  .isv-main-note { color: var(--text-secondary); font-size: 12px; margin: 0 0 12px; }

  .isv-legend { display: flex; flex-wrap: wrap; gap: 14px; margin-bottom: 10px; }
  .isv-legend-item { display: inline-flex; align-items: center; gap: 6px; font-size: 11.5px; color: var(--text-secondary); }
  .isv-swatch { width: 11px; height: 11px; border-radius: 2px; flex: none; }

  .isv-hscroll { overflow-x: auto; border: 1px solid var(--border); border-radius: 8px; }
  canvas.isv-head { display: block; }
  .isv-vscroll { overflow-y: auto; max-height: 520px; border-top: 1px solid var(--border); }
  canvas.isv-grid { display: block; }

  .isv-placeholder, .isv-empty-state {
    padding: 40px 20px; color: var(--text-secondary); font-size: 13px; text-align: center;
    border: 1px dashed var(--border); border-radius: 10px;
  }
  .isv-empty-state h3 { margin: 0 0 8px; color: var(--text-primary); font-size: 15px; }

  @media print {
    .filters, header.top select, header.top button, .isv-sidebar-controls { display: none !important; }
  }
</style>
<script>""" + SKIN_BOOT_JS + r"""</script>
</head>
<body>
<div class="wrap">
  """ + breadcrumb_nav_html() + r"""
  <header class="top">
    """ + LOGO_IMG_HTML + r"""
    <div class="titles">
      <h1 id="title"></h1>
      <p class="sub" id="subtitle"></p>
    </div>
""" + SKIN_PICKER_HTML + r"""
  </header>

  <section class="card">
    <h2 class="card-title">Accessory island synteny</h2>
    <p class="card-note" id="summary-note"></p>
    <div class="tiles">
      <div>
        <div class="tile-label">Islands shown</div>
        <div class="tile-value" id="t-shown"></div>
      </div>
      <div>
        <div class="tile-label">Located islands (total)</div>
        <div class="tile-value" id="t-total"></div>
      </div>
      <div>
        <div class="tile-label">Excluded</div>
        <div class="tile-value" id="t-excluded"></div>
      </div>
      <div>
        <div class="tile-label">Strains</div>
        <div class="tile-value" id="t-strains"></div>
      </div>
    </div>
  </section>

  <div id="isv-body">
    <div class="isv-explorer" id="isv-explorer">
      <aside class="isv-sidebar">
        <div class="isv-sidebar-controls">
          <input type="search" id="f-search" placeholder="Search locus or family ID…" aria-label="Search islands">
          <select id="f-sidebar-sort" aria-label="Sort islands by">
            <option value="size">Sort: size (families)</option>
            <option value="haplotypes">Sort: haplotype count</option>
            <option value="strains">Sort: carrying strains</option>
            <option value="name">Sort: locus ID</option>
          </select>
          <span class="count" id="sidebar-count" role="status" aria-live="polite"></span>
        </div>
        <div class="isv-list" id="isv-list" role="listbox" aria-label="Selected accessory islands"></div>
      </aside>

      <section class="card" id="isv-main" style="padding:14px">
        <div class="isv-main-head">
          <h2 id="isv-title"></h2>
          <div class="spacer"></div>
          <select id="f-row-sort" aria-label="Sort haplotype rows">
            <option value="pattern">Sort rows: presence pattern</option>
            <option value="count">Sort rows: strain count</option>
            <option value="strain">Sort rows: strain name</option>
          </select>
        </div>
        <p class="isv-main-note" id="isv-main-note"></p>
        <div class="isv-legend" id="isv-legend"></div>

        <div id="isv-canvas-area">
          <div class="isv-hscroll" id="hscroll">
            <canvas class="isv-head" id="glyphs"></canvas>
            <div class="isv-vscroll" id="vscroll">
              <canvas class="isv-grid" id="grid"></canvas>
            </div>
          </div>
        </div>
        <p class="isv-placeholder hidden" id="isv-no-selection">Select an island from the list to view its synteny grid.</p>
      </section>
    </div>

    <div class="isv-empty-state hidden" id="isv-empty-state">
      <h3>No accessory islands to show</h3>
      <p id="isv-empty-text"></p>
    </div>
  </div>
""" + FOOTER_HTML + r"""
</div>

<div class="tip" id="tip" role="tooltip" style="position:fixed;z-index:20;max-width:340px;
  background:var(--surface-1);border:1px solid var(--border);border-radius:8px;
  box-shadow:var(--shadow);padding:10px 12px;font-size:12px;pointer-events:none;display:none;"></div>

<script type="application/json" id="payload">/*__PAYLOAD__*/</script>
<script>
(function () {
  "use strict";

  var DATA = JSON.parse(document.getElementById("payload").textContent);
  var ISLANDS = DATA.islands || [];
  var CLASSES = DATA.classes || {};

""" + EL_HELPER_JS + r"""

  // Fixed vocabulary (lib/pfam_classes.py's CLASS_ORDER) -> an existing skin
  // token. See this module's docstring for why the mapping is fixed rather
  // than data-driven: the page has exactly two data-series colours to work
  // with, both already meaning something else, so five classes reuse tokens
  // whose *role* is "a distinguishable accent", not literally "series 1/2".
  var CLASS_TOKENS = {
    nlr: "--series-1",
    secondary_metabolite: "--series-2",
    transporter: "--warn",
    other: "--text-secondary",
    unannotated: "--muted"
  };

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
      wash: css("--wash"),
      hover: css("--hover-wash")
    };
  }
  function classColor(key) {
    return css(CLASS_TOKENS[key] || "--muted");
  }
  function classLabel(key) {
    return CLASSES[key] || key || "Unclassified";
  }

  // ---- sidebar state --------------------------------------------------------
  var state = {
    search: "",
    sidebarSort: "size",
    rowSort: "pattern",
    selected: ISLANDS.length ? 0 : -1
  };
  var sidebarView = [];

  function islandHaystack(isl) {
    return (isl.locus_id + " " + isl.families.join(" ")).toLowerCase();
  }

  function applySidebarFilter() {
    var q = state.search.trim().toLowerCase();
    sidebarView = [];
    for (var i = 0; i < ISLANDS.length; i++) {
      if (q && islandHaystack(ISLANDS[i]).indexOf(q) === -1) continue;
      sidebarView.push(i);
    }
    var s = state.sidebarSort;
    var cmp;
    if (s === "haplotypes") cmp = function (a, b) { return ISLANDS[b].haplotypes.length - ISLANDS[a].haplotypes.length; };
    else if (s === "strains") cmp = function (a, b) { return ISLANDS[b].n_strains - ISLANDS[a].n_strains; };
    else if (s === "name") cmp = function (a, b) { return ISLANDS[a].locus_id < ISLANDS[b].locus_id ? -1 : 1; };
    else cmp = function (a, b) { return ISLANDS[b].size - ISLANDS[a].size; };
    sidebarView.sort(cmp);
  }

  function renderSidebar() {
    var list = document.getElementById("isv-list");
    list.textContent = "";
    if (!sidebarView.length) {
      list.appendChild(el("p", "isv-empty-list",
        ISLANDS.length ? "No islands match the current search." : "No islands selected for this run."));
    } else {
      sidebarView.forEach(function (idx) {
        var isl = ISLANDS[idx];
        var btn = el("button", "isv-item");
        btn.type = "button";
        btn.setAttribute("role", "option");
        btn.setAttribute("aria-selected", idx === state.selected ? "true" : "false");
        if (idx === state.selected) btn.classList.add("sel");
        btn.appendChild(el("div", "isv-item-id", isl.locus_id));
        btn.appendChild(el("div", "isv-item-stats",
          isl.size + " families · " + isl.haplotypes.length + " haplotypes · " +
          isl.n_strains + " strains"));
        var chip = el("span", "isv-chip", classLabel(isl.dominant_class));
        chip.style.borderColor = classColor(isl.dominant_class);
        chip.style.color = classColor(isl.dominant_class);
        btn.appendChild(chip);
        btn.addEventListener("click", function () { selectIsland(idx); });
        list.appendChild(btn);
      });
    }
    document.getElementById("sidebar-count").textContent =
      sidebarView.length.toLocaleString() + " of " + ISLANDS.length.toLocaleString() + " shown";
  }

  function selectIsland(idx) {
    state.selected = idx;
    renderSidebar();
    renderMain();
  }

  // ---- row sorting ------------------------------------------------------
  function sortedHaplotypes(isl) {
    var haps = isl.haplotypes.slice();
    if (state.rowSort === "count") {
      haps.sort(function (a, b) { return b.count - a.count || (a.pattern < b.pattern ? -1 : 1); });
    } else if (state.rowSort === "strain") {
      haps.sort(function (a, b) {
        var an = a.strains[0] || "", bn = b.strains[0] || "";
        return an < bn ? -1 : an > bn ? 1 : 0;
      });
    } // "pattern" -- the payload already sorts haplotypes lexicographically
      // on the 0/1 pattern (lib/island_synteny.py's collapse_haplotypes),
      // which is what lines deletion breakpoints up into a shared edge.
    return haps;
  }

  // ---- canvas geometry ----------------------------------------------------
  var ROW_H = 20;
  var CELL_W = 22;
  var GUTTER = 170;
  var GLYPH_H = 14;
  var LABEL_ANGLE = Math.PI / 3;
  var HEAD_H = 74;

  function colX(i) { return GUTTER + i * CELL_W; }

  var glyphCanvas = document.getElementById("glyphs");
  var gridCanvas = document.getElementById("grid");
  var hctx = glyphCanvas.getContext("2d");
  var gctx = gridCanvas.getContext("2d");

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

  var LABEL_FONT = "600 10px system-ui, -apple-system, 'Segoe UI', sans-serif";

  function drawGlyphStrip(isl, totalW) {
    var P = palette();
    sizeCanvas(glyphCanvas, hctx, totalW, HEAD_H);
    hctx.clearRect(0, 0, totalW, HEAD_H);
    hctx.fillStyle = P.surface;
    hctx.fillRect(0, 0, totalW, HEAD_H);

    var tickColor = classColor(isl.dominant_class);
    isl.families.forEach(function (fam, i) {
      var x = colX(i);
      hctx.fillStyle = tickColor;
      hctx.fillRect(x + 3, 6, CELL_W - 6, GLYPH_H);

      hctx.save();
      hctx.translate(x + CELL_W / 2, HEAD_H - 8);
      hctx.rotate(-LABEL_ANGLE);
      hctx.fillStyle = P.primary;
      hctx.font = LABEL_FONT;
      hctx.textAlign = "left";
      hctx.textBaseline = "middle";
      hctx.fillText(ellipsize(hctx, fam, 110), 0, 0);
      hctx.restore();
    });

    hctx.fillStyle = P.axis;
    hctx.fillRect(0, GLYPH_H + 10, totalW, 1);
  }

  function drawGrid(isl, haps, totalW) {
    var P = palette();
    var h = Math.max(ROW_H, haps.length * ROW_H);
    sizeCanvas(gridCanvas, gctx, totalW, h);
    gctx.clearRect(0, 0, totalW, h);
    gctx.fillStyle = P.surface;
    gctx.fillRect(0, 0, totalW, h);

    var presentColor = classColor(isl.dominant_class);

    haps.forEach(function (hap, ri) {
      var y = ri * ROW_H;

      gctx.textBaseline = "middle";
      gctx.textAlign = "left";
      gctx.font = "600 11px ui-monospace, SFMono-Regular, Menlo, monospace";
      gctx.fillStyle = P.primary;
      gctx.fillText("×" + hap.count, 8, y + ROW_H / 2);

      if (hap.strains.length) {
        gctx.font = "10px system-ui, -apple-system, 'Segoe UI', sans-serif";
        gctx.fillStyle = P.secondary;
        var note = hap.count === 1 ? hap.strains[0] : hap.strains[0] + " +" + (hap.count - 1);
        gctx.fillText(ellipsize(gctx, note, GUTTER - 60), 50, y + ROW_H / 2);
      }

      for (var ci = 0; ci < hap.pattern.length; ci++) {
        var on = hap.pattern.charAt(ci) === "1";
        gctx.fillStyle = on ? presentColor : P.grid;
        gctx.fillRect(colX(ci) + 1, y + 1, CELL_W - 2, ROW_H - 2);
      }
    });

    // Deletion breakpoints: a vertical edge wherever consecutive columns do
    // not carry the same presence pattern across every drawn haplotype row.
    if (haps.length > 1) {
      for (var b = 1; b < isl.families.length; b++) {
        var differs = false;
        for (var r2 = 0; r2 < haps.length; r2++) {
          if (haps[r2].pattern.charAt(b - 1) !== haps[r2].pattern.charAt(b)) { differs = true; break; }
        }
        if (differs) {
          gctx.fillStyle = P.axis;
          gctx.fillRect(colX(b), 0, 1, h);
        }
      }
    }
  }

  function totalWidth(isl) {
    return colX(isl.families.length) + 8;
  }

  // ---- glyph-strip hover tooltip -------------------------------------------
  var tipEl = document.getElementById("tip");
  function colAtGlyph(clientX, isl) {
    var rect = glyphCanvas.getBoundingClientRect();
    var x = clientX - rect.left;
    for (var i = 0; i < isl.families.length; i++) {
      var cx = colX(i);
      if (x >= cx && x < cx + CELL_W) return i;
    }
    return -1;
  }
  glyphCanvas.addEventListener("mousemove", function (e) {
    var isl = ISLANDS[state.selected];
    if (!isl) return;
    var ci = colAtGlyph(e.clientX, isl);
    if (ci < 0) { tipEl.style.display = "none"; return; }
    tipEl.textContent = "";
    tipEl.appendChild(el("div", "tip-id", isl.families[ci]));
    tipEl.appendChild(el("div", null, classLabel(isl.dominant_class) + " (island dominant class)"));
    if (isl.domains.length) {
      tipEl.appendChild(el("div", null, "Domains in this island: " + isl.domains.join(", ")));
    }
    tipEl.style.display = "block";
    var w = tipEl.offsetWidth, hgt = tipEl.offsetHeight;
    var left = e.clientX + 14, top = e.clientY + 14;
    if (left + w > window.innerWidth - 8) left = e.clientX - w - 14;
    if (top + hgt > window.innerHeight - 8) top = e.clientY - hgt - 14;
    tipEl.style.left = Math.max(8, left) + "px";
    tipEl.style.top = Math.max(8, top) + "px";
  });
  glyphCanvas.addEventListener("mouseleave", function () { tipEl.style.display = "none"; });

  // ---- legend ---------------------------------------------------------------
  function renderLegend(isl) {
    var legend = document.getElementById("isv-legend");
    legend.textContent = "";
    var item = el("span", "isv-legend-item");
    var sw = el("span", "isv-swatch");
    sw.style.background = classColor(isl.dominant_class);
    item.appendChild(sw);
    item.appendChild(el("span", null, "Present (" + classLabel(isl.dominant_class) + ")"));
    legend.appendChild(item);

    var item2 = el("span", "isv-legend-item");
    var sw2 = el("span", "isv-swatch");
    sw2.style.background = css("--grid");
    item2.appendChild(sw2);
    item2.appendChild(el("span", null, "Absent"));
    legend.appendChild(item2);
  }

  // ---- top-level render -----------------------------------------------------
  function renderMain() {
    var isl = ISLANDS[state.selected];
    var canvasArea = document.getElementById("isv-canvas-area");
    var noSel = document.getElementById("isv-no-selection");
    if (!isl) {
      canvasArea.classList.add("hidden");
      noSel.classList.remove("hidden");
      document.getElementById("isv-title").textContent = "";
      document.getElementById("isv-main-note").textContent = "";
      document.getElementById("isv-legend").textContent = "";
      return;
    }
    canvasArea.classList.remove("hidden");
    noSel.classList.add("hidden");

    document.getElementById("isv-title").textContent = isl.locus_id;
    document.getElementById("isv-main-note").textContent =
      isl.size + " families in locus order · " + isl.n_strains + " strains carry this island · " +
      isl.haplotypes.length + " distinct presence patterns" +
      (isl.locus_contig ? " · " + isl.locus_contig +
        (isl.locus_start >= 0 ? ":" + isl.locus_start + "-" + isl.locus_end : "") : "");

    renderLegend(isl);

    var haps = sortedHaplotypes(isl);
    var w = totalWidth(isl);
    drawGlyphStrip(isl, w);
    drawGrid(isl, haps, w);
  }

  window.onSkinChange = function () {
    if (state.selected >= 0) renderMain();
  };

  // ---- empty state (no islands at all in this payload) -----------------------
  function renderEmptyState() {
    var explorer = document.getElementById("isv-explorer");
    var emptyState = document.getElementById("isv-empty-state");
    if (!ISLANDS.length) {
      explorer.classList.add("hidden");
      emptyState.classList.remove("hidden");
      document.getElementById("isv-empty-text").textContent =
        DATA.n_islands_total
          ? DATA.n_islands_total + " accessory island" + (DATA.n_islands_total === 1 ? "" : "s") +
            " were located, but " + DATA.n_islands_excluded + " " +
            (DATA.n_islands_excluded === 1 ? "was" : "were") +
            " excluded (carried by too few strains to compare) and none remain to draw."
          : "No accessory islands were located for this run.";
    } else {
      explorer.classList.remove("hidden");
      emptyState.classList.add("hidden");
    }
  }

  // ---- wiring -----------------------------------------------------------
  document.getElementById("f-search").addEventListener("input", function (e) {
    state.search = e.target.value;
    applySidebarFilter();
    renderSidebar();
  });
  document.getElementById("f-sidebar-sort").addEventListener("change", function (e) {
    state.sidebarSort = e.target.value;
    applySidebarFilter();
    renderSidebar();
  });
  document.getElementById("f-row-sort").addEventListener("change", function (e) {
    state.rowSort = e.target.value;
    renderMain();
  });

""" + SKIN_PICKER_JS + r"""

  // ---- init ---------------------------------------------------------------
  document.getElementById("title").textContent = DATA.project + " — island synteny";
  document.getElementById("subtitle").textContent =
    ISLANDS.length + " of " + DATA.n_islands_total + " located accessory islands shown, in locus order.";
  document.title = DATA.project + " — NovInvenio island synteny";

  document.getElementById("summary-note").textContent =
    DATA.n_islands_excluded
      ? DATA.n_islands_excluded + " located island" + (DATA.n_islands_excluded === 1 ? "" : "s") +
        " excluded from this view for being carried by too few strains to show a meaningful " +
        "presence pattern (a single-strain island is one filled row with no breakpoint in it)."
      : "Every located accessory island met the minimum strain-count filter for this view.";
  document.getElementById("t-shown").textContent = ISLANDS.length.toLocaleString();
  document.getElementById("t-total").textContent = (DATA.n_islands_total || 0).toLocaleString();
  document.getElementById("t-excluded").textContent = (DATA.n_islands_excluded || 0).toLocaleString();
  document.getElementById("t-strains").textContent = (DATA.strains || []).length.toLocaleString();

  renderEmptyState();
  if (ISLANDS.length) {
    applySidebarFilter();
    renderSidebar();
    renderMain();
  }
})();
</script>
</body>
</html>
"""
