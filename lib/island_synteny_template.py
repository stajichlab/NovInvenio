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

The glyph strip is PER FAMILY (per column): each tick is coloured by that
column's own `family_classes[i]` (lib.island_synteny.build_payload(), keyed
off bin/pangenome_domain_enrichment.py's per-family Pfam scan), not the
island's single `dominant_class` -- a family absent from that scan, or a run
with no domain-enrichment data supplied at all, is "unannotated" rather than
an error. Hovering a tick reveals that column's own family ID and domain
string (`families[i]` / `family_domains[i]`).

Domain annotation stays OUT of the grid cells -- this is a spec requirement
(docs/superpowers/specs/2026-09-19-pangenome-gainloss-visualization-design.md),
not a style choice: "one colored tick per column keyed to a Pfam class ...
not in cells". The main grid's "present" cell fill is always a single,
uniform colour (`P.primary`); it never varies by family class, both because
colour must encode evidence type rather than let one hue carry presence AND
annotation class at once (CLAUDE.md's report colour rule), and because the
view's whole purpose -- spotting a deletion breakpoint, the vertical edge
where a run of filled cells stops across every haplotype row -- is easiest
to read against a uniform field. The island-level `dominant_class` is
carried in the payload and used only for the sidebar's one-chip-per-island
summary, where a single colour per island is what that list is for.

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
<!--__DIAGNOSTICS_BANNER__-->
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
        <div class="tile-label">Truncated</div>
        <div class="tile-value" id="t-truncated"></div>
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
        <!-- Issue #119: a "species" <option> is appended to #f-row-sort by
             JS, only when DATA.species has entries (see the init block
             below) -- mirrors novelties.html's f-category filter, which
             stays absent unless the payload actually carries category
             data. A pairwise/mmseqs run with no --config never gets an
             affordance for a sort it cannot perform. -->

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
  // {Short: Species}, issue #119 -- optional, only present when
  // bin/pangenome_island_synteny.py was run with --config. Empty ({}) means
  // "no species data for this run", not an error.
  var SPECIES = DATA.species || {};
  function speciesOf(strain) {
    return SPECIES[strain] || "";
  }

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
  // The currently-drawn (sorted) haplotype rows, kept for the grid-row hover
  // tooltip below -- rowAtGrid() needs to map a pixel Y back to the SAME
  // array drawGrid() just rendered, not re-sort on every mousemove.
  var currentHaps = [];

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

  // The three functions below are written to take `speciesMap` as an
  // explicit argument, deliberately NOT closing over the module-level
  // SPECIES var, so tests/test_island_synteny_species_sort.py can extract
  // them verbatim (via a brace-matching slice, not a JS parser) and run
  // them under plain `node`, with no jsdom/DOM dependency -- these
  // functions never touch the DOM, so there is no reason to require it just
  // to test them.

  // {species-or-"" : count} across a haplotype's strains. A strain absent
  // from speciesMap (an incomplete --config) counts toward the "" bucket,
  // same as an explicitly unknown one -- it must never silently become some
  // OTHER species just because it's missing.
  function speciesCounts(hap, speciesMap) {
    var counts = {};
    hap.strains.forEach(function (s) {
      var sp = speciesMap[s] || "";
      counts[sp] = (counts[sp] || 0) + 1;
    });
    return counts;
  }

  // A haplotype row can carry strains of more than one species (NOT a rare
  // edge case for this feature's own target study -- Coccidioides immitis
  // vs posadasii commonly share a presence pattern, so a mixed row is the
  // common, scientifically interesting case, not the exception). Its
  // "species" for BANDING purposes is the MODE among its KNOWN-species
  // strains, tied broken by the lexicographically smallest species name --
  // deterministic, and cheap since a haplotype's strain list is already
  // small. Unknown ("" -- a strain missing from speciesMap, i.e. absent
  // from the config CSV) is an ABSENCE of data and must never outrank an
  // actually observed species: "" is excluded from the modal contest
  // entirely and is only ever returned when EVERY strain in the haplotype
  // is unmapped. This function only decides which band a row sorts into;
  // it does NOT claim the row is pure -- see speciesBandLabel() below for
  // the part that surfaces a mixed row (known or unknown) rather than
  // hiding it. Unknown-species rows still band last at the sort site (the
  // "￿" sentinel in sortedHaplotypes below), so an incomplete --config
  // doesn't shove unlabelled rows to the top.
  function haplotypeSpecies(hap, speciesMap) {
    var counts = speciesCounts(hap, speciesMap);
    var best = "", bestCount = -1;
    Object.keys(counts).sort().forEach(function (sp) {
      if (sp === "") return; // "no data" never wins a tie against real data
      if (counts[sp] > bestCount) { bestCount = counts[sp]; best = sp; }
    });
    if (best === "" && counts[""]) return ""; // every strain unmapped
    return best;
  }

  // The row's own label for the species column/tooltip -- the one place
  // that must never let the band claim a purity the row does not have. A
  // pure-species (or all-unknown) row is just its species name (or "Unknown
  // species"); a row spanning N>1 distinct species-or-unknown groups gets
  // "<modal species> (+<N-1> other[s])" so a reader sorting by species can
  // see, at a glance, which rows actually cross the immitis/posadasii
  // boundary rather than being told -- by a band that looks pure -- that
  // they don't.
  function speciesBandLabel(hap, speciesMap) {
    var counts = speciesCounts(hap, speciesMap);
    var groups = Object.keys(counts).length;
    var modal = haplotypeSpecies(hap, speciesMap);
    var modalLabel = modal || "Unknown species";
    var others = groups - 1;
    if (others <= 0) return modalLabel;
    return modalLabel + " (+" + others + (others === 1 ? " other" : " others") + ")";
  }

  function sortedHaplotypes(isl) {
    var haps = isl.haplotypes.slice();
    if (state.rowSort === "count") {
      haps.sort(function (a, b) { return b.count - a.count || (a.pattern < b.pattern ? -1 : 1); });
    } else if (state.rowSort === "strain") {
      haps.sort(function (a, b) {
        var an = a.strains[0] || "", bn = b.strains[0] || "";
        return an < bn ? -1 : an > bn ? 1 : 0;
      });
    } else if (state.rowSort === "species") {
      // Issue #119, spec's "by species (immitis/posadasii band)" sort --
      // GROUPING strains of the same species together, not a phylogeny
      // ordering (this pipeline has no strain tree; see
      // todo/pangenome-phylogeny-aware-gain-loss.md). Primary key: species
      // name, unknown-species rows sorted last via the "￿" sentinel.
      // Secondary key (within a species band): strain count descending, so
      // the haplotype carried by the most strains of that species reads
      // first. Tertiary key: first strain name, purely for a stable,
      // reproducible order among same-count haplotypes.
      haps.sort(function (a, b) {
        var sa = haplotypeSpecies(a, SPECIES) || "￿";
        var sb = haplotypeSpecies(b, SPECIES) || "￿";
        if (sa !== sb) return sa < sb ? -1 : 1;
        if (b.count !== a.count) return b.count - a.count;
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
  // Label header height is computed per island from the longest label
  // (labelLayout), not fixed: a fixed 74 px header cut off any label longer
  // than ~70 px, and family IDs ("Short|protein-id") are usually longer.
  var LABEL_MAX_W = 230;

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

  var LABEL_FONT = "600 11px system-ui, -apple-system, 'Segoe UI', sans-serif";

  // Rotated labels rise LABEL_ANGLE from the column centre. Adjacent labels
  // are CELL_W * sin(LABEL_ANGLE) ~ 19 px apart perpendicular to the text,
  // so 11 px text does not overlap. The header is as tall as the longest
  // label needs, and as wide as the last label's overhang to the right.
  function labelLayout(isl) {
    hctx.font = LABEL_FONT;
    var maxW = 0;
    isl.families.forEach(function (fam) {
      maxW = Math.max(maxW, Math.min(LABEL_MAX_W, hctx.measureText(fam).width));
    });
    return {
      head: Math.ceil(maxW * Math.sin(LABEL_ANGLE)) + GLYPH_H + 26,
      extraW: Math.ceil(maxW * Math.cos(LABEL_ANGLE)) + 12
    };
  }

  // A family's own class, guarding against an older/short payload (no
  // per-family arrays) rather than throwing mid-render.
  function familyClassAt(isl, i) {
    return (isl.family_classes && isl.family_classes[i]) || "unannotated";
  }
  function familyDomainsAt(isl, i) {
    return (isl.family_domains && isl.family_domains[i]) || "";
  }

  function drawGlyphStrip(isl, totalW) {
    var P = palette();
    var L = labelLayout(isl);
    var H = L.head, W = totalW + L.extraW;
    var glyphY = H - GLYPH_H - 6;
    sizeCanvas(glyphCanvas, hctx, W, H);
    hctx.clearRect(0, 0, W, H);
    hctx.fillStyle = P.surface;
    hctx.fillRect(0, 0, W, H);

    isl.families.forEach(function (fam, i) {
      var x = colX(i);
      // Class-coloured glyph sits just above the grid; the label rises from
      // above the glyph, so the two never overlap.
      hctx.fillStyle = classColor(familyClassAt(isl, i));
      hctx.fillRect(x + 3, glyphY, CELL_W - 6, GLYPH_H);

      hctx.save();
      hctx.translate(x + CELL_W / 2, glyphY - 6);
      hctx.rotate(-LABEL_ANGLE);
      hctx.fillStyle = P.primary;
      hctx.font = LABEL_FONT;
      hctx.textAlign = "left";
      hctx.textBaseline = "middle";
      hctx.fillText(ellipsize(hctx, fam, LABEL_MAX_W), 0, 0);
      hctx.restore();
    });

    hctx.fillStyle = P.axis;
    hctx.fillRect(0, H - 1, W, 1);
  }

  function drawGrid(isl, haps, totalW) {
    var P = palette();
    var h = Math.max(ROW_H, haps.length * ROW_H);
    sizeCanvas(gridCanvas, gctx, totalW, h);
    gctx.clearRect(0, 0, totalW, h);
    gctx.fillStyle = P.surface;
    gctx.fillRect(0, 0, totalW, h);

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
        // While sorted by species, the row label IS the species band --
        // speciesBandLabel() names every row spanning more than one
        // species-or-unknown group ("<modal> (+N other[s])") instead of
        // letting a mixed row render identically to a pure one (issue #119
        // review). Both branches draw straight to canvas via fillText, the
        // same mechanism this grid already uses for every other per-row/
        // per-column label, including untrusted strings (strain and family
        // IDs) -- canvas text has no HTML-injection surface, so this is
        // exactly as safe as the existing strain-name note it replaces.
        var note = (state.rowSort === "species" && Object.keys(SPECIES).length)
          ? speciesBandLabel(hap, SPECIES)
          : (hap.count === 1 ? hap.strains[0] : hap.strains[0] + " +" + (hap.count - 1));
        gctx.fillText(ellipsize(gctx, note, GUTTER - 60), 50, y + ROW_H / 2);
      }

      // Deliberately ONE colour for every present cell, never the per-column
      // Pfam class -- see this module's docstring. Annotation lives in the
      // glyph strip above, not in the cells: a multicoloured grid would make
      // the deletion-breakpoint edge (a run of filled cells stopping across
      // every haplotype row) harder to see, which is the one thing this
      // grid exists to show.
      for (var ci = 0; ci < hap.pattern.length; ci++) {
        var on = hap.pattern.charAt(ci) === "1";
        gctx.fillStyle = on ? P.primary : P.grid;
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
  // Column position in the island's example strain (payload
  // family_locations, from gene_positions + the tier-1 cluster table). An
  // older payload without it, or a column with no annotated protein in that
  // strain (e.g. a rescued genome_only call), says so instead of guessing.
  function appendLocation(isl, ci) {
    var locs = isl.family_locations;
    if (!locs) return;
    var loc = locs[ci];
    var who = isl.example_strain || "the example strain";
    if (!loc) {
      tipEl.appendChild(el("div", null, "No annotated gene in " + who));
      return;
    }
    if (loc.rescued) {
      tipEl.appendChild(el("div", null, "No annotated gene in " + who +
        "; TBLASTN rescue hit (genome_only) at " + loc.contig + ":" + loc.start.toLocaleString()));
      return;
    }
    tipEl.appendChild(el("div", null, "Location in " + who + ": " + loc.contig + ":" +
      loc.start.toLocaleString() + "-" + loc.end.toLocaleString() +
      " (" + (loc.end - loc.start + 1).toLocaleString() + " bp)"));
    tipEl.appendChild(el("div", null, "Protein: " + loc.protein +
      (loc.n_copies > 1 ? " (" + loc.n_copies + " copies in this strain)" : "")));
  }
  function appendColumnInfo(isl, ci) {
    tipEl.appendChild(el("div", "tip-id", isl.families[ci]));
    tipEl.appendChild(el("div", null, "Column " + (ci + 1) + " of " + isl.families.length));
    appendLocation(isl, ci);
    tipEl.appendChild(el("div", null, classLabel(familyClassAt(isl, ci))));
    var doms = familyDomainsAt(isl, ci);
    tipEl.appendChild(el("div", null, doms ? "Domains: " + doms : "No annotated Pfam domain"));
  }
  glyphCanvas.addEventListener("mousemove", function (e) {
    var isl = ISLANDS[state.selected];
    if (!isl) return;
    var ci = colAtGlyph(e.clientX, isl);
    if (ci < 0) { tipEl.style.display = "none"; return; }
    tipEl.textContent = "";
    appendColumnInfo(isl, ci);
    positionTip(e);
  });
  glyphCanvas.addEventListener("mouseleave", function () { tipEl.style.display = "none"; });

  function positionTip(e) {
    tipEl.style.display = "block";
    var w = tipEl.offsetWidth, hgt = tipEl.offsetHeight;
    var left = e.clientX + 14, top = e.clientY + 14;
    if (left + w > window.innerWidth - 8) left = e.clientX - w - 14;
    if (top + hgt > window.innerHeight - 8) top = e.clientY - hgt - 14;
    tipEl.style.left = Math.max(8, left) + "px";
    tipEl.style.top = Math.max(8, top) + "px";
  }

  // ---- grid-row hover tooltip (issue #119 review) ----------------------
  // The canvas row label (drawGrid, above) already names a mixed row
  // ("<modal> (+N other[s])") so the band never silently claims a purity it
  // doesn't have. This tooltip is the accessible-text counterpart: hovering
  // any row spells out the FULL per-species strain breakdown via
  // el()/textContent (never innerHTML) rather than leaving the reader to
  // infer it from an abbreviated canvas label. Only wired up when the
  // payload actually carries species data -- a pairwise/mmseqs run with no
  // --config has nothing to show here.
  function rowAtGrid(clientY) {
    var rect = gridCanvas.getBoundingClientRect();
    var ri = Math.floor((clientY - rect.top) / ROW_H);
    return ri;
  }
  // Cell hover: the row (strains, and their species breakdown when the
  // payload has species data) plus the column's family, presence in this
  // row, and location in the example strain.
  gridCanvas.addEventListener("mousemove", function (e) {
    var isl = ISLANDS[state.selected];
    var hap = currentHaps[rowAtGrid(e.clientY)];
    if (!isl || !hap) { tipEl.style.display = "none"; return; }
    var ci = colAtGlyph(e.clientX, isl);
    tipEl.textContent = "";
    var rowName = Object.keys(SPECIES).length ? speciesBandLabel(hap, SPECIES)
      : (hap.count === 1 ? hap.strains[0] : hap.strains[0] + " +" + (hap.count - 1) + " more");
    tipEl.appendChild(el("div", "tip-id", rowName + " (" + hap.count +
      (hap.count === 1 ? " strain)" : " strains)")));
    if (Object.keys(SPECIES).length) {
      var counts = speciesCounts(hap, SPECIES);
      Object.keys(counts).sort().forEach(function (sp) {
        tipEl.appendChild(el("div", null, (sp || "Unknown species") + ": " + counts[sp] +
          (counts[sp] === 1 ? " strain" : " strains")));
      });
    }
    if (ci >= 0) {
      var on = hap.pattern.charAt(ci) === "1";
      tipEl.appendChild(el("div", null, (on ? "Present" : "Absent") + " in this row"));
      appendColumnInfo(isl, ci);
    }
    positionTip(e);
  });
  gridCanvas.addEventListener("mouseleave", function () { tipEl.style.display = "none"; });

  // ---- legend ---------------------------------------------------------------
  // One swatch per DISTINCT class actually present among this island's
  // columns (not a fixed five-entry key), since the strip below is coloured
  // per family -- a legend fixed to `dominant_class` would silently mislabel
  // every other-coloured column.
  function renderLegend(isl) {
    var legend = document.getElementById("isv-legend");
    legend.textContent = "";
    var seen = {};
    isl.families.forEach(function (fam, i) {
      var key = familyClassAt(isl, i);
      if (seen[key]) return;
      seen[key] = true;
      var item = el("span", "isv-legend-item");
      var sw = el("span", "isv-swatch");
      sw.style.background = classColor(key);
      item.appendChild(sw);
      item.appendChild(el("span", null, classLabel(key)));
      legend.appendChild(item);
    });

    var absent = el("span", "isv-legend-item");
    var swAbsent = el("span", "isv-swatch");
    swAbsent.style.background = css("--grid");
    absent.appendChild(swAbsent);
    absent.appendChild(el("span", null, "Absent"));
    legend.appendChild(absent);
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
        (isl.locus_start >= 0 ? ":" + isl.locus_start + "-" + isl.locus_end : "") : "") +
      ". Columns: consecutive genes along this contig in " + (isl.example_strain || "the example strain") +
      ", left to right by position (a run of non-core families with no core gene between them). " +
      "Rows: strains with the same presence pattern; a filled cell means the family is present " +
      "in that strain somewhere in the genome, not necessarily at this locus.";

    renderLegend(isl);

    var haps = sortedHaplotypes(isl);
    currentHaps = haps;
    var w = totalWidth(isl);
    drawGlyphStrip(isl, w);
    drawGrid(isl, haps, w);
  }

  window.onSkinChange = function () {
    if (state.selected >= 0) renderMain();
  };

  // ---- shared "why is nothing drawn" reasons list -----------------------
  // A located island is missing from the view for one of two independent
  // reasons -- excluded by --min_strains, or truncated past --top_islands --
  // and both the summary-note (always visible) and the empty-state box
  // (visible only when ISLANDS is empty) need to say so without disagreeing.
  // This is the SAME array both call sites build a sentence from, so they
  // are structurally unable to drift back apart the way they did before:
  // the empty-state box used to name only the excluded count, silently
  // dropping the truncated one whenever --top_islands did the cutting.
  function absenceReasons() {
    var reasons = [];
    if (DATA.n_islands_excluded) {
      reasons.push(
        DATA.n_islands_excluded + " located island" + (DATA.n_islands_excluded === 1 ? "" : "s") +
        " excluded for being carried by too few strains to show a meaningful presence " +
        "pattern (a single-strain island is one filled row with no breakpoint in it)"
      );
    }
    if (DATA.n_islands_truncated) {
      reasons.push(
        DATA.n_islands_truncated + " qualifying island" + (DATA.n_islands_truncated === 1 ? "" : "s") +
        " truncated past the --top_islands limit"
      );
    }
    return reasons;
  }

  // ---- empty state (no islands at all in this payload) -----------------------
  function renderEmptyState() {
    var explorer = document.getElementById("isv-explorer");
    var emptyState = document.getElementById("isv-empty-state");
    if (!ISLANDS.length) {
      explorer.classList.add("hidden");
      emptyState.classList.remove("hidden");
      var text;
      if (!DATA.n_islands_total) {
        text = "No accessory islands were located for this run.";
      } else {
        var reasons = absenceReasons();
        text = DATA.n_islands_total + " accessory island" + (DATA.n_islands_total === 1 ? "" : "s") +
          " were located, but none remain to draw: " + reasons.join("; ") + ".";
      }
      document.getElementById("isv-empty-text").textContent = text;
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

  // Issue #119: only offer the "species" row sort when the payload actually
  // carries species data -- mirrors novelties.html's f-category filter
  // (lib/report_template.py), which appends its <option>s only when
  // DATA.novelty_categories is non-empty rather than always shipping a
  // static option that then does nothing on a run with no data for it.
  if (Object.keys(SPECIES).length) {
    document.getElementById("f-row-sort").appendChild(
      new Option("Sort rows: species", "species"));
  }

  // ---- init ---------------------------------------------------------------
  document.getElementById("title").textContent = DATA.project + " — island synteny";
  document.getElementById("subtitle").textContent =
    ISLANDS.length + " of " + DATA.n_islands_total + " located accessory islands shown, in locus order.";
  document.title = DATA.project + " — NovInvenio island synteny";

  (function () {
    var reasons = absenceReasons();
    document.getElementById("summary-note").textContent = reasons.length
      ? reasons.join("; ") + "."
      : "Every located accessory island met the minimum strain-count filter for this view.";
  })();
  document.getElementById("t-shown").textContent = ISLANDS.length.toLocaleString();
  document.getElementById("t-total").textContent = (DATA.n_islands_total || 0).toLocaleString();
  document.getElementById("t-excluded").textContent = (DATA.n_islands_excluded || 0).toLocaleString();
  document.getElementById("t-truncated").textContent = (DATA.n_islands_truncated || 0).toLocaleString();
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
