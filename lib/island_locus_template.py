"""Locus view fragments for island_synteny.html (spec
docs/superpowers/specs/2026-09-24-island-locus-view-design.md, section 7).

lib/island_synteny_template.py inserts three fragments:
  LOCUS_VIEW_CSS   into the page <style>
  LOCUS_VIEW_HTML  before the island view, inside .wrap
  LOCUS_VIEW_JS    inside the page's one IIFE, after the island-view code,
                   so it reuses el(), css(), palette(), classColor(),
                   classLabel(), sizeCanvas(), ellipsize(), positionTip(),
                   tipEl, appendLocation(), strainPopupLines(),
                   speciesCounts(), haplotypeSpecies() and SPECIES.

The page shows the locus view by default when the payload has `loci`; the
island (presence) view stays one click away (Ruling R1). Colours come from
skin tokens only: --series-1 = in place (protein-search presence),
--series-2 = TBLASTN rescue, --grid = absent, --warn = contig break.
Species are told apart by position and label, never by hue.
"""

LOCUS_VIEW_CSS = r"""
  .lv-switch { display: inline-flex; gap: 0; margin: 0 0 14px; border: 1px solid var(--border); border-radius: 8px; overflow: hidden; }
  .lv-switch button { border: none; border-radius: 0; padding: 6px 14px; background: var(--surface-1); }
  .lv-switch button[aria-pressed="true"] { background: var(--wash); box-shadow: inset 0 -2px 0 var(--series-1); font-weight: 600; }
  .lv-badge { display: inline-block; margin-left: 6px; padding: 1px 7px; border-radius: 999px; font-size: 10.5px; border: 1px solid var(--warn); color: var(--text-primary); }
  .lv-track-wrap { border-bottom: 1px solid var(--border); }
  canvas.lv-track, canvas.lv-head { display: block; }
  .lv-confound { margin: 0 0 12px; font-size: 12px; color: var(--text-secondary); }
"""

LOCUS_VIEW_HTML = r"""
  <p class="lv-confound hidden" id="lv-confound"></p>
  <div class="lv-switch hidden" id="lv-switch" role="group" aria-label="Choose view">
    <button type="button" id="lv-btn-loci" aria-pressed="true">Loci (flank-anchored)</button>
    <button type="button" id="lv-btn-islands" aria-pressed="false">Islands (presence)</button>
  </div>
  <div id="locus-view" class="hidden">
    <div class="isv-explorer">
      <aside class="isv-sidebar">
        <div class="isv-sidebar-controls">
          <input type="search" id="lv-search" placeholder="Search locus, strain or family ID…" aria-label="Search loci">
          <select id="lv-sort" aria-label="Sort loci by">
            <option value="informative">Sort: informative polymorphism</option>
            <option value="strains">Sort: strains with the locus</option>
            <option value="size">Sort: size (families)</option>
            <option value="name">Sort: locus ID</option>
          </select>
          <span class="count" id="lv-count" role="status" aria-live="polite"></span>
        </div>
        <div class="isv-list" id="lv-list" role="listbox" aria-label="Accessory loci"></div>
      </aside>
      <section class="card" id="lv-main" style="padding:14px">
        <div class="isv-main-head">
          <h2 id="lv-title"></h2>
          <span class="lv-badge hidden" id="lv-tier"></span>
        </div>
        <p class="isv-main-note" id="lv-note"></p>
        <div class="isv-legend" id="lv-legend"></div>
        <div class="isv-hscroll" id="lv-hscroll">
          <div class="lv-track-wrap"><canvas class="lv-track" id="lv-track"></canvas></div>
          <canvas class="lv-head" id="lv-head"></canvas>
          <div class="isv-vscroll"><canvas class="isv-grid" id="lv-grid"></canvas></div>
        </div>
      </section>
    </div>
  </div>
"""

LOCUS_VIEW_JS = r"""
  // ==== locus view ==========================================================
  var LOCI = DATA.loci || [];
  var LMETA = DATA.locus_meta || {};
  var LPARAMS = LMETA.locus_params || {};
  var lstate = { search: "", sort: "informative", selected: LOCI.length ? 0 : -1 };
  var lSidebar = [];
  var lRows = [];
  var lSlots = [];
  var L_ROW_H = 20, L_CELL_W = 22, L_TRACK_H = 64, L_GUTTER = 190;
  var lTrack = document.getElementById("lv-track");
  var lHead = document.getElementById("lv-head");
  var lGrid = document.getElementById("lv-grid");
  var ltctx = lTrack.getContext("2d");
  var lhctx = lHead.getContext("2d");
  var lgctx = lGrid.getContext("2d");

  // Pure helpers: no DOM, no module state, so tests can run them under node.
  function locusStateStyle(code) {
    if (code === "1") return { token: "--series-1", alpha: 1, hatch: false };
    if (code === "2") return { token: "--series-1", alpha: 0.35, hatch: false };
    if (code === "3") return { token: "--series-2", alpha: 1, hatch: true };
    if (code === "4") return { token: "--series-2", alpha: 0.35, hatch: true };
    if (code === "5") return { token: "--warn", alpha: 0.55, hatch: true };
    if (code === "6") return { token: "--text-secondary", alpha: 0.45, hatch: true };
    if (code === "7") return { token: "--text-primary", alpha: 0.7, hatch: false };
    return { token: "--grid", alpha: 1, hatch: false };
  }
  function locusStateLabel(code) {
    var labels = { "0": "absent", "1": "in place", "2": "elsewhere",
      "3": "rescue, in place (TBLASTN hit, no annotated gene)",
      "4": "rescue, elsewhere (TBLASTN hit, no annotated gene)", "5": "contig break",
      "6": "absent, DNA present (gene-model or annotation difference)",
      "7": "absent, DNA absent" };
    return labels[code] || "unknown";
  }
  function locusClassLabel(cls) {
    var labels = { full: "full locus", partial: "partial", empty: "empty site",
      model_difference: "model difference", uninformative: "uninformative" };
    return labels[cls] || cls;
  }
  function locusSortedRows(rows, speciesMap) {
    var order = ["full", "partial", "empty", "model_difference", "uninformative"];
    return rows.slice().sort(function (a, b) {
      var ca = order.indexOf(a.row_class), cb = order.indexOf(b.row_class);
      if (ca !== cb) return ca - cb;
      var sa = haplotypeSpecies(a, speciesMap) || "￿";
      var sb = haplotypeSpecies(b, speciesMap) || "￿";
      if (sa !== sb) return sa < sb ? -1 : 1;
      if (b.count !== a.count) return b.count - a.count;
      return a.codes < b.codes ? -1 : a.codes > b.codes ? 1 : 0;
    });
  }
  function locusSidebarStats(locus, params) {
    var c = locus.counts;
    var dnaOn = params && params.dna_check === true && locus.dna;
    var emptyCount = dnaOn ? locus.dna.empty_confirmed : c.empty;
    var text = locus.size + " families · " + locus.n_variants + " variants · " +
      "empty " + emptyCount + " · full " + c.full + " · partial " + c.partial;
    if (c.model_difference !== undefined) text += " · model difference " + c.model_difference;
    text += " · uninformative " + c.uninformative;
    if (dnaOn && locus.dna.unchecked > 0) text += " · empty site, not checked " + locus.dna.unchecked;
    return text;
  }
  function locusDnaNote(params, dna) {
    if (params.dna_check === true) {
      // Ruling R22: only DNA-confirmed empty sites count as "empty site";
      // unchecked strains (no DNA call at all) must not read as confirmed.
      var opening = "Empty site is DNA-confirmed:";
      if (dna && dna.unchecked > 0) {
        var m = dna.empty_confirmed, n = dna.unchecked;
        opening = "Empty site is DNA-confirmed for " + m + " of " + (m + n) +
          " checked strains; " + n + " had no DNA call.";
      }
      return opening + " blastn (megablast) of the exemplar's locus DNA " +
        "against the strain's DNA from its left to its right flank gene; a gene is DNA present at >= " +
        params.dna_min_id + "% identity over >= " + params.dna_min_cov + "% of its length. " +
        "Hatched grey = absent, DNA present (model difference, not a deletion); " +
        "dark = DNA absent (confirmed).";
    }
    return "Empty site is not DNA-confirmed (the DNA presence check did not run), so it can " +
      "be a gene-model or annotation difference.";
  }
  function locusSidebarOrder(loci, idxs, key) {
    var cmp;
    if (key === "strains") cmp = function (a, b) { return loci[b].n_carriers - loci[a].n_carriers || a - b; };
    else if (key === "size") cmp = function (a, b) { return loci[b].size - loci[a].size || a - b; };
    else if (key === "name") cmp = function (a, b) { return loci[a].locus_id < loci[b].locus_id ? -1 : loci[a].locus_id > loci[b].locus_id ? 1 : 0; };
    else cmp = function (a, b) { return loci[b].informative_score - loci[a].informative_score || loci[b].n_carriers - loci[a].n_carriers || a - b; };
    return idxs.slice().sort(cmp);
  }
  function locusSpeciesList(locus) {
    return Object.keys(locus.counts_by_species || {}).sort();
  }
  function breakpointBars(bp, speciesList) {
    var bars = speciesList.map(function (sp) {
      return { kind: "indel", species: sp, n: (bp.indel && bp.indel[sp]) || 0 };
    });
    bars.push({ kind: "contig_break", species: "", n: bp.contig_break || 0 });
    return bars;
  }
  function cellReasonLines(locus, row, ci, k) {
    var code = row.codes.charAt(ci);
    var lines = ["State: " + locusStateLabel(code)];
    var who = row.strains[0] + (row.count > 1 ? " (first of " + row.count + " strains)" : "");
    if (!row.rep) { lines.push("Position detail not stored for this row"); return lines; }
    var d = row.rep.d[ci];
    var dnaLine = code === "6"
      ? "The exemplar gene's DNA is at this site in " + row.strains[0] + " (blastn): a gene-model or annotation difference, not a deletion"
      : (code === "7" ? "The site lacks the exemplar gene's DNA in " + row.strains[0] + " (blastn)" : "");
    if (!d) {
      lines.push(code === "5"
        ? "No copy at the locus in " + who + "; the nearest placed column is within " + k + " genes of a contig end, so this absence is not evidence"
        : (code === "2" ? "Present in " + who + " but no position recorded" : "No copy in " + who));
      if (dnaLine) lines.push(dnaLine);
      return lines;
    }
    lines.push("In " + who + ": " + row.rep.c[d[0]] + ", gene rank " + d[1]);
    if (d.length > 2) {
      lines.push("Neighbour " + locus.families[d[2]] + " at rank " + (d[1] + d[3]) +
        " (" + (d[3] > 0 ? "+" : "") + d[3] + ")");
    } else {
      lines.push("No other column within " + k + " genes on this contig");
    }
    if (dnaLine) lines.push(dnaLine);
    return lines;
  }
  function locusTitle(locus) {
    var span = locus.exemplar_span;
    return locus.exemplar + ":" + locus.exemplar_contig +
      (span ? ":" + span.start + "-" + span.end : "");
  }
  function locusTierText(tier, params) {
    if (tier === "short_flanks") return "short flanks (" + params.flank_min + " genes per side)";
    if (tier === "contig_end") return "exemplar at contig end";
    return "";
  }

  function lcolX(i) { return L_GUTTER + i * L_CELL_W; }
  function lTotalWidth(locus) { return lcolX(locus.families.length) + 8; }
  function isFlank(locus, i) { return i < locus.n_left || i >= locus.n_left + locus.n_locus; }

  function locusHaystack(locus) {
    return (locus.locus_id + " " + locus.exemplar + " " + locus.families.join(" ")).toLowerCase();
  }
  function applyLocusFilter() {
    var q = lstate.search.trim().toLowerCase();
    var idxs = [];
    for (var i = 0; i < LOCI.length; i++) {
      if (q && locusHaystack(LOCI[i]).indexOf(q) === -1) continue;
      idxs.push(i);
    }
    lSidebar = locusSidebarOrder(LOCI, idxs, lstate.sort);
  }
  function renderLocusSidebar() {
    var list = document.getElementById("lv-list");
    list.textContent = "";
    if (!lSidebar.length) {
      list.appendChild(el("p", "isv-empty-list", "No loci match the current search."));
    }
    lSidebar.forEach(function (idx) {
      var locus = LOCI[idx];
      var btn = el("button", "isv-item");
      btn.type = "button";
      btn.setAttribute("role", "option");
      btn.setAttribute("aria-selected", idx === lstate.selected ? "true" : "false");
      if (idx === lstate.selected) btn.classList.add("sel");
      btn.appendChild(el("div", "isv-item-id", locus.locus_id));
      btn.appendChild(el("div", "isv-item-stats", locusSidebarStats(locus, LPARAMS)));
      var chip = el("span", "isv-chip", classLabel(locus.dominant_class));
      chip.style.borderColor = classColor(locus.dominant_class);
      chip.style.color = classColor(locus.dominant_class);
      btn.appendChild(chip);
      var tier = locusTierText(locus.tier, LPARAMS);
      if (tier) btn.appendChild(el("span", "lv-badge", tier));
      btn.addEventListener("click", function () { selectLocus(idx); });
      list.appendChild(btn);
    });
    document.getElementById("lv-count").textContent =
      lSidebar.length.toLocaleString() + " of " + LOCI.length.toLocaleString() + " shown";
  }
  function selectLocus(idx) {
    lstate.selected = idx;
    renderLocusSidebar();
    renderLocusMain();
  }

  function lRowLabel(row) {
    return (row.count === 1 ? row.strains[0] : row.strains[0] + " +" + (row.count - 1));
  }
  function locusGridSlots(rows) {
    var slots = [];
    var prevClass = null;
    rows.forEach(function (row) {
      if (row.row_class !== prevClass) {
        var count = rows.reduce(function (sum, r) {
          return r.row_class === row.row_class ? sum + r.count : sum;
        }, 0);
        slots.push({ header: true, cls: row.row_class, count: count });
        prevClass = row.row_class;
      }
      slots.push({ header: false, row: row });
    });
    return slots;
  }
  function locusRowAtSlot(slots, index) {
    var slot = slots[index];
    return slot && !slot.header ? slot.row : null;
  }
  function locusGutter(rows) {
    lgctx.font = ROW_LABEL_FONT;
    var maxW = 0;
    rows.forEach(function (r) { maxW = Math.max(maxW, lgctx.measureText(lRowLabel(r)).width); });
    return Math.min(GUTTER_MAX, Math.max(GUTTER_MIN, maxW + 60 + 12));
  }

  function hatch(ctx, x, y, w, h, color) {
    ctx.save();
    ctx.beginPath();
    ctx.rect(x, y, w, h);
    ctx.clip();
    ctx.strokeStyle = color;
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (var off = -h; off < w; off += 5) {
      ctx.moveTo(x + off, y + h);
      ctx.lineTo(x + off + h, y);
    }
    ctx.stroke();
    ctx.restore();
  }

  function drawLocusTrack(locus, totalW) {
    var P = palette();
    sizeCanvas(lTrack, ltctx, totalW, L_TRACK_H);
    ltctx.clearRect(0, 0, totalW, L_TRACK_H);
    ltctx.fillStyle = P.surface;
    ltctx.fillRect(0, 0, totalW, L_TRACK_H);
    var species = locusSpeciesList(locus);
    var maxN = 1;
    (locus.breakpoints || []).forEach(function (bp) {
      breakpointBars(bp, species).forEach(function (b) { maxN = Math.max(maxN, b.n); });
    });
    var base = L_TRACK_H - 14;
    ltctx.font = ROW_LABEL_FONT;
    ltctx.fillStyle = P.secondary;
    ltctx.textBaseline = "middle";
    ltctx.fillText("breakpoints (max " + maxN + ")", 8, 10);
    ltctx.fillText("species: " + (species.map(function (s) { return s || "unknown"; }).join(" | ")), 8, L_TRACK_H - 6);
    var nBars = species.length + 1;
    var barW = Math.max(2, Math.floor((L_CELL_W - 4) / nBars));
    (locus.breakpoints || []).forEach(function (bp) {
      var bars = breakpointBars(bp, species);
      var x0 = lcolX(bp.b) - (barW * nBars) / 2;
      bars.forEach(function (b, i) {
        if (!b.n) return;
        var h = Math.max(1, Math.round((base - 16) * b.n / maxN));
        ltctx.fillStyle = b.kind === "contig_break" ? css("--warn") : P.secondary;
        ltctx.fillRect(x0 + i * barW, base - h, barW - 1, h);
      });
    });
    ltctx.fillStyle = P.axis;
    ltctx.fillRect(0, base, totalW, 1);
  }

  function drawLocusHead(locus, totalW) {
    var P = palette();
    lhctx.font = LABEL_FONT;
    var maxW = 0;
    locus.families.forEach(function (f) { maxW = Math.max(maxW, Math.min(LABEL_MAX_W, lhctx.measureText(f).width)); });
    var H = Math.ceil(maxW * Math.sin(LABEL_ANGLE)) + GLYPH_H + 32;
    var W = totalW + Math.ceil(maxW * Math.cos(LABEL_ANGLE)) + 12;
    var glyphY = H - GLYPH_H - 10;
    sizeCanvas(lHead, lhctx, W, H);
    lhctx.clearRect(0, 0, W, H);
    lhctx.fillStyle = P.surface;
    lhctx.fillRect(0, 0, W, H);
    locus.families.forEach(function (fam, i) {
      var x = lcolX(i);
      lhctx.fillStyle = classColor((locus.family_classes && locus.family_classes[i]) || "unannotated");
      lhctx.fillRect(x + 3, glyphY, L_CELL_W - 6, GLYPH_H);
      if (isFlank(locus, i)) {
        lhctx.fillStyle = P.axis;
        lhctx.fillRect(x + 1, H - 6, L_CELL_W - 2, 3);
      }
      lhctx.save();
      lhctx.translate(x + L_CELL_W / 2, glyphY - 6);
      lhctx.rotate(-LABEL_ANGLE);
      lhctx.fillStyle = isFlank(locus, i) ? P.secondary : P.primary;
      lhctx.font = LABEL_FONT;
      lhctx.textAlign = "left";
      lhctx.textBaseline = "middle";
      lhctx.fillText(ellipsize(lhctx, fam, LABEL_MAX_W), 0, 0);
      lhctx.restore();
    });
    lhctx.font = ROW_LABEL_FONT;
    lhctx.fillStyle = P.secondary;
    lhctx.textBaseline = "middle";
    lhctx.fillText("flank = anchor (bar under column)", 8, H - 16);
  }

  function drawLocusGrid(locus, slots, totalW) {
    var P = palette();
    var h = Math.max(L_ROW_H, slots.length * L_ROW_H);
    sizeCanvas(lGrid, lgctx, totalW, h);
    lgctx.clearRect(0, 0, totalW, h);
    lgctx.fillStyle = P.surface;
    lgctx.fillRect(0, 0, totalW, h);
    slots.forEach(function (slot, si) {
      var y = si * L_ROW_H;
      lgctx.textBaseline = "middle";
      lgctx.textAlign = "left";
      if (slot.header) {
        lgctx.fillStyle = P.axis;
        lgctx.fillRect(0, y, totalW, 1);
        lgctx.font = "600 11px ui-monospace, SFMono-Regular, Menlo, monospace";
        lgctx.fillStyle = css("--text-primary");
        lgctx.fillText(locusClassLabel(slot.cls) + " · " + slot.count.toLocaleString() +
          " strain" + (slot.count === 1 ? "" : "s"), 8, y + L_ROW_H / 2);
        return;
      }
      var row = slot.row;
      lgctx.font = "600 11px ui-monospace, SFMono-Regular, Menlo, monospace";
      lgctx.fillStyle = P.primary;
      lgctx.fillText("×" + row.count, 8, y + L_ROW_H / 2);
      lgctx.font = ROW_LABEL_FONT;
      lgctx.fillStyle = P.secondary;
      lgctx.fillText(ellipsize(lgctx, lRowLabel(row), L_GUTTER - 60), 50, y + L_ROW_H / 2);
      for (var ci = 0; ci < row.codes.length; ci++) {
        var st = locusStateStyle(row.codes.charAt(ci));
        var x = lcolX(ci) + 1, cw = L_CELL_W - 2, ch = L_ROW_H - 2;
        lgctx.fillStyle = P.grid;
        lgctx.fillRect(x, y + 1, cw, ch);
        lgctx.globalAlpha = st.alpha;
        lgctx.fillStyle = css(st.token);
        lgctx.fillRect(x, y + 1, cw, ch);
        lgctx.globalAlpha = 1;
        if (st.hatch) hatch(lgctx, x, y + 1, cw, ch, P.surface);
      }
    });
    lgctx.fillStyle = P.axis;
    lgctx.fillRect(lcolX(locus.n_left), 0, 1, h);
    lgctx.fillRect(lcolX(locus.n_left + locus.n_locus), 0, 1, h);
  }

  function legendCodes(dnaCheck) {
    return dnaCheck === true ? ["1", "2", "3", "4", "0", "5", "6", "7"] : ["1", "2", "3", "4", "0", "5"];
  }
  function legendLabel(code, dnaCheck) {
    // With the DNA check on, code 0 (no call at all) and code 7 (checked,
    // DNA absent) share the same swatch (--grid, no hatch); only the label
    // tells them apart, so 0 must not read as if it were itself confirmed.
    if (code === "0" && dnaCheck === true) return "absent (not checked)";
    return locusStateLabel(code).split(" (")[0];
  }
  function renderLocusLegend() {
    var legend = document.getElementById("lv-legend");
    legend.textContent = "";
    var dnaCheck = LPARAMS.dna_check === true;
    var codes = legendCodes(dnaCheck);
    codes.forEach(function (code) {
      var st = locusStateStyle(code);
      var item = el("span", "isv-legend-item");
      var sw = el("span", "isv-swatch");
      sw.style.background = css(st.token);
      sw.style.opacity = String(st.alpha);
      if (st.hatch) sw.style.backgroundImage = "repeating-linear-gradient(45deg, transparent 0 3px, " + css("--surface-1") + " 3px 4px)";
      item.appendChild(sw);
      item.appendChild(el("span", null, legendLabel(code, dnaCheck)));
      legend.appendChild(item);
    });
    var bp = el("span", "isv-legend-item",
      "Bars above the columns: flank-intact strains changing between in place and " +
      (LPARAMS.dna_check === true ? "DNA absent" : "absent") +
      ", one bar per species (left to right as named); the last bar counts contig breaks");
    legend.appendChild(bp);
  }

  function renderLocusMain() {
    var locus = LOCI[lstate.selected];
    if (!locus) return;
    document.getElementById("lv-title").textContent = locusTitle(locus);
    var tierEl = document.getElementById("lv-tier");
    var tierText = locusTierText(locus.tier, LPARAMS);
    tierEl.textContent = tierText;
    tierEl.classList.toggle("hidden", !tierText);
    var c = locus.counts;
    var dnaOn = LPARAMS.dna_check === true && locus.dna;
    var emptyCount = dnaOn ? locus.dna.empty_confirmed : c.empty;
    document.getElementById("lv-note").textContent =
      "Exemplar " + locus.exemplar + ": carries the locus's largest variant with at least " +
      LPARAMS.flank + " genes on both sides (else " + LPARAMS.flank_min + "), ties by N50 then name. " +
      locus.n_left + " + " + locus.n_locus + " + " + locus.n_right +
      " columns (flank + locus + flank) in the exemplar's gene order. A cell is in place when " +
      "the strain has another column's gene within k = " + LPARAMS.k + " genes on the same contig. " +
      "Rows: strains with identical states, grouped by row class, then species. " +
      "Full " + c.full + ", partial " + c.partial + ", empty site " + emptyCount +
      (c.model_difference !== undefined ? ", model difference " + c.model_difference : "") +
      ", uninformative " + c.uninformative +
      (dnaOn && locus.dna.unchecked > 0 ? ", empty site, not checked " + locus.dna.unchecked : "") +
      " strains. " + locusDnaNote(LPARAMS, locus.dna);
    renderLocusLegend();
    lRows = locusSortedRows(locus.rows, SPECIES);
    lSlots = locusGridSlots(lRows);
    L_GUTTER = locusGutter(lRows);
    var w = lTotalWidth(locus);
    drawLocusTrack(locus, w);
    drawLocusHead(locus, w);
    drawLocusGrid(locus, lSlots, w);
    if (typeof renderClinkerPanel === "function") renderClinkerPanel(locus);
  }

  function lColAt(canvas, clientX, locus) {
    var x = clientX - canvas.getBoundingClientRect().left;
    for (var i = 0; i < locus.families.length; i++) {
      if (x >= lcolX(i) && x < lcolX(i) + L_CELL_W) return i;
    }
    return -1;
  }
  function appendLocusColumn(locus, ci) {
    var role = ci < locus.n_left ? "left flank (anchor)"
      : (ci >= locus.n_left + locus.n_locus ? "right flank (anchor)" : "locus gene");
    tipEl.appendChild(el("div", "tip-id", locus.families[ci]));
    tipEl.appendChild(el("div", null, "Column " + (ci + 1) + " of " + locus.families.length + ", " + role));
    tipEl.appendChild(el("div", null, "Frequency bin: " + ((locus.family_bins && locus.family_bins[ci]) || "unknown")));
    appendLocation({ family_locations: locus.family_locations, example_strain: locus.exemplar }, ci);
    tipEl.appendChild(el("div", null, classLabel((locus.family_classes && locus.family_classes[ci]) || "unannotated")));
    var doms = (locus.family_domains && locus.family_domains[ci]) || "";
    tipEl.appendChild(el("div", null, doms ? "Domains: " + doms : "No annotated Pfam domain"));
  }
  lHead.addEventListener("mousemove", function (e) {
    var locus = LOCI[lstate.selected];
    if (!locus) return;
    var ci = lColAt(lHead, e.clientX, locus);
    if (ci < 0) { tipEl.style.display = "none"; return; }
    tipEl.textContent = "";
    appendLocusColumn(locus, ci);
    positionTip(e);
  });
  lHead.addEventListener("mouseleave", function () { tipEl.style.display = "none"; });
  lGrid.addEventListener("mousemove", function (e) {
    var locus = LOCI[lstate.selected];
    var row = locusRowAtSlot(lSlots, Math.floor((e.clientY - lGrid.getBoundingClientRect().top) / L_ROW_H));
    if (!locus || !row) { tipEl.style.display = "none"; return; }
    tipEl.textContent = "";
    tipEl.appendChild(el("div", "tip-id", row.count + (row.count === 1 ? " strain" : " strains") +
      " · " + locusClassLabel(row.row_class)));
    strainPopupLines(row).slice(1).forEach(function (line) { tipEl.appendChild(el("div", null, line)); });
    if (Object.keys(SPECIES).length) {
      var counts = speciesCounts(row, SPECIES);
      Object.keys(counts).sort().forEach(function (sp) {
        tipEl.appendChild(el("div", null, (sp || "Unknown species") + ": " + counts[sp]));
      });
    }
    var ci = lColAt(lGrid, e.clientX, locus);
    if (ci >= 0) {
      tipEl.appendChild(el("div", "tip-id", locus.families[ci]));
      cellReasonLines(locus, row, ci, LPARAMS.k).forEach(function (line) {
        tipEl.appendChild(el("div", null, line));
      });
    }
    positionTip(e);
  });
  lGrid.addEventListener("mouseleave", function () { tipEl.style.display = "none"; });

  function setView(view) {
    var loci = view === "loci";
    document.getElementById("locus-view").classList.toggle("hidden", !loci);
    document.getElementById("island-view").classList.toggle("hidden", loci);
    document.getElementById("lv-btn-loci").setAttribute("aria-pressed", loci ? "true" : "false");
    document.getElementById("lv-btn-islands").setAttribute("aria-pressed", loci ? "false" : "true");
    if (loci) renderLocusMain();
    else if (state.selected >= 0) renderMain();
  }
  document.getElementById("lv-btn-loci").addEventListener("click", function () { setView("loci"); });
  document.getElementById("lv-btn-islands").addEventListener("click", function () { setView("islands"); });
  document.getElementById("lv-search").addEventListener("input", function (e) {
    lstate.search = e.target.value;
    applyLocusFilter();
    renderLocusSidebar();
  });
  document.getElementById("lv-sort").addEventListener("change", function (e) {
    lstate.sort = e.target.value;
    applyLocusFilter();
    renderLocusSidebar();
  });
  var islandSkinChange = window.onSkinChange;
  window.onSkinChange = function () {
    if (islandSkinChange) islandSkinChange();
    if (LOCI.length) renderLocusMain();
  };
  if (DATA.assembly_confound) {
    var conf = document.getElementById("lv-confound");
    conf.textContent = "Accessory content is confounded with assembly quality on this run " +
      "(assembly_quality_confound triggered), so uninformative rows are expected.";
    conf.classList.remove("hidden");
  }
  if (LOCI.length) {
    document.getElementById("lv-switch").classList.remove("hidden");
    applyLocusFilter();
    renderLocusSidebar();
    setView("loci");
  }
"""
