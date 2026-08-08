"""
Self-contained HTML/CSS/JS template for the NovInvenio LOSSES report —
candidate lineage-specific gene losses (present in the outgroup, absent from
the ingroup).

bin/make_losses_report.py substitutes two tokens:
  __PROJECT_TITLE__   plain-text project name (HTML-escaped by the caller)
  /*__PAYLOAD__*/     the JSON payload from lib/report_data.build_losses_payload()

Same file:// / no-network / textContent-only constraints as
lib/report_template.py. Structurally this is lib/core_report_template.py's
twin — a single sortable/filterable table built on lib/report_common.py's
shared fragments — with loss-specific columns (outgroup conservation
fraction, TBLASTN-vs-ingroup-genome flag) in place of the presence fraction.
"""

from report_common import (
    BASE_PAGE_CSS,
    DOWNLOAD_JS,
    EL_HELPER_JS,
    LINKOUT_HELPERS_JS,
    THEME_TOGGLE_JS,
    THEME_VARS_CSS,
)

LOSSES_HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__PROJECT_TITLE__ — NovInvenio candidate losses</title>
<style>
""" + THEME_VARS_CSS + BASE_PAGE_CSS + r"""
  .badge {
    display: inline-block; padding: 1px 6px; border-radius: 999px; font-size: 10.5px;
    border: 1px solid var(--border); color: var(--text-secondary);
  }
  .badge.warn { color: #a15c00; border-color: #a15c00; }
  @media (prefers-color-scheme: dark) {
    :root:where(:not([data-theme="light"])) .badge.warn { color: #e0a030; border-color: #e0a030; }
  }
  :root[data-theme="dark"] .badge.warn { color: #e0a030; border-color: #e0a030; }
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

  <section class="card">
    <h2 class="card-title">Run summary</h2>
    <p class="card-note" id="summary-note"></p>
    <div class="tiles">
      <div>
        <div class="tile-label">Loss candidates</div>
        <div class="tile-value" id="t-total"></div>
      </div>
      <div>
        <div class="tile-label">Outgroup species</div>
        <div class="tile-value" id="t-out"></div>
      </div>
      <div>
        <div class="tile-label">Gene families</div>
        <div class="tile-value" id="t-fam"></div>
      </div>
      <div>
        <div class="tile-label">Flagged by ingroup TBLASTN</div>
        <div class="tile-value" id="t-flagged"></div>
      </div>
    </div>
  </section>

  <div class="filters" role="group" aria-label="Filter loss candidates">
    <input type="search" id="f-search" placeholder="Search ID, gene, product, Pfam…" aria-label="Search proteins">
    <select id="f-src" aria-label="Source outgroup proteome"></select>
    <select id="f-fsrc" aria-label="Annotation source"></select>
    <select id="f-family" aria-label="Gene family"></select>
    <label class="check" style="display:inline-flex;align-items:center;gap:6px;font-size:13px;cursor:pointer">
      <input type="checkbox" id="f-notb"> No ingroup TBLASTN hit
    </label>
    <select id="f-sort" aria-label="Sort by">
      <option value="priority">Sort: priority (clean loss, broad in outgroup, no TBLASTN hit first)</option>
      <option value="breadth">Sort: outgroup family breadth</option>
      <option value="frac">Sort: outgroup presence fraction</option>
      <option value="id">Sort: protein ID</option>
      <option value="src">Sort: source proteome</option>
      <option value="pos">Sort: genomic position (chrom, start)</option>
    </select>
    <button id="f-reset" type="button">Reset</button>
    <div class="spacer"></div>
    <span class="count" id="count"></span>
    <button id="dl-tsv" type="button">Download TSV</button>
  </div>

  <div class="explorer">
    <section class="card" style="padding:14px">
      <div class="tbl-scroll" id="tbl-scroll">
        <table class="data">
          <caption class="sr-only">Candidate gene losses matching the current filters</caption>
          <thead id="tbl-head"></thead>
          <tbody id="tbl-body"></tbody>
        </table>
      </div>
      <p class="empty hidden" id="tbl-empty">No proteins match the current filters.</p>
    </section>

    <aside class="card detail" id="detail" aria-live="polite">
      <p class="placeholder">Select a protein in the table to see its annotation and database links.</p>
    </aside>
  </div>
</div>

<script type="application/json" id="payload">/*__PAYLOAD__*/</script>
<script>
(function () {
  "use strict";

  var DATA = JSON.parse(document.getElementById("payload").textContent);
  var F = {};
  DATA.fields.forEach(function (name, i) { F[name] = i; });

  var ROWS = DATA.rows;
  var PROTEOMES = DATA.proteomes;
  var FAMILIES = DATA.families || [];
  var nRows = ROWS.length;
  var N_OUT = DATA.n_outgroup != null ? DATA.n_outgroup :
    PROTEOMES.filter(function (p) { return p.group === "OUT"; }).length;
  var N_IN = DATA.n_ingroup != null ? DATA.n_ingroup :
    PROTEOMES.filter(function (p) { return p.group === "IN"; }).length;

""" + EL_HELPER_JS + LINKOUT_HELPERS_JS + DOWNLOAD_JS + r"""

  function familyLabel(fam) {
    return fam.rep + " (" + fam.size + " in " + fam.species.length +
      (fam.species.length === 1 ? " species)" : " species)");
  }

  var HAY = new Array(nRows);
  for (var h = 0; h < nRows; h++) {
    var row = ROWS[h];
    var famRep = row[F.fam] >= 0 ? FAMILIES[row[F.fam]].rep : "";
    HAY[h] = (row[F.id] + " " + row[F.gene] + " " + row[F.prod] + " " +
              row[F.pfam_n] + " " + row[F.sprot] + " " + famRep).toLowerCase();
  }

  var state = {
    search: "",
    src: "",
    fsrc: "",
    family: -1,
    noTb: false,
    sort: "priority",
    selected: -1
  };
  var view = [];

  function applyFilters() {
    var q = state.search.trim().toLowerCase();
    var terms = q ? q.split(/\s+/) : [];
    var srcIdx = state.src ? PROTEOMES.findIndex(function (p) { return p.short === state.src; }) : -1;
    var fsrcIdx = state.fsrc ? DATA.fsources.indexOf(state.fsrc) : -2;

    view = [];
    for (var i = 0; i < nRows; i++) {
      var row = ROWS[i];
      if (state.family >= 0 && row[F.fam] !== state.family) continue;
      if (srcIdx >= 0 && row[F.src] !== srcIdx) continue;
      if (state.fsrc && row[F.fsrc] !== fsrcIdx) continue;
      if (state.noTb && row[F.tb_hit]) continue;
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
    else if (s === "frac") cmp = function (a, b) { return (ROWS[b][F.frac] - ROWS[a][F.frac]) || cmpId(a, b); };
    else if (s === "breadth") cmp = function (a, b) { return (ROWS[b][F.out_breadth] - ROWS[a][F.out_breadth]) || cmpId(a, b); };
    else if (s === "pos") cmp = function (a, b) {
      var ca = ROWS[a][F.chrom] || "", cb = ROWS[b][F.chrom] || "";
      if (ca !== cb) return ca < cb ? -1 : 1;
      var sa = ROWS[a][F.start], sb = ROWS[b][F.start];
      if (sa == null && sb == null) return cmpId(a, b);
      if (sa == null) return 1;
      if (sb == null) return -1;
      return (sa - sb) || cmpId(a, b);
    };
    // priority: the strongest loss call is a gene family retained in the fewest ingroup
    // species (cleanest loss), conserved across the most outgroup species (broadest
    // ortholog set), with no ingroup TBLASTN hit — a genomic hit means "might just be a
    // missed gene model", so it ranks below rows with none.
    else cmp = function (a, b) {
      return (ROWS[a][F.in_retained] - ROWS[b][F.in_retained]) ||
             (ROWS[b][F.out_breadth] - ROWS[a][F.out_breadth]) ||
             (ROWS[a][F.tb_hit] - ROWS[b][F.tb_hit]) || cmpId(a, b);
    };
    view.sort(cmp);
  }

  // ---- table --------------------------------------------------------------
  var TBL_COLS = [
    { label: "Protein ID", get: function (r) { return ROWS[r][F.id]; }, cls: "mono", sortKey: "id" },
    { label: "Outgroup source", get: function (r) { return ROWS[r][F.src] >= 0 ? PROTEOMES[ROWS[r][F.src]].short : ""; }, sortKey: "src" },
    { label: "Chrom", get: function (r) { return ROWS[r][F.chrom] || ""; }, cls: "mono", sortKey: "pos" },
    { label: "Start", get: function (r) { return ROWS[r][F.start] != null ? ROWS[r][F.start] : ""; }, cls: "num", sortKey: "pos" },
    { label: "Outgroup breadth", get: function (r) { return ROWS[r][F.out_breadth] + " / " + N_OUT + " species"; }, cls: "num", sortKey: "breadth" },
    { label: "Ingroup retained", get: function (r) { return ROWS[r][F.in_retained] + " / " + N_IN + " species"; }, cls: "num" },
    { label: "Outgroup presence", get: function (r) { return Math.round(ROWS[r][F.frac] * 100) + "%"; }, cls: "num", sortKey: "frac" },
    { label: "Ingroup TBLASTN", get: function (r) { return ROWS[r][F.tb_hit] ? ROWS[r][F.tb_genomes] : "none"; }, cls: "wrap-cell" },
    {
      label: "Gene family", cls: "wrap-cell",
      get: function (r) {
        var fi = ROWS[r][F.fam];
        return fi >= 0 ? FAMILIES[fi].rep + " (" + FAMILIES[fi].size + ")" : "";
      }
    },
    { label: "Gene", get: function (r) { return ROWS[r][F.gene]; } },
    { label: "Product", get: function (r) { return ROWS[r][F.prod]; }, cls: "wrap-cell" },
    { label: "Source of annotation", get: function (r) { return ROWS[r][F.fsrc] >= 0 ? DATA.fsources[ROWS[r][F.fsrc]] : ""; } },
    {
      label: "Pfam domains", cls: "wrap-cell",
      get: function (r) { return ROWS[r][F.pfam_n]; },
      render: function (td, r) { td.appendChild(pfamLinksInline(ROWS[r][F.pfam_n], ROWS[r][F.pfam_a])); }
    }
  ];

  var TBL_PAGE = 300;
  var tblShown = TBL_PAGE;
  var tblBody = document.getElementById("tbl-body");
  var tblScroll = document.getElementById("tbl-scroll");

  function renderTableHead() {
    var thead = document.getElementById("tbl-head");
    thead.textContent = "";
    var tr = document.createElement("tr");
    TBL_COLS.forEach(function (c) {
      var th = el("th", c.cls || null, c.label);
      th.scope = "col";
      if (c.sortKey) {
        th.classList.add("sortable");
        th.tabIndex = 0;
        th.setAttribute("role", "button");
        th.setAttribute("aria-label", "Sort by " + c.label);
        (function (key) {
          th.addEventListener("click", function () {
            state.sort = key;
            document.getElementById("f-sort").value = key;
            refresh(true);
          });
        })(c.sortKey);
        th.addEventListener("keydown", function (e) {
          if (e.key === "Enter" || e.key === " ") { e.preventDefault(); this.click(); }
        });
      }
      tr.appendChild(th);
    });
    thead.appendChild(tr);
    markSortHeader();
  }

  function markSortHeader() {
    var ths = document.querySelectorAll("#tbl-head th");
    TBL_COLS.forEach(function (c, i) {
      if (ths[i]) ths[i].classList.toggle("sort-active", !!c.sortKey && c.sortKey === state.sort);
    });
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
      if (ri === state.selected) tr.className = "sel";
      TBL_COLS.forEach(function (c) {
        var td = el("td", c.cls || null, c.render ? null : (c.get(ri) || ""));
        if (c.render) c.render(td, ri);
        tr.appendChild(td);
      });
      frag.appendChild(tr);
    }
    tblBody.appendChild(frag);
    document.getElementById("tbl-empty").classList.toggle("hidden", view.length > 0);
  }

  function markTableSelection() {
    Array.prototype.forEach.call(tblBody.childNodes, function (tr) {
      tr.className = Number(tr.dataset.ri) === state.selected ? "sel" : "";
    });
  }

  tblScroll.addEventListener("scroll", function () {
    if (tblScroll.scrollTop + tblScroll.clientHeight > tblScroll.scrollHeight - 200 && tblShown < view.length) {
      tblShown = Math.min(tblShown + TBL_PAGE, view.length);
      renderTable(false);
    }
  });
  tblBody.addEventListener("click", function (e) {
    var tr = e.target.closest("tr");
    if (tr && tr.dataset.ri !== undefined) select(Number(tr.dataset.ri));
  });

  // ---- detail panel ---------------------------------------------------------
  var detailEl = document.getElementById("detail");

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
        "Select a protein in the table to see its annotation and database links."));
      return;
    }
    var row = ROWS[ri];
    var sp = row[F.src] >= 0 ? PROTEOMES[row[F.src]] : null;

    detailEl.appendChild(el("h3", null, row[F.id]));
    if (sp) {
      detailEl.appendChild(el("div", "species",
        sp.species + (sp.strain ? " " + sp.strain : "") + " · " + sp.short + " · outgroup"));
    }

    if (row[F.chrom]) {
      detailEl.appendChild(field("Location",
        row[F.chrom] + (row[F.start] != null ? ":" + row[F.start] : "")));
    }

    var retained = row[F.in_retained];
    var statusText = (retained === 0
        ? "Absent from every ingroup proteome"
        : "Retained in " + retained + " of " + N_IN + " ingroup species") +
      "; the gene family is present in " + row[F.out_breadth] + " of " + N_OUT +
      " outgroup species (this protein hits " + Math.round(row[F.frac] * 100) +
      "% of the outgroup).";
    var statusBox = el("div");
    statusBox.appendChild(el("div", "field-value", statusText));
    if (row[F.tb_hit]) {
      var warn = el("span", "badge warn",
        "TBLASTN hit in ingroup genome: " + row[F.tb_genomes]);
      warn.title = "The protein search found no ortholog, but a translated nucleotide " +
        "search did — this may be a missed gene model rather than a true loss.";
      statusBox.appendChild(el("div", null, warn));
    } else if (DATA.tblastn_genomes && DATA.tblastn_genomes.length) {
      statusBox.appendChild(el("div", null, el("span", "badge", "No ingroup TBLASTN hit")));
    }
    detailEl.appendChild(field("Status", statusBox));

    if (row[F.fam] >= 0) {
      var fam = FAMILIES[row[F.fam]];
      var famBox = el("div");
      famBox.appendChild(el("div", "field-value",
        fam.size + " members across " + fam.species.length + " species (" + fam.species.join(", ") + ")"));
      var famLinks = el("div", "links");
      var famBtn = el("button", null, "Show family members (" + fam.size + ")");
      famBtn.type = "button";
      famBtn.addEventListener("click", function () { setFamilyFilter(row[F.fam]); });
      famLinks.appendChild(famBtn);
      famBox.appendChild(famLinks);
      detailEl.appendChild(field("Gene family — independently recovered in multiple outgroup species", famBox));
    }

    if (row[F.gene]) detailEl.appendChild(field("Gene name (outgroup)", row[F.gene]));
    if (row[F.prod]) detailEl.appendChild(field("Product", row[F.prod]));
    if (row[F.fsrc] >= 0) detailEl.appendChild(field("Annotation source", DATA.fsources[row[F.fsrc]]));
    if (row[F.sprot]) detailEl.appendChild(field("Best SwissProt hit", uniprotLinkNode(row[F.sprot])));

    if (row[F.pfam_n]) {
      var pfamCount = row[F.pfam_n].split(",").length;
      detailEl.appendChild(field("Pfam domains (" + pfamCount + ")",
        pfamChipsNode(row[F.pfam_n], row[F.pfam_a])));
    }

    // external links — resolved against the *outgroup* protein, since that is
    // where the gene actually is (e.g. "this looks like S. cerevisiae's ERG3").
    var links = el("div", "links");
    var acc = uniprotAcc(row[F.sprot]);
    var fsrcName = row[F.fsrc] >= 0 ? DATA.fsources[row[F.fsrc]] : "";
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
    if (fsrcName.indexOf("ModelOrg_") === 0) {
      var fd = el("a", null, "FungiDB gene");
      fd.href = "https://fungidb.org/fungidb/app/record/gene/" + encodeURIComponent(geneIdFromProteinId(row[F.id]));
      fd.target = "_blank"; fd.rel = "noopener noreferrer";
      links.appendChild(fd);
    }
    if (!acc && fsrcName.indexOf("ModelOrg_") !== 0) {
      var term = row[F.gene] || geneIdFromProteinId(row[F.id]);
      var nc = el("a", null, "Search NCBI Protein");
      nc.href = "https://www.ncbi.nlm.nih.gov/protein/?term=" + encodeURIComponent(term);
      nc.target = "_blank"; nc.rel = "noopener noreferrer";
      nc.title = "No model-organism or SwissProt hit — worth a manual check";
      links.appendChild(nc);
    }
    if (links.childNodes.length) detailEl.appendChild(field("External resources (outgroup protein)", links));
  }

  function select(ri) {
    state.selected = ri;
    renderDetail();
    markTableSelection();
  }

  // ---- downloads ------------------------------------------------------------
  document.getElementById("dl-tsv").addEventListener("click", function () {
    var lines = [TBL_COLS.map(function (c) { return c.label; }).join("\t")];
    view.forEach(function (ri) {
      lines.push(TBL_COLS.map(function (c) { return String(c.get(ri) || "").replace(/[\t\n\r]/g, " "); }).join("\t"));
    });
    download(DATA.project + ".losses.tsv", lines.join("\n") + "\n", "text/tab-separated-values");
  });

  // ---- summary ----------------------------------------------------------
  function renderSummary() {
    var flagged = 0;
    for (var i = 0; i < nRows; i++) { if (ROWS[i][F.tb_hit]) flagged++; }
    var nOut = PROTEOMES.filter(function (p) { return p.group === "OUT"; }).length;

    var absenceText = (DATA.loss_ingroup_max_frac > 0)
      ? "retained in at most " + Math.round(DATA.loss_ingroup_max_frac * 100) + "% of the ingroup"
      : "absent from every ingroup proteome";
    document.getElementById("summary-note").textContent =
      "Gene families conserved in ≥" + Math.round((DATA.outgroup_min_frac || 0) * 100) +
      "% of the outgroup but " + absenceText + " — candidate lineage-specific gene losses, " +
      "prioritized by ingroup retention (cleanest loss first), outgroup family breadth, and " +
      "absence of ingroup genomic (TBLASTN) evidence.";
    document.getElementById("t-total").textContent = nRows.toLocaleString();
    document.getElementById("t-out").textContent = nOut + (nOut === 1 ? " species" : " species");
    document.getElementById("t-fam").textContent = FAMILIES.length.toLocaleString();
    document.getElementById("t-flagged").textContent =
      nRows ? flagged.toLocaleString() + " (" + Math.round((flagged / nRows) * 100) + "%)" : "—";
  }

  // ---- wiring -------------------------------------------------------------
  function refresh(resetScroll) {
    applyFilters();
    markSortHeader();
    if (state.selected >= 0 && view.indexOf(state.selected) === -1) {
      state.selected = -1;
      renderDetail();
    }
    document.getElementById("count").textContent =
      "Showing " + view.length.toLocaleString() + " of " + nRows.toLocaleString() + " proteins";
    renderTable(true);
  }

  function populateSelects() {
    var src = document.getElementById("f-src");
    src.appendChild(new Option("All outgroup proteomes", ""));
    PROTEOMES.filter(function (p) { return p.group === "OUT"; }).forEach(function (p) {
      src.appendChild(new Option(p.short + " — " + p.species, p.short));
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
  document.getElementById("f-notb").addEventListener("change", function (e) { state.noTb = e.target.checked; refresh(true); });
  document.getElementById("f-sort").addEventListener("change", function (e) { state.sort = e.target.value; refresh(true); });
  document.getElementById("f-reset").addEventListener("click", function () {
    state.search = ""; state.src = ""; state.fsrc = ""; state.family = -1;
    state.noTb = false; state.sort = "priority";
    document.getElementById("f-search").value = "";
    document.getElementById("f-src").value = "";
    document.getElementById("f-fsrc").value = "";
    document.getElementById("f-family").value = "";
    document.getElementById("f-notb").checked = false;
    document.getElementById("f-sort").value = "priority";
    refresh(true);
  });

""" + THEME_TOGGLE_JS + r"""

  // ---- init ---------------------------------------------------------------
  document.getElementById("title").textContent = DATA.project + " — candidate gene losses";
  var inNames = PROTEOMES.filter(function (p) { return p.group === "IN"; })
    .map(function (p) { return p.species; }).join(", ");
  var outNames = PROTEOMES.filter(function (p) { return p.group === "OUT"; })
    .map(function (p) { return p.species; }).join(", ");
  document.getElementById("subtitle").textContent =
    "Ingroup: " + inNames + "  ·  Outgroup: " + outNames;
  document.title = DATA.project + " — NovInvenio candidate losses";

  populateSelects();
  renderSummary();
  renderTableHead();
  refresh(true);
  renderDetail();
})();
</script>
</body>
</html>
"""
