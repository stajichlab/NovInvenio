"""
Small reusable HTML/CSS/JS string fragments shared by the NovInvenio report
pages (lib/report_template.py, lib/core_report_template.py,
lib/losses_report_template.py, bin/make_index_report.py, view/generate_index.py).

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

from skins import skin_boot_js, skin_picker_html, skin_picker_js, skins_css

# The full skin registry as CSS. Named for what it is now; the old
# THEME_VARS_CSS spelling is gone along with the light/dark-only model.
SKIN_VARS_CSS = skins_css()

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

  // ---- genome-database gene records --------------------------------------
  // Driven by the config CSV's optional SourceDB column, carried through to
  // each payload proteome entry as `source_db` (see lib/config_parser.py).
  // Before this, a gene-record link was gated on the annotation source being
  // a ModelOrg_* entry, so only the handful of model organisms got one even
  // though most of the sample pool lives in FungiDB or MycoCosm.
  //
  // Accepted forms: "fungidb", "mycocosm:<portal>", "ensemblfungi:<species>",
  // "veupathdb:<project>", or any URL template containing "{gene}".
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

  // ---- the whole "External resources" block ------------------------------
  // One builder for all three reports so they can no longer drift apart.
  //   o.id        protein ID
  //   o.gene      gene_name ('' if none)
  //   o.sprot     Best_Swissprot ('' if none)
  //   o.pfam      Pfam_Names ('' if none)
  //   o.fsrcName  annotation source label ('' if none)
  //   o.seq       protein sequence ('' when the payload carries no sequences)
  //   o.proteome  payload proteomes[] entry for the row's source species
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
    if (!db && o.fsrcName && o.fsrcName.indexOf("ModelOrg_") === 0) {
      // Legacy behaviour for configs with no SourceDB column: the model
      // organisms in configs/modelorgs.yaml are all FungiDB-backed.
      db = genomeDbLink("fungidb", o.id);
    }
    if (db) links.appendChild(db);
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
