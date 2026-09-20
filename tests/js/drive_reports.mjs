/**
 * Behavioural checks for the report pages' JavaScript, driven through jsdom.
 *
 * tests/test_report_js_behaviour.py generates real novelties.html / core.html
 * from a fixture and runs this against them. The Python suite can only check
 * the payload and the page's static structure -- everything below (the skin
 * picker's three states and persistence, and the whole external-link fallback
 * chain in lib/report_common.py's externalLinksNode) only exists at runtime.
 *
 * Usage:  node drive_reports.mjs <fixture-dir> <jsdom-path>
 *
 * <jsdom-path> may be a package directory or a resolved entry file; it is
 * loaded through createRequire, since an ESM `import` cannot resolve a bare
 * package directory.
 * Prints one PASS/FAIL line per check; exits non-zero if any failed.
 */
import fs from 'node:fs';
import path from 'node:path';
import { createRequire } from 'node:module';

const [, , FX, JSDOM_PATH] = process.argv;
const { JSDOM } = createRequire(import.meta.url)(JSDOM_PATH);

let failures = 0;
function check(name, cond, extra) {
  if (cond) {
    console.log('PASS ' + name);
  } else {
    console.log('FAIL ' + name + (extra ? '  <- ' + extra : ''));
    failures++;
  }
}

function boot(file) {
  return new JSDOM(fs.readFileSync(file, 'utf8'), {
    runScripts: 'dangerously',
    resources: 'usable',
    pretendToBeVisual: true,
    // A real http(s) URL matters: a file:// document gets an opaque origin in
    // jsdom and every localStorage access throws, so the skin-persistence
    // checks below would silently exercise the catch branch instead.
    url: 'https://example.org/reports/' + path.basename(file),
    beforeParse(window) {
      // jsdom implements neither of these, and the page calls both on load.
      window.matchMedia = () => ({ matches: false, addEventListener() {}, removeEventListener() {} });
      // Canvas pixel output isn't observable through jsdom, so drawn text is
      // recorded on the window instead -- island_synteny.html has no table
      // twin of its grid, and this is the only way to prove a redraw
      // actually happened and in what order (used by the row-sort check
      // below; harmless no-op bookkeeping for the other pages).
      window.__fillTextCalls = [];
      const ctx = {
        setTransform() {}, clearRect() {}, fillRect() {}, strokeRect() {},
        measureText: (t) => ({ width: (t || '').length * 6 }),
        fillText(t) { window.__fillTextCalls.push(t); },
        save() {}, restore() {}, translate() {}, rotate() {},
        beginPath() {}, moveTo() {}, lineTo() {}, stroke() {},
      };
      window.HTMLCanvasElement.prototype.getContext = () => ctx;
      // jsdom doesn't implement <dialog>'s modal methods.
      window.HTMLDialogElement.prototype.showModal = function () { this.setAttribute('open', ''); };
      window.HTMLDialogElement.prototype.close = function () { this.removeAttribute('open'); };
    },
  });
}

const ev = (w, t) => new w.Event(t, { bubbles: true });
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const hrefs = (el) => [...el.querySelectorAll('a')].map((a) => a.href);
const btns = (el) => [...el.querySelectorAll('button')].map((b) => b.textContent);

// ---------------------------------------------------------------- novelties
{
  const dom = boot(path.join(FX, 'novelties.html'));
  const w = dom.window, d = w.document;
  const errors = [];
  w.addEventListener('error', (e) => errors.push(String(e.error)));
  await sleep(60);

  check('novelties: loads without error', errors.length === 0, errors.join('; '));

  // ---- skin picker: three states + persistence ----
  const sel = d.getElementById('skin');
  check('novelties: skin picker present', !!sel);
  check('novelties: picker offers follow-system + every skin', sel.options.length >= 3,
        sel.options.length);
  check('novelties: canvas repaint hook wired', typeof w.onSkinChange === 'function');

  sel.value = 'neuromancer';
  sel.dispatchEvent(ev(w, 'change'));
  check('novelties: choosing a skin stamps data-skin',
        d.documentElement.getAttribute('data-skin') === 'neuromancer',
        d.documentElement.getAttribute('data-skin'));
  check('novelties: choosing a skin persists it',
        w.localStorage.getItem('novinvenio.skin') === 'neuromancer',
        w.localStorage.getItem('novinvenio.skin'));

  sel.value = '';
  sel.dispatchEvent(ev(w, 'change'));
  check('novelties: follow-system removes the stamp',
        !d.documentElement.hasAttribute('data-skin'));
  check('novelties: follow-system clears the stored choice',
        w.localStorage.getItem('novinvenio.skin') === null);

  // ---- external-link fallback chain ----
  const tblTab = [...d.querySelectorAll('.tab')].find((t) => /table/i.test(t.textContent));
  tblTab.dispatchEvent(ev(w, 'click'));
  await sleep(30);
  const rows = [...d.querySelectorAll('#tbl-body tr')];
  check('novelties: table renders rows', rows.length > 0, rows.length);

  const pick = (rs, id) => rs.find((r) => r.textContent.includes(id));
  const detail = () => d.getElementById('detail');

  // n1 -- SourceDB "fungidb", a taxid, a Pfam domain, and a 1803 aa sequence.
  pick(rows, 'n1').dispatchEvent(ev(w, 'click'));
  let h = hrefs(detail()), b = btns(detail());
  check('n1: SourceDB drives a FungiDB gene link', h.some((x) => x.includes('fungidb.org')), h.join(' '));
  const fungidbHrefs = h.filter((x) => x.includes('fungidb.org'));
  check('n1: xrefs-derived FungiDB link replaces genomeDbLink\'s, not both',
        fungidbHrefs.length === 1, fungidbHrefs.join(' '));
  check('n1: the surviving FungiDB link uses the xrefs gene ID, not the UniProt/protein_id',
        fungidbHrefs[0] && fungidbHrefs[0].includes('/gene/NCU10683'), fungidbHrefs[0]);
  check('n1: GeneID xref renders an NCBI Gene link',
        h.some((x) => x.includes('ncbi.nlm.nih.gov/gene/5847462')), h.join(' '));
  check('n1: KEGG xref renders a URL-encoded link (colon -> %3A)',
        h.some((x) => x.includes('ncr%3ANCU10683')), h.join(' '));
  check('n1: an unrecognized xref DB renders no link',
        !h.some((x) => x.includes('xyz')), h.join(' '));
  check('n1: NCBI_TaxID gives a direct taxid lookup',
        h.some((x) => x.includes('wwwtax.cgi?id=367110')), h.join(' '));
  // fmtEvalue() (lib/report_common.py) formats via Number(ev).toPrecision(4), which
  // reformats "4.5e-09" to "4.500e-9" (4 sig figs, no leading zero in the exponent) --
  // match the actual rendered form, not the raw fixture string (issue #81).
  check('n1: Pfam chip carries its E-value', /4\.500e-9/.test(detail().textContent));
  check('n1: annotated row gets no remote-homology cluster',
        !h.some((x) => x.includes('hhpred')));
  check('n1: a 1803 aa query is a POST button, not an over-long URL',
        b.some((x) => /BLASTP/.test(x)) && !h.some((x) => x.includes('blast.ncbi.nlm.nih.gov')),
        'buttons=' + b.join(',') + ' hrefs=' + h.filter((x) => x.includes('blast.ncbi')).join(','));

  // ---- presence-chip click popup (e-value + target protein) ----
  {
    const chips = [...detail().querySelectorAll('button.pm.on-pres')];
    const afumChip = chips.find((c) => c.textContent === 'Afum');
    check('n1: Afum presence chip is a clickable button', !!afumChip, chips.map((c) => c.textContent).join(','));
    afumChip.dispatchEvent(ev(w, 'click'));
    const hiDialog = d.getElementById('hitinfo-dialog');
    check('presence chip click: opens the hit-info dialog', hiDialog && hiDialog.hasAttribute('open'));
    const hiText = d.getElementById('hitinfo-body').textContent;
    check('presence chip click: shows Present status', /Present/.test(hiText), hiText);
    check('presence chip click: shows the e-value', /3\.200e-40/.test(hiText), hiText);
    check('presence chip click: shows the resolved target name (gene_name + description)',
          /afuA/.test(hiText) && /Some Aspergillus protein/.test(hiText), hiText);
    d.getElementById('hitinfo-close').dispatchEvent(ev(w, 'click'));
    check('presence chip click: close button closes the dialog', !hiDialog.hasAttribute('open'));
  }

  // n2 -- SourceDB "mycocosm:<portal>", no annotation at all, 120 aa.
  pick(rows, 'n2').dispatchEvent(ev(w, 'click'));
  h = hrefs(detail()); b = btns(detail());
  check('n2: mycocosm:<portal> builds a JGI link',
        h.some((x) => x.includes('mycocosm.jgi.doe.gov')), h.join(' '));
  check('n2: unannotated row offers HHpred', h.some((x) => x.includes('hhpred')), h.join(' '));
  check('n2: unannotated row offers Foldseek', h.some((x) => x.includes('foldseek')));
  check('n2: unannotated row offers InterProScan', h.some((x) => x.includes('interpro')));
  check('n2: Copy FASTA is offered', b.some((x) => /Copy FASTA/.test(x)), b.join(','));
  check('n2: Copy FASTA comes before the tools that need the clipboard',
        detail().textContent.indexOf('Copy FASTA') < detail().textContent.indexOf('HHpred'));
  check('n2: a 120 aa query stays a plain GET link',
        h.some((x) => x.includes('blast.ncbi.nlm.nih.gov')) && !b.some((x) => /BLASTP/.test(x)),
        'buttons=' + b.join(','));

  // shared -- a SwissProt hit, and no sequence (non-novelty under the default
  // --report_sequences novelties).
  const nov = d.getElementById('f-nov');
  nov.checked = false;
  nov.dispatchEvent(ev(w, 'change'));
  await sleep(30);
  const allRows = [...d.querySelectorAll('#tbl-body tr')];
  check('novelties: clearing the novelty filter reveals more rows',
        allRows.length > rows.length, rows.length + ' -> ' + allRows.length);
  pick(allRows, 'shared').dispatchEvent(ev(w, 'click'));
  h = hrefs(detail());
  check('shared: SwissProt accession resolves to UniProt',
        h.some((x) => x.includes('uniprot.org/uniprotkb/P12345')), h.join(' '));
  check('shared: SwissProt accession resolves to AlphaFold',
        h.some((x) => x.includes('alphafold.ebi.ac.uk/entry/P12345')));
  check('shared: a row with no sequence offers no sequence tools',
        !h.some((x) => x.includes('blast.ncbi.nlm.nih.gov')) &&
        !h.some((x) => x.includes('uniprot.org/blast')) &&
        !btns(detail()).some((x) => /Copy FASTA/.test(x)));

  const q = d.getElementById('f-search');
  q.value = 'n2';
  q.dispatchEvent(ev(w, 'input'));
  await sleep(220);
  check('novelties: search box filters', d.getElementById('count').textContent.includes('1'),
        d.getElementById('count').textContent);
}

// --------------------------------------------------------------------- core
{
  const dom = boot(path.join(FX, 'core.html'));
  const w = dom.window, d = w.document;
  await sleep(60);
  check('core: skin picker present', !!d.getElementById('skin'));
  const rows = [...d.querySelectorAll('#tbl-body tr')];
  check('core: table renders rows', rows.length > 0, rows.length);
  rows[0].dispatchEvent(ev(w, 'click'));
  const det = d.getElementById('detail');
  check('core: external links render', hrefs(det).length > 0);
  // fmtEvalue() reformats "1e-20" to "1.000e-20" (Number(ev).toPrecision(4)) --
  // match the actual rendered form (issue #81).
  check('core: Pfam E-value reaches the page', /1\.000e-20/.test(det.textContent),
        det.textContent.slice(0, 200));
  check('core: embeds no sequences, so offers no sequence tools',
        !btns(det).some((x) => /Copy FASTA/.test(x)));
}

// ------------------------------------------------------- hostile SourceDB
// A config CSV travels between users, so a SourceDB value is only
// semi-trusted. The free-form "{gene}" template goes straight into an href,
// and must be rejected unless it is http(s).
{
  const dom = boot(path.join(FX, 'hostile.html'));
  const w = dom.window, d = w.document;
  await sleep(60);
  const tblTab = [...d.querySelectorAll('.tab')].find((t) => /table/i.test(t.textContent));
  tblTab.dispatchEvent(ev(w, 'click'));
  await sleep(30);
  const rows = [...d.querySelectorAll('#tbl-body tr')];
  rows.find((r) => r.textContent.includes('n1')).dispatchEvent(ev(w, 'click'));
  const det = d.getElementById('detail');
  const all = [...det.querySelectorAll('a')].map((a) => a.getAttribute('href') || '');
  check('hostile: a javascript: SourceDB template yields no link',
        !all.some((x) => /^\s*javascript:/i.test(x)), all.join(' | '));
  check('hostile: a data: SourceDB template yields no link',
        !all.some((x) => /^\s*data:/i.test(x)), all.join(' | '));
  // The row must still get its other links -- rejecting the DB link is not a
  // reason to render nothing.
  check('hostile: the rest of the link row still renders',
        all.some((x) => x.includes('ncbi.nlm.nih.gov')), all.join(' | '));
  // And a legitimate https template must still work.
  rows.find((r) => r.textContent.includes('n2')).dispatchEvent(ev(w, 'click'));
  const det2 = d.getElementById('detail');
  check('hostile fixture: an https {gene} template still builds a link',
        hrefs(det2).some((x) => x.includes('custom.example.org/gene/')),
        hrefs(det2).join(' | '));
}

// ---------------------------------------------------------- island synteny
// The most JS-heavy of the four pages (issue #120): a canvas grid, a glyph
// strip, a stateful sortable/filterable sidebar, and a row-sort select that
// re-renders. It has no table twin of the grid, so the row-sort check below
// asserts on the recorded canvas fillText() calls (window.__fillTextCalls,
// wired in boot() above) rather than on any DOM element -- there is no DOM
// element that carries per-row haplotype identity to assert against.
{
  const dom = boot(path.join(FX, 'island_synteny.html'));
  const w = dom.window, d = w.document;
  const errors = [];
  w.addEventListener('error', (e) => errors.push(String(e.error)));
  await sleep(60);

  check('island synteny: loads without error', errors.length === 0, errors.join('; '));

  // ---- sidebar selection redraws the main panel ----
  // renderSidebar() rebuilds .isv-item as fresh nodes on every render (see
  // its `list.textContent = ""` in lib/island_synteny_template.py), so the
  // NodeList must be re-queried after each click -- a node captured before a
  // click is a detached element afterwards and its classList is stale.
  const isvItems = () => [...d.querySelectorAll('.isv-item')];
  check('island synteny: sidebar lists both islands', isvItems().length === 2, isvItems().length);
  // Default sidebar sort is "size" desc: island A (2 families) before B (1).
  check('island synteny: first item selected by default',
        d.getElementById('isv-title').textContent === 'S1:c1:1-2',
        d.getElementById('isv-title').textContent);
  check('island synteny: selected item is marked in the sidebar',
        isvItems()[0].classList.contains('sel') && isvItems()[0].getAttribute('aria-selected') === 'true');

  isvItems()[1].dispatchEvent(ev(w, 'click'));
  check('island synteny: selecting the second island updates the title',
        d.getElementById('isv-title').textContent === 'S2:c2:5-6',
        d.getElementById('isv-title').textContent);
  check('island synteny: selecting the second island moves the sidebar highlight',
        isvItems()[1].classList.contains('sel') && !isvItems()[0].classList.contains('sel'));
  check('island synteny: main note reflects the newly selected island',
        /1 families in locus order/.test(d.getElementById('isv-main-note').textContent),
        d.getElementById('isv-main-note').textContent);

  // Back to island A, which has two haplotypes: S1+S2 (pattern "11", count 2)
  // and S3 (pattern "01", count 1). Default row-sort is "pattern"
  // (lexicographic), so "01" (S3) draws before "11" (S1+S2).
  w.__fillTextCalls.length = 0;
  isvItems()[0].dispatchEvent(ev(w, 'click'));
  const countCalls = () => w.__fillTextCalls.filter((t) => /^×/.test(t));
  check('island synteny: default row sort is pattern order (S3 before S1+S2)',
        countCalls().join(',') === '×1,×2', countCalls().join(','));

  const rowSort = d.getElementById('f-row-sort');
  w.__fillTextCalls.length = 0;
  rowSort.value = 'count';
  rowSort.dispatchEvent(ev(w, 'change'));
  check('island synteny: choosing "count" reorders rows (S1+S2 before S3)',
        countCalls().join(',') === '×2,×1', countCalls().join(','));

  w.__fillTextCalls.length = 0;
  rowSort.value = 'strain';
  rowSort.dispatchEvent(ev(w, 'change'));
  check('island synteny: choosing "strain" reorders rows (S1+S2 sorts before S3)',
        countCalls().join(',') === '×2,×1', countCalls().join(','));

  // ---- skin picker repaints the grid ----
  const sel = d.getElementById('skin');
  check('island synteny: skin picker present', !!sel);
  w.__fillTextCalls.length = 0;
  sel.value = 'neuromancer';
  sel.dispatchEvent(ev(w, 'change'));
  check('island synteny: choosing a skin stamps data-skin',
        d.documentElement.getAttribute('data-skin') === 'neuromancer',
        d.documentElement.getAttribute('data-skin'));
  check('island synteny: choosing a skin persists it',
        w.localStorage.getItem('novinvenio.skin') === 'neuromancer',
        w.localStorage.getItem('novinvenio.skin'));
  check('island synteny: window.onSkinChange is wired', typeof w.onSkinChange === 'function');
  check('island synteny: a skin change triggers a repaint (canvas redrawn)',
        w.__fillTextCalls.length > 0, w.__fillTextCalls.length);
}

// ---------------------------------------------- island synteny: empty state
{
  const dom = boot(path.join(FX, 'island_synteny_empty.html'));
  const w = dom.window, d = w.document;
  const errors = [];
  w.addEventListener('error', (e) => errors.push(String(e.error)));
  await sleep(60);

  check('island synteny empty: loads without error', errors.length === 0, errors.join('; '));
  check('island synteny empty: explorer is hidden',
        d.getElementById('isv-explorer').classList.contains('hidden'));
  check('island synteny empty: empty-state panel is shown',
        !d.getElementById('isv-empty-state').classList.contains('hidden'));

  // 1 located island excluded (single-strain) + 2 qualifying islands
  // truncated by --top_islands 0. lib/island_synteny_template.py's
  // absenceReasons() helper feeds BOTH #summary-note (the always-visible
  // card above #isv-body) and #isv-empty-text (the text inside the
  // empty-state panel itself) from the same array, so they must agree --
  // this is exactly the check that would have caught the original bug,
  // where #isv-empty-text was built from a second, independent, incomplete
  // sentence and silently dropped the truncated count.
  const summaryText = d.getElementById('summary-note').textContent;
  const emptyText = d.getElementById('isv-empty-text').textContent;
  check('island synteny empty: summary-note names the excluded count',
        /1.*excluded/.test(summaryText), summaryText);
  check('island synteny empty: summary-note names the truncated count',
        /2 qualifying island/.test(summaryText), summaryText);
  check('island synteny empty: isv-empty-text names the excluded count',
        /1.*excluded/.test(emptyText), emptyText);
  check('island synteny empty: isv-empty-text names the truncated count',
        /2 qualifying island/.test(emptyText), emptyText);
  const marker = 'none remain to draw: ';
  const reasonsFromEmptyText = emptyText.slice(emptyText.indexOf(marker) + marker.length);
  check('island synteny empty: summary-note and isv-empty-text agree on both reasons',
        emptyText.includes(marker) && summaryText === reasonsFromEmptyText,
        'summary=' + summaryText + ' | empty=' + emptyText);
}

// ---------------------------------- island synteny: truncation-only empty
// excluded=0, truncated=2 -- the case the original bug would have gotten
// backwards (it would have named 0 excluded and said nothing about
// truncation). Both islands qualify; --top_islands 0 truncates both.
{
  const dom = boot(path.join(FX, 'island_synteny_empty_truncated_only.html'));
  const w = dom.window, d = w.document;
  const errors = [];
  w.addEventListener('error', (e) => errors.push(String(e.error)));
  await sleep(60);

  check('island synteny empty (truncation-only): loads without error',
        errors.length === 0, errors.join('; '));
  const emptyText = d.getElementById('isv-empty-text').textContent;
  const summaryText = d.getElementById('summary-note').textContent;
  check('island synteny empty (truncation-only): isv-empty-text names the truncated count',
        /2 qualifying island/.test(emptyText), emptyText);
  check('island synteny empty (truncation-only): isv-empty-text does not claim any exclusion',
        !/excluded/.test(emptyText), emptyText);
  check('island synteny empty (truncation-only): summary-note agrees (no exclusion claimed)',
        !/excluded/.test(summaryText), summaryText);
}

console.log(failures === 0 ? 'ALL PASSED' : failures + ' FAILED');
process.exit(failures === 0 ? 0 : 1);
