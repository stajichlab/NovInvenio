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
      const ctx = {
        setTransform() {}, clearRect() {}, fillRect() {}, strokeRect() {},
        measureText: (t) => ({ width: (t || '').length * 6 }),
        fillText() {}, save() {}, restore() {}, translate() {}, rotate() {},
        beginPath() {}, moveTo() {}, lineTo() {}, stroke() {},
      };
      window.HTMLCanvasElement.prototype.getContext = () => ctx;
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
  check('n1: NCBI_TaxID gives a direct taxid lookup',
        h.some((x) => x.includes('wwwtax.cgi?id=367110')), h.join(' '));
  check('n1: Pfam chip carries its E-value', /4\.5e-09/.test(detail().textContent));
  check('n1: annotated row gets no remote-homology cluster',
        !h.some((x) => x.includes('hhpred')));
  check('n1: a 1803 aa query is a POST button, not an over-long URL',
        b.some((x) => /BLASTP/.test(x)) && !h.some((x) => x.includes('blast.ncbi.nlm.nih.gov')),
        'buttons=' + b.join(',') + ' hrefs=' + h.filter((x) => x.includes('blast.ncbi')).join(','));

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
  check('core: Pfam E-value reaches the page', /1e-20/.test(det.textContent),
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

console.log(failures === 0 ? 'ALL PASSED' : failures + ' FAILED');
process.exit(failures === 0 ? 0 : 1);
