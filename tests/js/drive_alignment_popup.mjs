/**
 * Behavioural checks for lib/report_common.py's ALIGNMENT_POPUP_JS, driven
 * through jsdom. This fragment isn't wired into any real report template yet
 * (that lands in issue #75), so this drives it standalone: a minimal page
 * carrying just ALIGNMENT_POPUP_HTML + ALIGNMENT_POPUP_JS, with `fetch` and
 * the DecompressionStream/Blob/Response/TextDecoder Web APIs stubbed from
 * Node's own globals (jsdom does not implement most of these itself).
 *
 * Usage: node drive_alignment_popup.mjs <page.html> <jsdom-path>
 * Prints one PASS/FAIL line per check; exits non-zero if any failed.
 */
import fs from 'node:fs';
import zlib from 'node:zlib';
import { createRequire } from 'node:module';

const [, , PAGE, JSDOM_PATH] = process.argv;
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

const HIT_MEMBER = {
  memberA: [{
    genome: 'Afum', sseqid: 'scaffold_1', evalue: 1e-40, bitscore: 150, pident: 75.0,
    length: 10, qstart: 1, qend: 10, sstart: 1000, send: 1120, sframe: 1,
    qseq: 'MKVL-ACDEFGHI', sseq: 'MKVLQACDEFGHX', aligned_as: 'rep1',
  }],
};
const SHARD_JSON = JSON.stringify(HIT_MEMBER);
const GZIP_BYTES = zlib.gzipSync(Buffer.from(SHARD_JSON, 'utf8'));

function boot(fetchImpl, { withDecompression = true } = {}) {
  const dom = new JSDOM(fs.readFileSync(PAGE, 'utf8'), {
    runScripts: 'dangerously',
    resources: 'usable',
    pretendToBeVisual: true,
    url: 'https://example.test/report.html',
    beforeParse(window) {
      window.Blob = globalThis.Blob;
      window.Response = globalThis.Response;
      window.TextDecoder = globalThis.TextDecoder;
      if (withDecompression) window.DecompressionStream = globalThis.DecompressionStream;
      window.fetch = fetchImpl;
      window.HTMLDialogElement.prototype.showModal = function () { this.setAttribute('open', ''); };
      window.HTMLDialogElement.prototype.close = function () { this.removeAttribute('open'); };
    },
  });
  return dom;
}

function fetchOk(bytes) {
  return async () => ({ ok: true, status: 200, arrayBuffer: async () => bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength) });
}

async function run() {
  // --- happy path: gzip shard, gap-aware midline, aligned_as note ---
  {
    const dom = boot(fetchOk(new Uint8Array(GZIP_BYTES)));
    await dom.window.NIAlignments.open('alignments/', 'Afum', 'memberA');
    await new Promise((r) => setTimeout(r, 20));
    const doc = dom.window.document;
    check('happy path: title set', doc.getElementById('alignment-title').textContent === 'memberA vs Afum');
    const stats = doc.getElementById('alignment-stats').textContent;
    check('happy path: stats include evalue/pident/aligned_as', /evalue=1e-40/.test(stats) && /pident=75.0%/.test(stats) && /rep1/.test(stats), stats);
    const body = doc.getElementById('alignment-body').textContent;
    // qseq='MKVL-ACDEFGHI' sseq='MKVLQACDEFGHX' -> midline: match except at
    // the gap column (pos 4, space) and the final X/I mismatch (space).
    check('happy path: body contains query line', body.includes('MKVL-ACDEFGHI'), body);
    check('happy path: body contains subject line', body.includes('MKVLQACDEFGHX'), body);
    // qseq='MKVL-ACDEFGHI' sseq='MKVLQACDEFGHX': gap at position 4 (blank),
    // mismatch I/X at position 12 (blank), '|' everywhere else.
    check('happy path: midline blanks the gap and the mismatch, matches elsewhere', body.includes('|||| ||||||| '), JSON.stringify(body));
    check('happy path: dialog opened', doc.getElementById('alignment-dialog').hasAttribute('open'));
  }

  // --- protein not present in the shard ---
  {
    const dom = boot(fetchOk(new Uint8Array(GZIP_BYTES)));
    await dom.window.NIAlignments.open('alignments/', 'Afum', 'nobodyHome');
    await new Promise((r) => setTimeout(r, 20));
    const body = dom.window.document.getElementById('alignment-body').textContent;
    check('missing protein: shows no-archive message', /No archived alignment/.test(body), body);
  }

  // --- DecompressionStream unsupported: plain fallback message, no throw ---
  {
    const dom = boot(fetchOk(new Uint8Array(GZIP_BYTES)), { withDecompression: false });
    await dom.window.NIAlignments.open('alignments/', 'Afum', 'memberA');
    await new Promise((r) => setTimeout(r, 20));
    const body = dom.window.document.getElementById('alignment-body').textContent;
    check('no DecompressionStream: shows browser-support message', /modern browser/.test(body), body);
  }

  // --- host already transparently decoded gzip (no magic bytes) ---
  {
    const dom = boot(fetchOk(new TextEncoder().encode(SHARD_JSON)));
    await dom.window.NIAlignments.open('alignments/', 'Afum', 'memberA');
    await new Promise((r) => setTimeout(r, 20));
    const title = dom.window.document.getElementById('alignment-title').textContent;
    check('plain (non-gzip) response parses directly', title === 'memberA vs Afum', title);
  }

  // --- close button + backdrop click wiring ---
  {
    const dom = boot(fetchOk(new Uint8Array(GZIP_BYTES)));
    await dom.window.NIAlignments.open('alignments/', 'Afum', 'memberA');
    await new Promise((r) => setTimeout(r, 20));
    const doc = dom.window.document;
    doc.getElementById('alignment-close').dispatchEvent(new dom.window.Event('click', { bubbles: true }));
    check('close button closes the dialog', !doc.getElementById('alignment-dialog').hasAttribute('open'));
  }

  // --- openOrNewTab: plain click opens the dialog, modifier-click opens a new tab ---
  {
    const dom = boot(fetchOk(new Uint8Array(GZIP_BYTES)));
    let openedUrl = null;
    dom.window.open = (url) => { openedUrl = url; return null; };

    await dom.window.NIAlignments.openOrNewTab('alignments/', 'Afum', 'memberA', { ctrlKey: false });
    await new Promise((r) => setTimeout(r, 20));
    check('openOrNewTab: plain click opens the in-page dialog',
      dom.window.document.getElementById('alignment-title').textContent === 'memberA vs Afum');
    check('openOrNewTab: plain click does not call window.open', openedUrl === null);

    dom.window.NIAlignments.openOrNewTab('alignments/', 'Afum', 'memberA', { ctrlKey: true });
    check('openOrNewTab: Ctrl-click opens alignment.html in a new tab',
      openedUrl === 'alignment.html?dir=alignments%2F&genome=Afum&protein=memberA', openedUrl);

    openedUrl = null;
    dom.window.NIAlignments.openOrNewTab('loss_alignments/', 'Ncra', 'memberB', { metaKey: true }, 2);
    check('openOrNewTab: Cmd-click includes a non-zero hit index',
      openedUrl === 'alignment.html?dir=loss_alignments%2F&genome=Ncra&protein=memberB&hit=2', openedUrl);
  }

  console.log(failures === 0 ? 'ALL PASSED' : (failures + ' FAILURE(S)'));
  process.exit(failures === 0 ? 0 : 1);
}

run();
