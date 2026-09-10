/**
 * Behavioural check for lib/index_page.py's render_alignment_viewer_page()
 * (issue #86): given ?dir=&genome=&protein= in the URL, the page should
 * auto-call window.NIAlignments.open() on load and render the alignment,
 * with no click needed -- this is the standalone "open in new tab" page,
 * not the in-page dialog driven by drive_alignment_popup.mjs.
 *
 * Usage: node drive_alignment_viewer_page.mjs <page.html> <jsdom-path>
 */
import fs from 'node:fs';
import zlib from 'node:zlib';
import { createRequire } from 'node:module';

const [, , PAGE, JSDOM_PATH] = process.argv;
const { JSDOM } = createRequire(import.meta.url)(JSDOM_PATH);

let failures = 0;
function check(name, cond, extra) {
  if (cond) { console.log('PASS ' + name); }
  else { console.log('FAIL ' + name + (extra ? '  <- ' + extra : '')); failures++; }
}

const SHARD_JSON = JSON.stringify({
  memberA: { hits: [{ genome: 'Afum', sseqid: 's1', evalue: 1e-30, bitscore: 100, pident: 80.0,
              length: 5, qstart: 1, qend: 5, sstart: 10, send: 15, sframe: 1,
              qseq: 'MKVLA', sseq: 'MKVLA' }] },
});
const GZIP_BYTES = new Uint8Array(zlib.gzipSync(Buffer.from(SHARD_JSON, 'utf8')));

async function run() {
  const dom = new JSDOM(fs.readFileSync(PAGE, 'utf8'), {
    runScripts: 'dangerously',
    resources: 'usable',
    pretendToBeVisual: true,
    url: 'https://example.test/alignment.html?dir=alignments/&genome=Afum&protein=memberA',
    beforeParse(window) {
      window.Blob = globalThis.Blob;
      window.Response = globalThis.Response;
      window.TextDecoder = globalThis.TextDecoder;
      window.DecompressionStream = globalThis.DecompressionStream;
      window.fetch = async () => ({
        ok: true, status: 200,
        arrayBuffer: async () => GZIP_BYTES.buffer.slice(GZIP_BYTES.byteOffset, GZIP_BYTES.byteOffset + GZIP_BYTES.byteLength),
      });
      window.HTMLDialogElement.prototype.showModal = function () { this.setAttribute('open', ''); };
      window.HTMLDialogElement.prototype.close = function () { this.removeAttribute('open'); };
    },
  });

  await new Promise((r) => setTimeout(r, 30));
  const doc = dom.window.document;
  check('viewer page: auto-opens the dialog from the URL query string',
    doc.getElementById('alignment-dialog').hasAttribute('open'));
  check('viewer page: title reflects the query params, no click needed',
    doc.getElementById('alignment-title').textContent === 'memberA vs s1 (Afum)',
    doc.getElementById('alignment-title').textContent);
  check('viewer page: sets document.title', dom.window.document.title.includes('memberA vs Afum'));

  console.log(failures === 0 ? 'ALL PASSED' : (failures + ' FAILURE(S)'));
  process.exit(failures === 0 ? 0 : 1);
}

run();
