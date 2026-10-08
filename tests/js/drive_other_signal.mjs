/** Behavioural checks for the Outgroup signal controls (issue #208). Usage: node drive_other_signal.mjs <novelties.html> <jsdom-path> */
import fs from 'node:fs';
import { createRequire } from 'node:module';

const [, , FILE, JSDOM_PATH] = process.argv;
const { JSDOM } = createRequire(import.meta.url)(JSDOM_PATH);
let failures = 0;
const check = (n, c, x) => { console.log((c ? 'PASS ' : 'FAIL ') + n + (c || !x ? '' : '  <- ' + x)); if (!c) failures++; };
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

const dom = new JSDOM(fs.readFileSync(FILE, 'utf8'), {
  runScripts: 'dangerously', resources: 'usable', pretendToBeVisual: true,
  url: 'https://example.org/reports/novelties.html',
  beforeParse(w) {
    w.matchMedia = () => ({ matches: false, addEventListener() {}, removeEventListener() {} });
    const ctx = { setTransform() {}, clearRect() {}, fillRect() {}, strokeRect() {}, measureText: (t) => ({ width: (t || '').length * 6 }),
      fillText() {}, save() {}, restore() {}, translate() {}, rotate() {}, beginPath() {}, moveTo() {}, lineTo() {}, stroke() {} };
    w.HTMLCanvasElement.prototype.getContext = () => ctx;
  },
});
const w = dom.window, d = w.document;
const errors = [];
w.addEventListener('error', (e) => errors.push(String(e.error)));
await sleep(80);
check('loads without error', errors.length === 0, errors.join('; '));

const wrap = d.getElementById('f-osig-wrap');
check('signal controls are shown', wrap && !wrap.classList.contains('hidden'));
const tags = () => [...d.querySelectorAll('span.osig')].map((s) => s.textContent);
const ev = (t) => new w.Event(t, { bubbles: true });
check('table shows classes for the novelties', tags().length >= 2, tags().join(','));
const before = tags().slice().sort().join(',');

const sel = d.getElementById('f-osig');
sel.value = 'domain_only'; sel.dispatchEvent(ev('change')); await sleep(30);
check('filter domain_only keeps only domain_only', tags().length > 0 && tags().every((t) => t === 'domain_only'), tags().join(','));

sel.value = ''; sel.dispatchEvent(ev('change')); await sleep(30);
const q = d.getElementById('f-oq');
q.value = '10'; q.dispatchEvent(ev('input')); await sleep(30);
check('lowering the threshold removes domain_only', !tags().includes('domain_only'), tags().join(','));
q.value = '100'; q.dispatchEvent(ev('input')); await sleep(30);
check('raising the threshold removes broad', !tags().includes('broad'), tags().join(','));
check('threshold change alters the visible classes', tags().slice().sort().join(',') !== before, before);

console.log(failures ? 'FAILED ' + failures : 'ALL PASSED');
process.exit(failures ? 1 : 0);
