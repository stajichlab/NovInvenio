"""
Runtime behaviour of lib/report_common.py's ALIGNMENT_POPUP_JS (issue #74),
driven through jsdom. Not yet wired into any real report template (that's
issue #75) -- this drives the fragment standalone against a minimal fixture
page, the same way tests/test_report_js_behaviour.py drives the full reports.

Skips cleanly unless jsdom resolves -- see that module's docstring for how to
install it (an npm dependency, not a pixi/conda one).
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
DRIVER = Path(__file__).parent / 'js' / 'drive_alignment_popup.mjs'

sys.path.insert(0, str(REPO / 'lib'))
from index_page import render_alignment_viewer_page  # noqa: E402
from report_common import ALIGNMENT_POPUP_HTML, ALIGNMENT_POPUP_JS  # noqa: E402

VIEWER_DRIVER = Path(__file__).parent / 'js' / 'drive_alignment_viewer_page.mjs'

_RESOLVE_JSDOM_JS = (
    'import {createRequire} from "node:module";'
    'const r = createRequire(process.cwd() + "/x.js");'
    'try { console.log(r.resolve("jsdom")); } catch (e) { process.exit(3); }'
)


def _find_jsdom() -> str | None:
    env = os.environ.get('NOVINVENIO_JSDOM')
    if env and Path(env).exists():
        return env
    node = shutil.which('node')
    if not node:
        return None
    proc = subprocess.run(
        [node, '--input-type=module', '-e', _RESOLVE_JSDOM_JS],
        capture_output=True, text=True, cwd=str(REPO), check=False,
    )
    return proc.stdout.strip() or None


def test_alignment_popup_js_syntax_parses():
    node = shutil.which('node')
    if not node:
        pytest.skip('node not available')
    proc = subprocess.run([node, '--check', '-'], input=ALIGNMENT_POPUP_JS,
                          capture_output=True, text=True, check=False)
    assert proc.returncode == 0, proc.stderr


def test_alignment_popup_js_behaviour(tmp_path):
    jsdom = _find_jsdom()
    if not jsdom:
        pytest.skip('jsdom not installed (see tests/test_report_js_behaviour.py docstring)')
    node = shutil.which('node')

    page = tmp_path / 'page.html'
    page.write_text(
        '<!doctype html><html><body>\n' + ALIGNMENT_POPUP_HTML +
        '\n<script>\n' + ALIGNMENT_POPUP_JS + '\n</script>\n</body></html>\n'
    )

    proc = subprocess.run([node, str(DRIVER), str(page), jsdom],
                          capture_output=True, text=True, check=False)
    report = proc.stdout + proc.stderr
    failed = [ln for ln in proc.stdout.splitlines() if ln.startswith('FAIL')]
    assert proc.returncode == 0 and not failed, 'jsdom behaviour checks failed:\n' + report
    assert proc.stdout.count('PASS ') >= 8, report


def test_alignment_viewer_page_syntax_parses():
    node = shutil.which('node')
    if not node:
        pytest.skip('node not available')
    doc = render_alignment_viewer_page()
    import re
    scripts = re.findall(r'<script>\n?(.*?)</script>', doc, re.S)
    assert len(scripts) == 2  # SKIN_BOOT_JS (head) + ALIGNMENT_POPUP_JS/boot (body)
    for script in scripts:
        proc = subprocess.run([node, '--check', '-'], input=script,
                              capture_output=True, text=True, check=False)
        assert proc.returncode == 0, proc.stderr


def test_alignment_viewer_page_behaviour(tmp_path):
    jsdom = _find_jsdom()
    if not jsdom:
        pytest.skip('jsdom not installed (see tests/test_report_js_behaviour.py docstring)')
    node = shutil.which('node')

    page = tmp_path / 'alignment.html'
    page.write_text(render_alignment_viewer_page())

    proc = subprocess.run([node, str(VIEWER_DRIVER), str(page), jsdom],
                          capture_output=True, text=True, check=False)
    report = proc.stdout + proc.stderr
    failed = [ln for ln in proc.stdout.splitlines() if ln.startswith('FAIL')]
    assert proc.returncode == 0 and not failed, 'jsdom behaviour checks failed:\n' + report
    assert proc.stdout.count('PASS ') >= 3, report
