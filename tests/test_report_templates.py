"""
Structural guards for the four report page templates.

CLAUDE.md documents ``node --check`` as step 1 of verifying template JS edits,
but nothing ran it. Now that the templates are assembled by concatenating
fragments from lib/report_common.py and lib/skins.py, a mistake in a *fragment*
breaks every page at once -- so the check belongs in the suite rather than in a
procedure someone has to remember. It skips cleanly where node is unavailable.

The behavioural half of that procedure (driving a generated report under jsdom)
still lives in CLAUDE.md: it needs an npm install, which does not belong in a
pytest run.
"""
import re
import shutil
import subprocess
import sys

import pytest
from core_report_template import CORE_HTML_TEMPLATE
from index_page import render_gallery_page, render_project_page
from losses_report_template import LOSSES_HTML_TEMPLATE
from report_template import HTML_TEMPLATE
from skins import SKINS

_INDEX_PAGE = render_project_page(
    project='demo',
    reports=[{'href': 'novelties.html', 'title': 'Novelty candidates', 'desc': 'd'}],
    tiles=[('Ingroup proteomes', '2')],
    ingroup=[{'short': 'Ncra', 'species': 'Neurospora crassa', 'strain': '', 'taxon': 'Pez'}],
    outgroup=[],
)
_GALLERY_PAGE = render_gallery_page(
    projects=[{'title': 'demo', 'href': 'demo/report.html', 'desc': 'd', 'meta': 'm',
               'subpages': [{'href': 'demo/core.html', 'title': 'core'}]}],
)

PAGES = {
    'novelties': HTML_TEMPLATE,
    'core': CORE_HTML_TEMPLATE,
    'losses': LOSSES_HTML_TEMPLATE,
    'report': _INDEX_PAGE,
    'gallery': _GALLERY_PAGE,
}

# The <script> bodies of a page: the payload block is `<script type=...>` and is
# skipped by this pattern, so what is left is executable JS.
_SCRIPT_RE = re.compile(r'<script>\n?(.*?)</script>', re.DOTALL)


def _scripts(page: str) -> list[str]:
    return _SCRIPT_RE.findall(page)


@pytest.mark.parametrize('name', sorted(PAGES))
def test_page_javascript_parses(name, tmp_path):
    node = shutil.which('node')
    if not node:
        pytest.skip('node not available')
    scripts = _scripts(PAGES[name])
    assert scripts, f'{name}: no executable <script> block found'
    for i, body in enumerate(scripts):
        f = tmp_path / f'{name}_{i}.js'
        f.write_text(body)
        proc = subprocess.run([node, '--check', str(f)],
                              capture_output=True, text=True, check=False)
        assert proc.returncode == 0, f'{name} script {i}:\n{proc.stderr}'


@pytest.mark.parametrize('name', sorted(PAGES))
def test_page_is_self_contained(name):
    """No CDN scripts, external stylesheets or webfonts -- the pages open offline."""
    page = PAGES[name]
    assert '<script src=' not in page
    assert '<link rel="stylesheet"' not in page
    assert 'fonts.googleapis.com' not in page
    assert '@import' not in page


@pytest.mark.parametrize('name', sorted(PAGES))
def test_page_offers_every_skin(name):
    page = PAGES[name]
    for skin_id in SKINS:
        assert f'value="{skin_id}"' in page, f'{name} is missing the {skin_id} skin option'
        assert f':root[data-skin="{skin_id}"]' in page


@pytest.mark.parametrize('name', sorted(PAGES))
def test_page_applies_stored_skin_before_first_paint(name):
    """The boot snippet must be in <head>, or the page flashes the default palette."""
    page = PAGES[name]
    head = page.split('</head>', 1)[0]
    assert 'novinvenio.skin' in head, f'{name}: skin boot script is not in <head>'


@pytest.mark.parametrize('name', ['novelties', 'core', 'losses'])
def test_report_has_no_hardcoded_theme_colours(name):
    """Colours belong to a skin. A raw hex in a template is how the old
    light/dark-only pages drifted apart in the first place."""
    page = PAGES[name]
    style = page.split('<style>', 1)[1].split('</style>', 1)[0]
    # Strip the generated skin block -- that is where hex literals legitimately live.
    body_css = style.split('@media print', 1)[-1]
    body_css = body_css.split('}', 3)[-1] if '@media print' in style else style
    stray = re.findall(r'(?<![\w-])#[0-9a-fA-F]{3,8}(?![\w-])', body_css)
    assert not stray, f'{name}: hardcoded colours outside the skin registry: {stray}'


def test_report_template_shares_helpers_rather_than_copying_them():
    """report_template.py used to carry its own copy of the linkout helpers,
    which is how core.html and losses.html ended up resolving proteins
    slightly differently. Guard the consolidation."""
    src = (
        __import__('pathlib').Path(sys.modules['report_template'].__file__).read_text()
    )
    assert 'LINKOUT_HELPERS_JS' in src
    # The old inline definitions must be gone, not merely shadowed.
    assert src.count('function uniprotAcc(') == 0
    assert src.count('function geneIdFromProteinId(') == 0


@pytest.mark.parametrize('name', ['novelties', 'core', 'losses'])
def test_filter_count_is_announced(name):
    """The result count changes on every filter; a screen reader needs to hear it."""
    assert 'id="count" role="status" aria-live="polite"' in PAGES[name]
