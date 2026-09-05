"""
Runtime behaviour of the report pages' JavaScript, driven through jsdom.

`tests/test_report_templates.py` covers syntax and static structure; this
covers what only exists once the page runs -- the skin picker's three states
and its persistence, and the whole external-link fallback chain in
`lib/report_common.py`'s `externalLinksNode()`.

That chain is worth the machinery: it has four independent inputs (SwissProt
accession, the config's `SourceDB`, the presence of a Pfam domain, the presence
of a sequence) and its branches are how a protein reaches the right database.
A payload test cannot see any of it, and a typo in a URL template would ship
silently.

**Skipped unless jsdom is importable.** It is a Node dependency and this is a
pixi/conda project, so it is not a hard requirement -- install it wherever you
want this coverage (CI, or locally):

    npm install --prefix ~/.cache/novinvenio-jsdom jsdom

and point `NOVINVENIO_JSDOM` at the resulting `node_modules/jsdom`, or just run
`npm install jsdom` somewhere Node's normal resolution will find it.
"""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
DRIVER = Path(__file__).parent / 'js' / 'drive_reports.mjs'

# Ncra: a plain keyed SourceDB. Afum: the keyed form that takes an argument.
# Scer/Spom carry a taxid but no SourceDB, so the NCBI-search fallback is
# exercised too.
CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup,SourceDB,NCBI_TaxID
IN,Neurospora crassa,OR74A,Ncra.pep.fa,Ncra.dna.fa,Ncra,Pezizomycotina,fungidb,367110
IN,Aspergillus fumigatus,Af293,Afum.pep.fa,Afum.dna.fa,Afum,Pezizomycotina,mycocosm:Aspfu1,330879
OUT,Schizosaccharomyces pombe,972h,Spom.pep.fa,Spom.dna.fa,Spom,Taphrinomycotina,,284812
OUT,Saccharomyces cerevisiae,S288c,Scer.pep.fa,Scer.dna.fa,Scer,Saccharomycotina,,559292
"""

# Ncra's SourceDB is a hostile free-form template; Afum's is a legitimate one.
# A config CSV is copied between users and projects, so this is the realistic
# way an unchecked scheme would arrive.
HOSTILE_CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup,SourceDB,NCBI_TaxID
IN,Neurospora crassa,OR74A,Ncra.pep.fa,Ncra.dna.fa,Ncra,Pezizomycotina,javascript:alert(1)//{gene},367110
IN,Aspergillus fumigatus,Af293,Afum.pep.fa,Afum.dna.fa,Afum,Pezizomycotina,https://custom.example.org/gene/{gene},330879
OUT,Schizosaccharomyces pombe,972h,Spom.pep.fa,Spom.dna.fa,Spom,Taphrinomycotina,,284812
OUT,Saccharomyces cerevisiae,S288c,Scer.pep.fa,Scer.dna.fa,Scer,Saccharomycotina,,559292
"""

# n1: annotated (model-org gene name + Pfam) and long. n2: no annotation at all
# and short. shared: a SwissProt hit, present everywhere so it is not a novelty
# and therefore carries no sequence under the default --report_sequences.
MATRIX = (
    "protein_id\tsource_proteome\tNcra\tAfum\tSpom\tScer\tgene_name\t"
    "product_description\tfunction_source\tBest_Swissprot\tPfam_Names\t"
    "Pfam_Accessions\tPfam_Evalues\n"
    "n1\tNcra\t1\t1\t0\t0\tada-1\tall development altered-1\tModelOrg_Ncra\t\t"
    "bZIP_1\tPF00170.27\t4.5e-09\n"
    "n2\tAfum\t1\t1\t0\t0\t\t\t\t\t\t\t\n"
    "shared\tNcra\t1\t1\t1\t1\t\tconserved thing\tPfam\t"
    "sp|P12345|TEST_YEAST Some protein\tAAA\tPF00004.31\t1e-20\n"
)
TBLASTN = "protein_id\tSpom\tScer\nn1\t0\t0\nn2\t0\t1\n"


def _fasta() -> str:
    long_seq = 'MKV' + 'ACDEFGHIKLMNPQRSTVWY' * 90      # 1803 aa -> POST branch
    short_seq = 'MKVLLA' * 20                           # 120 aa  -> GET branch
    return (
        f'>n1 unannotated candidate\n{long_seq}\n'
        f'>n2 another\n{short_seq}\n'
        f'>shared thing\n{"MKQTA" * 30}\n'
    )


# Ask Node to resolve jsdom itself, so any normal resolution path (a global
# install, an npx cache, a node_modules above the repo) works without us
# guessing at directory layouts.
_RESOLVE_JSDOM_JS = (
    'import {createRequire} from "node:module";'
    'const r = createRequire(process.cwd() + "/x.js");'
    'try { console.log(r.resolve("jsdom")); } catch (e) { process.exit(3); }'
)


def _find_jsdom() -> str | None:
    """Locate a jsdom module path, or None to skip."""
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


@pytest.fixture(scope='module')
def fixture_dir(tmp_path_factory):
    """Generate real report pages the driver can open."""
    d = tmp_path_factory.mktemp('reports')
    (d / 'config.csv').write_text(CONFIG)
    (d / 'hostile_config.csv').write_text(HOSTILE_CONFIG)
    (d / 'matrix.tsv').write_text(MATRIX)
    (d / 'tblastn.tsv').write_text(TBLASTN)
    (d / 'candidates.fa').write_text(_fasta())

    def run(script, *args):
        proc = subprocess.run(
            [sys.executable, str(REPO / 'bin' / script), *args],
            capture_output=True, text=True, check=False,
        )
        assert proc.returncode == 0, f'{script} failed:\n{proc.stderr}'

    common = ['--matrix', str(d / 'matrix.tsv'),
              '--tblastn_summary', str(d / 'tblastn.tsv'),
              '--candidates_fa', str(d / 'candidates.fa')]
    run('make_report.py', '--config', str(d / 'config.csv'),
        *common, '--output', str(d / 'novelties.html'))
    run('make_report.py', '--config', str(d / 'hostile_config.csv'),
        *common, '--output', str(d / 'hostile.html'))
    run('make_core_report.py', '--matrix', str(d / 'matrix.tsv'),
        '--config', str(d / 'config.csv'), '--core_min_frac', '0.9',
        '--output', str(d / 'core.html'))
    return d


def test_report_javascript_behaviour(fixture_dir):
    jsdom = _find_jsdom()
    if not jsdom:
        pytest.skip('jsdom not installed (see this module\'s docstring)')
    proc = subprocess.run(
        [shutil.which('node'), str(DRIVER), str(fixture_dir), jsdom],
        capture_output=True, text=True, check=False,
    )
    report = proc.stdout + proc.stderr
    failed = [ln for ln in proc.stdout.splitlines() if ln.startswith('FAIL')]
    assert proc.returncode == 0 and not failed, (
        'jsdom behaviour checks failed:\n' + report
    )
    # Guard against the driver silently doing nothing.
    assert proc.stdout.count('PASS ') >= 25, report


def test_hostile_source_db_is_rejected_at_the_python_layer(fixture_dir):
    """The scheme check lives in JS, but the value must survive the payload
    intact -- so if the JS guard is ever removed, the jsdom test above is what
    catches it. This just pins that the payload really does carry the hostile
    value, i.e. that the JS test is exercising what it claims to."""
    page = (fixture_dir / 'hostile.html').read_text()
    payload = json.loads(
        page.split('<script type="application/json" id="payload">', 1)[1]
            .split('</script>', 1)[0]
    )
    ncra = next(p for p in payload['proteomes'] if p['short'] == 'Ncra')
    assert ncra['source_db'].startswith('javascript:')
