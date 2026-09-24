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

It also covers `island_synteny.html` (issue #120) -- the most JS-heavy of the
four report pages: a canvas grid, a glyph strip, a stateful sortable/
filterable sidebar, and a row-sort select that re-renders. That page has no
table twin of its grid, so the row-sort checks assert on
`window.__fillTextCalls` (a call-recording canvas stub wired in
`tests/js/drive_reports.mjs`'s `boot()`) instead of on any DOM element -- see
that file's island-synteny section for why.

Writing this coverage surfaced a real bug (since fixed, same commit series):
the empty-state panel's own text (`#isv-empty-text`) named only the excluded
count, never the count of islands truncated by `--top_islands`, even though
the always-visible `#summary-note` above it always named both via its own
`reasons` array. `lib/island_synteny_template.py` now has one
`absenceReasons()` helper both call sites build their sentence from, so they
cannot drift apart again; the three empty-state fixtures below (both counts
non-zero, truncation-only, and the original three-island one) and the
"summary-note and empty-state agree" check in the driver are what would have
caught the original bug and are what guard against a regression of it now.

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
# Drome: "ncbipep" -- a plain NCBI RefSeq proteome with no FungiDB/MycoCosm
# record, uses the protein ID as-is rather than a gene-ID transform. Scer/Spom
# carry a taxid but no SourceDB, so the NCBI-search fallback is exercised too.
CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup,SourceDB,NCBI_TaxID
IN,Neurospora crassa,OR74A,Ncra.pep.fa,Ncra.dna.fa,Ncra,Pezizomycotina,fungidb,367110
IN,Aspergillus fumigatus,Af293,Afum.pep.fa,Afum.dna.fa,Afum,Pezizomycotina,mycocosm:Aspfu1,330879
IN,Drosophila melanogaster,ISO1,Drome.pep.fa,,Drome,Arthropoda,ncbipep,7227
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
    "protein_id\tsource_proteome\tNcra\tAfum\tDrome\tSpom\tScer\tgene_name\t"
    "product_description\tfunction_source\tBest_Swissprot\tPfam_Names\t"
    "Pfam_Accessions\tPfam_Evalues\tuniprot_xrefs\n"
    # n1/n2 stay present across ALL THREE ingroup species (Ncra/Afum/Drome),
    # not just the original two -- with a 3rd ingroup species now in play,
    # the default ingroup_min_frac recomputation needs their ingroup coverage
    # unchanged (2/2 -> 3/3), or they'd silently drop below the default
    # novelty threshold (2/3 < 0.75) and the whole default-filtered table
    # would come up empty. n3 is deliberately Drome-only (1/3 ingroup
    # coverage) -- a real single-species-specific gene isn't expected to
    # pass the default "conserved across most of the ingroup" novelty test,
    # so it's checked after the novelty filter is cleared, same as `shared`.
    # `shared` also needs Drome=1: core.html's core_min_frac=0.9 needs its
    # presence fraction across ALL FIVE proteomes unchanged at 1.0 (it was
    # 4/4 before Drome existed) -- 4/5=0.8 would silently drop it below
    # 0.9 and empty the whole core report.
    "n1\tNcra\t1\t1\t1\t0\t0\tada-1\tall development altered-1\tModelOrg_Ncra\t\t"
    "bZIP_1\tPF00170.27\t4.5e-09\tVEuPathDB:FungiDB:NCU10683|GeneID:5847462|KEGG:ncr:NCU10683|UnknownDB:xyz\n"
    "n2\tAfum\t1\t1\t1\t0\t0\t\t\t\t\t\t\t\t\n"
    "n3\tDrome\t0\t0\t1\t0\t0\t\t\t\t\t\t\t\t\n"
    "shared\tNcra\t1\t1\t1\t1\t1\t\tconserved thing\tPfam\t"
    "sp|P12345|TEST_YEAST Some protein\tAAA\tPF00004.31\t1e-20\t\n"
)
TBLASTN = "protein_id\tSpom\tScer\nn1\t0\t0\nn2\t0\t1\nn3\t0\t0\n"

# island_synteny.html (issue #120). Two islands so sidebar selection is
# observable (locus_id changes); island A has two haplotypes with different
# counts/patterns so a row-sort change visibly reorders the recorded canvas
# draw calls (S1+S2 share pattern "11" count 2, S3 is "01" count 1 -- pattern
# order and count order disagree, which is what makes the reorder provable).
ISV_ISLANDS = (
    "n_strains\texample_strain\tisland_size\tmember_families\t"
    "n_supporting_pairs\tclassifications\tpfam_domains\tlocus_id\t"
    "locus_contig\tlocus_start\tlocus_end\tn_members_with_coordinates\t"
    "n_contigs_in_locus\n"
    "2\tS1\t2\tfamA,famB\t1\tunexplained_physical\t-\t"
    "S1:c1:1-2\tc1\t1\t2\t2\t1\n"
    "2\tS2\t1\tfamC\t1\tunexplained_physical\t-\t"
    "S2:c2:5-6\tc2\t5\t6\t1\t1\n"
)
ISV_MATRIX = (
    "family\tS1\tS2\tS3\n"
    "famA\tpresent\tpresent\tabsent\n"
    "famB\tpresent\tpresent\tpresent\n"
    "famC\tabsent\tpresent\tpresent\n"
)
ISV_POSITIONS = (
    "Short\tfamily\tcontig\trank\n"
    "S1\tfamA\tc1\t1\n"
    "S1\tfamB\tc1\t2\n"
    "S2\tfamC\tc2\t1\n"
)

# Empty-state fixture: one single-strain island (excluded) plus two
# qualifying islands truncated away by --top_islands 0, so the empty-state
# text must name both the excluded (1) and truncated (2) counts.
ISV_EMPTY_ISLANDS = (
    "n_strains\texample_strain\tisland_size\tmember_families\t"
    "n_supporting_pairs\tclassifications\tpfam_domains\tlocus_id\t"
    "locus_contig\tlocus_start\tlocus_end\tn_members_with_coordinates\t"
    "n_contigs_in_locus\n"
    "1\tS1\t1\tfamA\t1\tunexplained_physical\t-\tS1:c1:1-1\tc1\t1\t1\t1\t1\n"
    "2\tS1\t1\tfamB\t1\tunexplained_physical\t-\tS1:c1:2-2\tc1\t2\t2\t1\t1\n"
    "2\tS1\t1\tfamC\t1\tunexplained_physical\t-\tS1:c1:3-3\tc1\t3\t3\t1\t1\n"
)
ISV_EMPTY_MATRIX = "family\tS1\nfamA\tpresent\nfamB\tpresent\nfamC\tpresent\n"
ISV_EMPTY_POSITIONS = "Short\tfamily\tcontig\trank\n"

# Truncation-ONLY empty-state fixture: both islands qualify (n_strains=2,
# nothing excluded) but --top_islands 0 truncates both -- excluded=0,
# truncated=2. This is the case renderEmptyState() used to get wrong: its
# old wording only ever named the excluded count, so with excluded=0 it
# would have said "...were located, but 0 were excluded ... and none remain
# to draw" -- naming a reason that did not apply and omitting the one that
# did.
ISV_TRUNC_ONLY_ISLANDS = (
    "n_strains\texample_strain\tisland_size\tmember_families\t"
    "n_supporting_pairs\tclassifications\tpfam_domains\tlocus_id\t"
    "locus_contig\tlocus_start\tlocus_end\tn_members_with_coordinates\t"
    "n_contigs_in_locus\n"
    "2\tS1\t1\tfamA\t1\tunexplained_physical\t-\tS1:c1:1-1\tc1\t1\t1\t1\t1\n"
    "2\tS1\t1\tfamB\t1\tunexplained_physical\t-\tS1:c1:2-2\tc1\t2\t2\t1\t1\n"
)
ISV_TRUNC_ONLY_MATRIX = "family\tS1\nfamA\tpresent\nfamB\tpresent\n"
ISV_TRUNC_ONLY_POSITIONS = "Short\tfamily\tcontig\trank\n"


def _fasta() -> str:
    long_seq = 'MKV' + 'ACDEFGHIKLMNPQRSTVWY' * 90      # 1803 aa -> POST branch
    short_seq = 'MKVLLA' * 20                           # 120 aa  -> GET branch
    return (
        f'>n1 unannotated candidate\n{long_seq}\n'
        f'>n2 another\n{short_seq}\n'
        f'>n3 ncbipep candidate\n{short_seq}\n'
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
    # n1's Afum hit -- exercises the presence-chip click popup (e-value + target
    # protein name, resolved through descriptions.tsv).
    (d / 'evalues.tsv').write_text(
        'protein_id\tsource_proteome\tNcra\tAfum\tDrome\tSpom\tScer\n'
        'n1\tNcra\t\t3.2e-40\t\t\t\n'
    )
    (d / 'targets.tsv').write_text(
        'protein_id\tsource_proteome\tNcra\tAfum\tDrome\tSpom\tScer\n'
        'n1\tNcra\t\ttr|Q1|Q1_AFUM\t\t\t\n'
    )
    (d / 'descriptions.tsv').write_text(
        'protein_id\tgene_name\tdescription\n'
        'tr|Q1|Q1_AFUM\tafuA\tSome Aspergillus protein\n'
    )

    # issue #159: n2 rests on two narrow ingroup hits, n1 on none; n3/shared have
    # no sidecar row (not computed) and must not show as zeros.
    (d / 'query_lowcov.tsv').write_text(
        'protein_id\tsource_proteome\tqcov_threshold\tquery_hit_cells\t'
        'query_lowcov_cells\tquery_lowcov_proteomes\n'
        'n1\tNcra\t15.0\t2\t0\t\n'
        'n2\tAfum\t15.0\t2\t2\tDrome,Ncra\n'
    )

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
        *common, '--evalues', str(d / 'evalues.tsv'), '--targets', str(d / 'targets.tsv'),
        '--descriptions', str(d / 'descriptions.tsv'),
        '--query_lowcov', str(d / 'query_lowcov.tsv'), '--output', str(d / 'novelties.html'))
    run('make_report.py', '--config', str(d / 'hostile_config.csv'),
        *common, '--output', str(d / 'hostile.html'))
    run('make_core_report.py', '--matrix', str(d / 'matrix.tsv'),
        '--config', str(d / 'config.csv'), '--core_min_frac', '0.9',
        '--output', str(d / 'core.html'))

    (d / 'isv_islands.tsv').write_text(ISV_ISLANDS)
    (d / 'isv_matrix.tsv').write_text(ISV_MATRIX)
    (d / 'isv_positions.tsv').write_text(ISV_POSITIONS)
    run('pangenome_island_synteny.py',
        '--islands_with_domains', str(d / 'isv_islands.tsv'),
        '--presence_matrix', str(d / 'isv_matrix.tsv'),
        '--family_positions', str(d / 'isv_positions.tsv'),
        '--project', 'demo', '--output', str(d / 'island_synteny.html'))

    (d / 'isv_empty_islands.tsv').write_text(ISV_EMPTY_ISLANDS)
    (d / 'isv_empty_matrix.tsv').write_text(ISV_EMPTY_MATRIX)
    (d / 'isv_empty_positions.tsv').write_text(ISV_EMPTY_POSITIONS)
    run('pangenome_island_synteny.py',
        '--islands_with_domains', str(d / 'isv_empty_islands.tsv'),
        '--presence_matrix', str(d / 'isv_empty_matrix.tsv'),
        '--family_positions', str(d / 'isv_empty_positions.tsv'),
        '--project', 'demo', '--top_islands', '0',
        '--output', str(d / 'island_synteny_empty.html'))

    (d / 'isv_trunc_only_islands.tsv').write_text(ISV_TRUNC_ONLY_ISLANDS)
    (d / 'isv_trunc_only_matrix.tsv').write_text(ISV_TRUNC_ONLY_MATRIX)
    (d / 'isv_trunc_only_positions.tsv').write_text(ISV_TRUNC_ONLY_POSITIONS)
    run('pangenome_island_synteny.py',
        '--islands_with_domains', str(d / 'isv_trunc_only_islands.tsv'),
        '--presence_matrix', str(d / 'isv_trunc_only_matrix.tsv'),
        '--family_positions', str(d / 'isv_trunc_only_positions.tsv'),
        '--project', 'demo', '--top_islands', '0',
        '--output', str(d / 'island_synteny_empty_truncated_only.html'))
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
    assert proc.stdout.count('PASS ') >= 50, report


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
