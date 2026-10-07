"""Structural checks for issues #193 and #203.

#193: BUILD_ISLANDS (and the marker searches that feed it) must run without
--pangenome_island_pfam_hmm; only the Pfam-dependent steps stay gated.
#203: the three REPORT_TABLES inputs that can be EMPTY_EVALUES_STUB outputs
(always named empty_evalues.tsv) must each stage under a unique name.
"""
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
WF = (REPO / 'workflows' / 'pangenome_profile.nf').read_text()
REPORT = (REPO / 'modules' / 'pangenome' / 'report.nf').read_text()

GATE = 'if (params.pangenome_island_pfam_hmm) {'


def _gated_blocks():
    """Return the body of every top-level `if (params.pangenome_island_pfam_hmm) {` block."""
    blocks = []
    for m in re.finditer(re.escape(GATE), WF):
        depth, i = 1, m.end()
        while depth:
            depth += {'{': 1, '}': -1}.get(WF[i], 0)
            i += 1
        blocks.append(WF[m.end():i])
    return blocks


def test_build_islands_is_not_pfam_gated():
    assert 'BUILD_ISLANDS(' in WF
    assert not any('BUILD_ISLANDS(' in b for b in _gated_blocks())


def test_marker_search_is_not_pfam_gated():
    assert not any('MARKER_HMMSEARCH(' in b for b in _gated_blocks())


def test_pfam_steps_stay_gated():
    gated = '\n'.join(_gated_blocks())
    for step in ('HMMPRESS_PFAM(', 'FAMILY_PFAM_SCAN(', 'DOMAIN_ENRICHMENT(', 'MODULE_DOMAINS(', 'ISLAND_SYNTENY('):
        assert step in gated, step


def test_report_tables_takes_real_islands_unconditionally():
    call = WF.split('REPORT_TABLES(', 1)[1].split(')', 1)[0]
    assert call.split(',')[0].strip() == 'BUILD_ISLANDS.out.islands'
    assert 'EMPTY_SIGNIFICANT_ISLANDS_STUB' not in WF


def test_report_tables_stub_capable_inputs_have_unique_stage_names():
    block = REPORT.split('process REPORT_TABLES', 1)[1].split('output:', 1)[0]
    names = {}
    for arg in ('significant_islands', 'island_pfam_enrichment', 'domtblout', 'strain_inventory'):
        m = re.search(rf"path\({arg}, stageAs: '([^']+)'\)", block)
        assert m, f'{arg} has no stageAs'
        names[arg] = m.group(1)
    assert len(set(names.values())) == len(names)
    assert 'empty_evalues.tsv' not in names.values()


# --- real `nextflow -preview` graph checks ---------------------------------
import shutil  # noqa: E402

import pytest  # noqa: E402

from tests.test_bool_cli_params import minimal_study, _dag  # noqa: E402,F401

NEXTFLOW = shutil.which('nextflow')
needs_nextflow = pytest.mark.skipif(NEXTFLOW is None, reason='nextflow not on PATH')


@needs_nextflow
def test_graph_without_pfam_hmm_still_builds_islands(minimal_study):
    dag = _dag(minimal_study)
    assert 'BUILD_ISLANDS' in dag
    assert 'REPORT_TABLES' in dag
    assert 'FAMILY_PFAM_SCAN' not in dag
    assert 'DOMAIN_ENRICHMENT' not in dag


@needs_nextflow
def test_graph_with_pfam_hmm_keeps_pfam_steps(minimal_study, tmp_path):
    hmm = tmp_path / 'stub.hmm'
    hmm.touch()
    dag = _dag(minimal_study, '--pangenome_island_pfam_hmm', str(hmm))
    for step in ('BUILD_ISLANDS', 'FAMILY_PFAM_SCAN', 'DOMAIN_ENRICHMENT'):
        assert step in dag, step
