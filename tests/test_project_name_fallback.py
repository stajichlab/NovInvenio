"""Issue #194: Helpers.projectName(params) fallback order, and pangenome.nf's
--pangenome_project/--project warning.

Root cause (confirmed with a minimal local reproduction, see this test file's
`test_reassigning_a_config_declared_param_silently_no_ops` below and
lib/Helpers.groovy's projectName() doc comment): pangenome.nf used to try to
supply the samplesheet-basename default itself, at the top of its workflow{}
block, by assigning `params.project = ...` when neither --project nor
--pangenome_project was given. On Nextflow 26.04.6 that assignment silently
does not persist for a param that is declared (even with a null default) in
nextflow.config's `params {}` block -- `project` and `pangenome_project` both
are. So params.project stayed null, the warning's own `${params.project}`
interpolation printed the literal text "null", and every publishDir/storeDir
closure that calls Helpers.projectName(params) fell through to the 'output'
literal regardless -- verified on a real run (fusarium_FOL ni_smoke_norescue,
0dddb42): outputs landed in results/ni_smoke_norescue/output/pangenome/
instead of results/ni_smoke_norescue/config/pangenome/.

The fix moves the samplesheet-basename fallback into projectName() itself, so
it no longer depends on a script-side mutation of params.project. These
tests exercise the real Groovy code (no Python re-implementation to drift
out of sync) via `nextflow run` against a tiny throwaway probe script that
loads this repo's lib/Helpers.groovy the same way pangenome.nf does (Nextflow
auto-loads lib/*.groovy relative to the main script's directory).
"""
import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
NEXTFLOW = shutil.which('nextflow')

pytestmark = pytest.mark.skipif(NEXTFLOW is None, reason='nextflow not on PATH')

PROBE_SCRIPT = '''#!/usr/bin/env nextflow
nextflow.enable.dsl=2
workflow {
    println "PROJECTNAME=" + Helpers.projectName(params)
}
'''


@pytest.fixture
def probe(tmp_path):
    """A throwaway Nextflow project: probe.nf + lib/Helpers.groovy copied from REPO,
    so Nextflow's automatic lib/ loading picks up the real, unmodified class."""
    (tmp_path / 'lib').mkdir()
    shutil.copy(REPO / 'lib' / 'Helpers.groovy', tmp_path / 'lib' / 'Helpers.groovy')
    script = tmp_path / 'probe.nf'
    script.write_text(PROBE_SCRIPT)
    return script


def _project_name(probe, *extra_args):
    proc = subprocess.run(
        [NEXTFLOW, '-q', 'run', str(probe), *extra_args],
        cwd=probe.parent, capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    lines = [l for l in proc.stdout.splitlines() if l.startswith('PROJECTNAME=')]
    assert lines, proc.stdout + proc.stderr
    return lines[0].split('=', 1)[1]


def test_falls_back_to_output_literal_when_nothing_set(probe):
    assert _project_name(probe) == 'output'


def test_uses_params_project_when_set(probe):
    assert _project_name(probe, '--project', 'foo') == 'foo'


def test_uses_pangenome_project_when_project_unset(probe):
    assert _project_name(probe, '--pangenome_project', 'bar') == 'bar'


def test_uses_config_basename_extension_stripped(probe):
    assert _project_name(probe, '--config', '/tmp/whatever/myconf.csv') == 'myconf'


def test_uses_pangenome_samplesheet_basename_extension_stripped(probe):
    assert _project_name(probe, '--pangenome_samplesheet', '/tmp/whatever/sheet.csv') == 'sheet'


def test_project_beats_pangenome_project(probe):
    assert _project_name(probe, '--project', 'foo', '--pangenome_project', 'bar') == 'foo'


def test_pangenome_project_beats_config(probe):
    assert _project_name(
        probe, '--pangenome_project', 'bar', '--config', '/tmp/whatever/myconf.csv') == 'bar'


def test_config_beats_pangenome_samplesheet(probe):
    assert _project_name(
        probe, '--config', '/tmp/whatever/myconf.csv',
        '--pangenome_samplesheet', '/tmp/whatever/sheet.csv') == 'myconf'


def test_reassigning_a_config_declared_param_silently_no_ops(tmp_path):
    """Documents the underlying Nextflow behaviour that made pangenome.nf's old
    self-assignment approach unreliable: a param declared in nextflow.config's
    `params {}` block (even with a null default) cannot be usefully reassigned
    from inside the script/workflow body -- a read immediately afterward still
    sees the original value. This is not a projectName() test; it is the
    reproduction that explains why projectName() must do the fallback itself."""
    (tmp_path / 'nextflow.config').write_text('params {\n    project = null\n}\n')
    (tmp_path / 'repro.nf').write_text(
        '#!/usr/bin/env nextflow\n'
        'nextflow.enable.dsl=2\n'
        'workflow {\n'
        '    params.project = "world"\n'
        '    println "REASSIGNED=" + params.project\n'
        '}\n'
    )
    proc = subprocess.run(
        [NEXTFLOW, '-q', 'run', str(tmp_path / 'repro.nf')],
        cwd=tmp_path, capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    lines = [l for l in proc.stdout.splitlines() if l.startswith('REASSIGNED=')]
    assert lines, proc.stdout + proc.stderr
    # The whole point: this is "null" (the assignment did not take), not "world".
    assert lines[0] == 'REASSIGNED=null'


def test_pangenome_nf_warning_no_longer_prints_null_literal():
    """Regression for the exact symptom reported in #194: the warning text used
    to interpolate a params.project that had just been (silently, ineffectively)
    assigned, so it always printed the literal word 'null'. It must now report
    the value Helpers.projectName(params) actually resolves to."""
    text = (REPO / 'pangenome.nf').read_text()
    block = text.split('workflow {', 1)[1].split('samples_ch = Channel')[0]
    # Strip line comments before checking for a live assignment, so this
    # doesn't false-positive on the doc comment explaining WHY the old
    # `params.project = ...` self-assignment was removed.
    code_only = '\n'.join(
        line for line in block.splitlines() if not line.strip().startswith('//'))
    assert not re.search(r"params\.project\s*=\s*\S", code_only), (
        "pangenome.nf should no longer self-assign params.project -- "
        "Helpers.projectName(params) covers the fallback")
    assert 'Helpers.projectName(params)' in block
