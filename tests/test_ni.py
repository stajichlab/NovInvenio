"""Unit tests for bin/ni (issue #76): scaffolds a new analysis-deploy repo by
copying a fixed file manifest from a --reference checkout.

Uses a synthetic minimal reference tree (not the real NovInvenio_Investigations
checkout) so this test is hermetic and doesn't depend on that repo existing on
the test machine.
"""
import importlib.util
import subprocess
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
NI = REPO / 'bin' / 'ni'

# bin/ni has no .py extension (so it reads as `ni init ...` on the command
# line) -- load it as a module by explicit path rather than a bare `import`.
# spec_from_file_location can't infer a loader for an extension-less file, so
# an explicit SourceFileLoader is required.
_loader = SourceFileLoader('ni_module', str(NI))
_spec = importlib.util.spec_from_loader('ni_module', _loader)
ni_module = importlib.util.module_from_spec(_spec)
_loader.exec_module(ni_module)


def _make_reference(tmp_path):
    """A minimal reference tree covering every REFERENCE_FILES entry, plus
    an extra file NOT in the manifest (to prove ni doesn't copy everything)."""
    ref = tmp_path / 'reference'
    for rel in ni_module.REFERENCE_FILES:
        p = ref / rel
        if rel == 'assets/logo':
            p.mkdir(parents=True, exist_ok=True)
            (p / 'NI_logo_favicon.ico').write_bytes(b'fake-ico-bytes')
            (p / 'NI_logo_card-96.png').write_bytes(b'fake-png-bytes')
        else:
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(f'# reference content for {rel}\n')
            if rel.endswith('.sh'):
                p.chmod(0o755)
    (ref / 'bin' / 'not_in_manifest.py').parent.mkdir(parents=True, exist_ok=True)
    (ref / 'bin' / 'not_in_manifest.py').write_text('# should not be copied\n')
    return ref


def _run(*args, **env):
    proc = subprocess.run(
        [sys.executable, str(NI), *args],
        capture_output=True, text=True, env=env or None,
    )
    return proc


def test_init_copies_every_manifest_file(tmp_path):
    reference = _make_reference(tmp_path)
    target = tmp_path / 'new_repo'

    proc = _run('init', str(target), '--reference', str(reference))
    assert proc.returncode == 0, proc.stderr

    for rel in ni_module.REFERENCE_FILES:
        assert (target / rel).exists(), f'{rel} was not copied'
    # Not in the manifest -- must not be copied.
    assert not (target / 'bin' / 'not_in_manifest.py').exists()


def test_init_preserves_executable_bit(tmp_path):
    reference = _make_reference(tmp_path)
    target = tmp_path / 'new_repo'
    _run('init', str(target), '--reference', str(reference))

    src_mode = (reference / 'bin' / 'run_study.sh').stat().st_mode
    dst_mode = (target / 'bin' / 'run_study.sh').stat().st_mode
    assert dst_mode == src_mode


def test_init_renders_pixi_toml_with_target_name(tmp_path):
    reference = _make_reference(tmp_path)
    target = tmp_path / 'my_new_study_repo'
    _run('init', str(target), '--reference', str(reference))

    pixi = (target / 'pixi.toml').read_text()
    assert 'name = "my_new_study_repo"' in pixi
    assert '[workspace]' in pixi


def test_init_renders_fresh_domains_yaml_not_references_content(tmp_path):
    reference = _make_reference(tmp_path)
    # Give the reference a domains.yaml with populated, repo-specific content --
    # ni should NOT copy it, since conf/domains.yaml isn't in REFERENCE_FILES.
    (reference / 'conf').mkdir(exist_ok=True)
    (reference / 'conf' / 'domains.yaml').write_text('- slug: fungi\n  name: Fungal\n  desc: Real study data.\n')
    target = tmp_path / 'new_repo'
    _run('init', str(target), '--reference', str(reference))

    domains = (target / 'conf' / 'domains.yaml').read_text()
    assert 'Not yet populated' in domains
    assert 'Real study data' not in domains


def test_init_fails_without_reference(tmp_path):
    target = tmp_path / 'new_repo'
    proc = _run('init', str(target))
    assert proc.returncode != 0
    assert 'reference' in proc.stderr.lower()


def test_init_fails_on_nonexistent_reference(tmp_path):
    target = tmp_path / 'new_repo'
    proc = _run('init', str(target), '--reference', str(tmp_path / 'does_not_exist'))
    assert proc.returncode != 0


def test_init_fails_on_nonempty_target(tmp_path):
    reference = _make_reference(tmp_path)
    target = tmp_path / 'existing'
    target.mkdir()
    (target / 'something').write_text('already here\n')

    proc = _run('init', str(target), '--reference', str(reference))
    assert proc.returncode != 0
    assert 'already exists' in proc.stderr


def test_init_warns_but_succeeds_on_missing_reference_file(tmp_path):
    reference = _make_reference(tmp_path)
    # Delete one manifest file from the reference -- ni should warn, not crash.
    (reference / 'bin' / 'run_study.sh').unlink()
    target = tmp_path / 'new_repo'

    proc = _run('init', str(target), '--reference', str(reference))
    assert proc.returncode == 0, proc.stderr
    assert 'run_study.sh' in proc.stderr
    assert not (target / 'bin' / 'run_study.sh').exists()


def test_init_creates_empty_studies_dir(tmp_path):
    reference = _make_reference(tmp_path)
    target = tmp_path / 'new_repo'
    _run('init', str(target), '--reference', str(reference))

    assert (target / 'studies').is_dir()
