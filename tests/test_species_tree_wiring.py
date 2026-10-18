"""Static checks that the SPECIES_TREE process and its params are wired."""
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CONFIG = (REPO / 'nextflow.config').read_text()
MODULE = (REPO / 'modules' / 'species_tree.nf').read_text()
MAIN = (REPO / 'main.nf').read_text()


def test_params_defined_with_spec_defaults():
    assert re.search(r'^\s*species_tree\s*=\s*true\b', CONFIG, re.M)
    assert re.search(r'^\s*species_tree_markers\s*=\s*10\b', CONFIG, re.M)
    assert re.search(r"^\s*species_tree_markerset\s*=\s*'fungi_odb12'", CONFIG, re.M)
    assert re.search(r'^\s*phyling_markerset_dir\s*=\s*null\b', CONFIG, re.M)


def test_process_runs_the_three_phyling_steps():
    assert 'phyling align' in MODULE
    assert 'phyling filter' in MODULE
    assert '-n ${params.species_tree_markers}' in MODULE
    assert 'phyling tree' in MODULE
    assert '-M ft' in MODULE


def test_process_emits_the_four_contract_files():
    for name in ('species_tree.svg', 'species_tree.rooted.nwk',
                 'species_tree.markers.tsv', 'species_tree.status.txt'):
        assert name in MODULE
    for emit in ('emit: svg', 'emit: nwk', 'emit: markers', 'emit: status'):
        assert emit in MODULE


def test_skip_and_failure_statuses_present():
    assert 'skipped: fewer than 4 species' in MODULE
    assert 'skipped: disabled' in MODULE
    # the tip-label mismatch message comes from bin/species_tree.py (exit 3); the process
    # turns any non-zero exit into `failed: <last stderr line>`.
    assert 'tree tip labels do not match config Short names' in (REPO / 'bin' / 'species_tree.py').read_text()
    assert 'status "failed: \$(tail -1 species_tree.err' in MODULE


def test_main_includes_and_calls_species_tree():
    assert re.search(r"include \{ SPECIES_TREE \} from './modules/species_tree'", MAIN)
    assert re.search(r'^\s*SPECIES_TREE\(', MAIN, re.M)
    assert 'species_tree_svg' in MAIN and 'species_tree_status' in MAIN


REPORT_NF = (REPO / 'workflows' / 'report.nf').read_text()


def test_report_workflow_passes_tree_files_to_collate():
    assert re.search(r'^\s*species_tree_svg\s*//', REPORT_NF, re.M)
    assert re.search(r'^\s*species_tree_status\s*//', REPORT_NF, re.M)
    assert re.search(r'config_csv,\s*species_tree_svg,\s*species_tree_status,?\s*\)', REPORT_NF)
    assert '--species_tree_svg ${species_tree_svg}' in REPORT_NF
    assert '--species_tree_status ${species_tree_status}' in REPORT_NF


def test_main_passes_species_tree_outputs_to_report():
    report_call = MAIN[MAIN.index('REPORT(\n'):]
    assert 'species_tree_svg,' in report_call
    assert 'species_tree_status' in report_call


def test_status_file_carries_marker_and_rooting_lines():
    assert 'markers: ' in MODULE
    assert 'rooting: method=' in MODULE
