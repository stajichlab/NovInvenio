import csv
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'lib'))

from config_parser import get_ingroup, get_outgroup, parse_config  # noqa: E402

# Two sibling subphyla under Ascomycota (Pezizomycotina, Taphrinomycotina) plus
# a distant phylum (Basidiomycota), so tests can exercise sibling-default
# outgroup selection, explicit --outgroup-taxon, and stratified capping.
SAMPLES = """\
Species,Strain,Protein,DNA,Short,Lineage
Neurospora crassa,OR74A,Ncra.pep.fa,Ncra.dna.fa,Ncra,Fungi;Ascomycota;Pezizomycotina;Sordariomycetes;Neurospora
Aspergillus fumigatus,Af293,Afum.pep.fa,Afum.dna.fa,Afum,Fungi;Ascomycota;Pezizomycotina;Eurotiomycetes;Aspergillus
Zymoseptoria tritici,,Ztri.pep.fa,Ztri.dna.fa,Ztri,Fungi;Ascomycota;Pezizomycotina;Dothideomycetes;Zymoseptoria
Neolecta irregularis,DAH-3,Nirr.pep.fa,Nirr.dna.fa,Nirr,Fungi;Ascomycota;Taphrinomycotina;Neolectomycetes;Neolecta
Schizosaccharomyces pombe,972h,Spom.pep.fa,Spom.dna.fa,Spom,Fungi;Ascomycota;Taphrinomycotina;Schizosaccharomycetes;Schizosaccharomyces
Saccharomyces cerevisiae,S288c,Scer.pep.fa,Scer.dna.fa,Scer,Fungi;Ascomycota;Taphrinomycotina;Schizosaccharomycetes;Saccharomyces
Cryptococcus neoformans,H99,Cneo.pep.fa,Cneo.dna.fa,Cneo,Fungi;Basidiomycota;Agaricomycotina;Tremellomycetes;Cryptococcus
"""


@pytest.fixture
def samples_csv(tmp_path):
    p = tmp_path / 'samples.csv'
    p.write_text(SAMPLES)
    return p


def run_make_config(samples_csv, output, extra_args=()):
    return subprocess.run(
        [sys.executable, str(REPO / 'bin' / 'make_config.py'),
         '--samples', str(samples_csv),
         '--output', str(output),
         *extra_args],
        capture_output=True, text=True,
    )


def read_config(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh))


def test_taxon_split_with_sibling_default_outgroup(samples_csv, tmp_path):
    output = tmp_path / 'config.csv'
    result = run_make_config(samples_csv, output, ['--ingroup-taxon', 'Pezizomycotina'])
    assert result.returncode == 0, result.stderr

    rows = read_config(output)
    ingroup_shorts = {r['Short'] for r in rows if r['GROUP'] == 'IN'}
    outgroup_shorts = {r['Short'] for r in rows if r['GROUP'] == 'OUT'}

    assert ingroup_shorts == {'Ncra', 'Afum', 'Ztri'}
    # sibling default = other children of Ascomycota (Pezizomycotina's parent)
    assert outgroup_shorts == {'Nirr', 'Spom', 'Scer'}


def test_explicit_outgroup_taxon(samples_csv, tmp_path):
    output = tmp_path / 'config.csv'
    result = run_make_config(samples_csv, output, [
        '--ingroup-taxon', 'Pezizomycotina',
        '--outgroup-taxon', 'Basidiomycota',
    ])
    assert result.returncode == 0, result.stderr

    rows = read_config(output)
    outgroup_shorts = {r['Short'] for r in rows if r['GROUP'] == 'OUT'}
    assert outgroup_shorts == {'Cneo'}


def test_max_per_outgroup_taxon_caps_deterministically(samples_csv, tmp_path):
    output = tmp_path / 'config.csv'
    result = run_make_config(samples_csv, output, [
        '--ingroup-taxon', 'Pezizomycotina',
        '--max-per-outgroup-taxon', '1',
    ])
    assert result.returncode == 0, result.stderr

    rows = read_config(output)
    outgroup_shorts = {r['Short'] for r in rows if r['GROUP'] == 'OUT'}
    # Nirr/Spom/Scer all share Taphrinomycotina as the sibling segment (the
    # lineage token immediately under Ascomycota), so they form one group of
    # 3 that gets capped to 1 — the alphabetically-first Short.
    assert len(outgroup_shorts) == 1
    assert outgroup_shorts == {'Nirr'}

    # deterministic: re-running produces the same result
    output2 = tmp_path / 'config2.csv'
    run_make_config(samples_csv, output2, [
        '--ingroup-taxon', 'Pezizomycotina',
        '--max-per-outgroup-taxon', '1',
    ])
    assert read_config(output2) == rows


def test_ingroup_short_pins_focal_species(samples_csv, tmp_path):
    output = tmp_path / 'config.csv'
    result = run_make_config(samples_csv, output, [
        '--ingroup-short', 'Nirr',
        '--outgroup-taxon', 'Pezizomycotina',
    ])
    assert result.returncode == 0, result.stderr

    rows = read_config(output)
    ingroup_shorts = {r['Short'] for r in rows if r['GROUP'] == 'IN'}
    assert ingroup_shorts == {'Nirr'}


def test_duplicate_short_in_samples_file_errors(tmp_path):
    dup = tmp_path / 'samples.csv'
    dup.write_text(SAMPLES + "Neurospora crassa,OR74A,Ncra.pep.fa,Ncra.dna.fa,Ncra,Fungi;Ascomycota;Pezizomycotina;Sordariomycetes;Neurospora\n")
    output = tmp_path / 'config.csv'
    result = run_make_config(dup, output, ['--ingroup-taxon', 'Pezizomycotina'])
    assert result.returncode != 0
    assert 'Duplicate Short' in result.stderr


def test_output_round_trips_through_config_parser(samples_csv, tmp_path):
    output = tmp_path / 'config.csv'
    result = run_make_config(samples_csv, output, ['--ingroup-taxon', 'Pezizomycotina'])
    assert result.returncode == 0, result.stderr

    parsed = parse_config(output)
    assert {s.short for s in get_ingroup(parsed)} == {'Ncra', 'Afum', 'Ztri'}
    assert {s.short for s in get_outgroup(parsed)} == {'Nirr', 'Spom', 'Scer'}
