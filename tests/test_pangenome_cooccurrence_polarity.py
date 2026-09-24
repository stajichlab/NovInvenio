"""Outgroup-frequency polarity for bin/pangenome_cooccurrence.py (design decision Q4,
2026-09-20: "KEEP direction, ADD asymmetry_a, GATE ... with FIXED thresholds").

`direction_a` requires a family to be present in EVERY outgroup strain to call "loss".
With a many-strain outgroup (the Coccidioides reciprocal studies: 169 or 360 outgroup
strains) one missed call turns it "ambiguous", so loss is close to impossible. These
tests pin the additive replacement: `asymmetry_a` (outgroup carrier fraction) and
`direction_a_freq` (thresholded on it), reported ALONGSIDE an unchanged `direction_a`.
"""
import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'lib'))
from pangenome_matrix import PresenceMatrix  # noqa: E402

spec = importlib.util.spec_from_file_location(
    'pangenome_cooccurrence', REPO / 'bin' / 'pangenome_cooccurrence.py')
coocc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(coocc)

spec_pc = importlib.util.spec_from_file_location(
    'pangenome_pair_classification', REPO / 'bin' / 'pangenome_pair_classification.py')
pairclass = importlib.util.module_from_spec(spec_pc)
spec_pc.loader.exec_module(pairclass)


# --- polarize_by_frequency ---------------------------------------------------------

@pytest.mark.parametrize('present,total,expected', [
    (10, 10, 'loss'),
    (9, 10, 'loss'),        # 0.9 >= loss_min 0.9
    (8, 10, 'ambiguous'),
    (2, 10, 'ambiguous'),
    (1, 10, 'gain'),        # 0.1 <= gain_max 0.1
    (0, 10, 'gain'),
    (0, 0, 'ambiguous'),    # no outgroup: never a confident call
])
def test_polarize_by_frequency(present, total, expected):
    assert coocc.polarize_by_frequency(present, total, 0.9, 0.1) == expected


def test_single_outgroup_matches_the_strict_rule():
    # With one outgroup strain the fraction is 0 or 1, so both rules agree.
    for present in (0, 1):
        assert (coocc.polarize_by_frequency(present, 1, 0.9, 0.1)
                == coocc.polarize_direction(present, 1))


def test_asymmetry_fraction():
    assert coocc.outgroup_fraction(9, 10) == pytest.approx(0.9)
    assert coocc.outgroup_fraction(0, 0) is None


@pytest.mark.parametrize('loss_min,gain_max', [(0.1, 0.9), (0.5, 0.5), (1.1, 0.1), (0.9, -0.1)])
def test_invalid_thresholds_rejected(loss_min, gain_max):
    with pytest.raises(ValueError):
        coocc.validate_polarity_thresholds(loss_min, gain_max)


# --- find_cooccurring_pairs wiring -------------------------------------------------

INGROUP = ['A', 'B', 'C', 'D']
OUTGROUP = [f'O{i}' for i in range(10)]


def run(present_by_family, **kw):
    m = PresenceMatrix(families=list(present_by_family), strains=INGROUP + OUTGROUP)
    for fam, carriers in present_by_family.items():
        for s in carriers:
            m.set_call(fam, s, 'present')
    freq = [{'family': f, 'bin': 'shell'} for f in present_by_family]
    outgroup_presence = {
        f: (sum(1 for s in OUTGROUP if m.is_present(f, s)), len(OUTGROUP))
        for f in present_by_family
    }
    return coocc.find_cooccurring_pairs(
        m, freq, {s: 'c1' for s in INGROUP + OUTGROUP}, outgroup_presence,
        min_strain_count=2, fdr_alpha=1.0, strains=INGROUP, screen_alpha=1.0, **kw)


# X is carried by 9 of 10 outgroup strains: ancestral, lost in C,D. The strict rule
# calls it ambiguous (one outgroup miss); the frequency rule calls it loss.
NINE_OF_TEN = {'X': {'A', 'B'} | set(OUTGROUP[:9]), 'Y': {'A', 'B'}}


def test_rows_carry_fraction_and_frequency_call_alongside_strict_call():
    rows = run(NINE_OF_TEN)
    row = next(r for r in rows if r['family_a'] == 'X')
    assert row['direction_a'] == 'ambiguous'
    assert row['asymmetry_a'] == pytest.approx(0.9)
    assert row['direction_a_freq'] == 'loss'


def test_thresholds_are_passed_through():
    rows = run(NINE_OF_TEN, polarity_loss_min=0.95, polarity_gain_max=0.05)
    row = next(r for r in rows if r['family_a'] == 'X')
    assert row['direction_a_freq'] == 'ambiguous'


# --- main(): header and values -----------------------------------------------------

def write_inputs(tmp_path):
    strains = INGROUP + OUTGROUP
    cfg = tmp_path / 'config.csv'
    cfg.write_text('GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup\n' + ''.join(
        f"{'IN' if s in INGROUP else 'OUT'},sp,{s},{s}.pep,{s}.dna,{s}.gff3,{s},c1\n"
        for s in strains))
    mat = tmp_path / 'matrix.tsv'
    lines = ['family\t' + '\t'.join(strains)]
    for fam, carriers in NINE_OF_TEN.items():
        lines.append(fam + '\t' + '\t'.join('present' if s in carriers else 'absent' for s in strains))
    mat.write_text('\n'.join(lines) + '\n')
    freq = tmp_path / 'freq.tsv'
    freq.write_text('family\tfrequency\tstrain_count\tbin\nX\t0.5\t2\tshell\nY\t0.5\t2\tshell\n')
    return cfg, mat, freq


def test_main_writes_new_columns(tmp_path, monkeypatch):
    cfg, mat, freq = write_inputs(tmp_path)
    out = tmp_path / 'pairs.tsv'
    monkeypatch.setattr(sys, 'argv', [
        'pangenome_cooccurrence.py', '--matrix', str(mat), '--frequency_table', str(freq),
        '--config', str(cfg), '--min_strain_count', '2', '--fdr_alpha', '1.0',
        '--screen_alpha', '1.0', '--output', str(out)])
    coocc.main()
    header, *rows = out.read_text().splitlines()
    cols = header.split('\t')
    assert cols[:8] == ['family_a', 'family_b', 'jaccard', 'fisher_p', 'fdr_q',
                        'permutation_p', 'direction_a', 'clade_composition']
    assert cols[8:] == ['asymmetry_a', 'direction_a_freq']
    x = dict(zip(cols, next(r.split('\t') for r in rows if r.startswith('X\t'))))
    assert x['direction_a'] == 'ambiguous'
    assert x['asymmetry_a'] == '0.9000'
    assert x['direction_a_freq'] == 'loss'


def test_main_rejects_crossed_thresholds(tmp_path, monkeypatch):
    cfg, mat, freq = write_inputs(tmp_path)
    monkeypatch.setattr(sys, 'argv', [
        'pangenome_cooccurrence.py', '--matrix', str(mat), '--frequency_table', str(freq),
        '--config', str(cfg), '--polarity_loss_min_frac', '0.2',
        '--polarity_gain_max_frac', '0.8', '--output', str(tmp_path / 'o.tsv')])
    with pytest.raises(SystemExit):
        coocc.main()


# --- pair_classification passes the new columns through ---------------------------

def test_pair_classification_passthrough_columns():
    idx = {'family_a': 0, 'asymmetry_a': 8, 'direction_a_freq': 9}
    assert pairclass.passthrough_columns(idx) == ['asymmetry_a', 'direction_a_freq']
    assert pairclass.passthrough_columns({'family_a': 0}) == []
