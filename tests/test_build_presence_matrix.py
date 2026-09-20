import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / 'bin' / 'build_presence_matrix.py'

CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
IN,In one,,in1.pep.fa,,In1,X
IN,In two,,in2.pep.fa,,In2,X
OUT,Out one,,out1.pep.fa,,Out1,Y
"""

HIT_HEADER = 'query_id\ttarget_id\tevalue\tbitscore\tquery_proteome\ttarget_proteome\n'


@pytest.fixture
def run_dir(tmp_path):
    (tmp_path / 'config.csv').write_text(CONFIG)
    return tmp_path


PARALOG_HEADER = 'protein_ID\tparalog_protein_ID\tbitscore\tevalue\n'


def run(run_dir, hits_text, query_group='IN', min_frac='1.0', other_max_frac='0.0',
        paralog_text=None, competition_scope=None, rescue_evalue=None,
        rescue_delta=None):
    hits_path = run_dir / 'hits.tsv'
    hits_path.write_text(HIT_HEADER + hits_text)
    matrix_out = run_dir / 'matrix.tsv'
    candidates_out = run_dir / 'candidates.txt'
    cmd = [sys.executable, str(SCRIPT),
           '--hits', str(hits_path),
           '--config', str(run_dir / 'config.csv'),
           '--ingroup-min-frac', min_frac,
           '--query-group', query_group,
           '--other-max-frac', other_max_frac,
           '--output-matrix', str(matrix_out),
           '--output-candidates', str(candidates_out)]
    if paralog_text is not None:
        paralog_path = run_dir / 'paralog_cutoffs.tsv'
        paralog_path.write_text(PARALOG_HEADER + paralog_text)
        cmd += ['--paralog-cutoffs', str(paralog_path)]
    if competition_scope is not None:
        cmd += ['--paralog-competition-scope', competition_scope]
    if rescue_evalue is not None:
        cmd += ['--paralog-rescue-evalue', rescue_evalue]
    if rescue_delta is not None:
        cmd += ['--paralog-rescue-delta', rescue_delta]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    matrix = pd.read_csv(matrix_out, sep='\t')
    candidates = candidates_out.read_text().splitlines() if candidates_out.stat().st_size else []
    return matrix, candidates


def test_default_query_group_is_ingroup_novelty_direction(run_dir):
    # g1 (In1) hits In2 but never Out1 — a novelty candidate under the default
    # (unchanged) ingroup-as-query behaviour.
    matrix, candidates = run(run_dir, 'g1\tg1_in2\t1e-10\t100\tIn1\tIn2\n')
    assert candidates == ['In1::g1']
    row = matrix[matrix['protein_id'] == 'g1'].iloc[0]
    assert row['In1'] == 1 and row['In2'] == 1 and row['Out1'] == 0


def test_query_group_out_finds_a_loss_candidate(run_dir):
    # h1 (Out1) hits nothing in the ingroup — a candidate loss: conserved in
    # the (single-member) outgroup query set, absent from the ingroup.
    matrix, candidates = run(run_dir, 'h1\th1_self\t1e-10\t100\tOut1\tOut1\n',
                              query_group='OUT')
    assert candidates == ['Out1::h1']
    row = matrix[matrix['protein_id'] == 'h1'].iloc[0]
    assert row['Out1'] == 1 and row['In1'] == 0 and row['In2'] == 0


def test_query_group_out_excludes_hits_present_in_ingroup(run_dir):
    # h2 (Out1) also hits In1 — present in 1/2 of the ingroup, so not a loss
    # under the default strict absence (--other-max-frac 0.0).
    matrix, candidates = run(
        run_dir,
        'h2\th2_in1\t1e-10\t100\tOut1\tIn1\n',
        query_group='OUT',
    )
    assert candidates == []
    row = matrix[matrix['protein_id'] == 'h2'].iloc[0]
    assert row['Out1'] == 1 and row['In1'] == 1


def test_query_group_out_allows_nearly_missing_with_other_max_frac(run_dir):
    # Same h2 present in 1/2 of the ingroup (frac 0.5): excluded at the default
    # 0.0 above, but kept once --other-max-frac allows up to half the ingroup —
    # a "nearly missing" (lost from most, not all, of the ingroup) candidate.
    matrix, candidates = run(
        run_dir,
        'h2\th2_in1\t1e-10\t100\tOut1\tIn1\n',
        query_group='OUT',
        other_max_frac='0.5',
    )
    assert candidates == ['Out1::h2']


def test_flat_default_evalue_rejects_a_weak_hit(run_dir):
    # Filter 1 is a flat significance floor (--default-evalue, 1e-5 by default), not a
    # per-query paralog-derived one.
    matrix, candidates = run(run_dir, 'g1\tg1_in2\t1e-3\t20\tIn1\tIn2\n',
                             paralog_text='g1\tg1p\t42\t1e-8\n')
    assert (matrix['protein_id'] == 'g1').sum() == 0
    assert candidates == []


def test_supplied_paralog_cutoff_no_longer_rejects_a_strong_hit(run_dir):
    # Regression test for the 2026-09-03 fix: a strong hit (1e-10, well within the flat
    # default 1e-5) used to be rejected outright whenever the query's own in-genome
    # paralog e-value (here 1e-20) was tighter than the hit -- an absolute-magnitude
    # proxy, not an actual test of whether the paralog explains this hit. Filter 2 (the
    # real paralog-competition test) doesn't fire here since g1p was never itself
    # searched against In2.
    matrix, candidates = run(run_dir, 'g1\tg1_in2\t1e-10\t100\tIn1\tIn2\n',
                             paralog_text='g1\tg1p\t42\t1e-20\n')
    row = matrix[matrix['protein_id'] == 'g1'].iloc[0]
    assert row['In2'] == 1
    assert 'In1::g1' in candidates


# --- Filter 2 (paralog competition) scope -----------------------------------
#
# Mirrors the pezizo5_fungi hexA/hex-1 case (see docs/hexA_filtering.md): hexA
# (In1) is a true ortholog of hex1 (In2's HEX-1 gene) and matches it strongly
# (1e-69), but hexA is a derivative of eIF5A. In2's genome carries both hex1 and
# its own eIF5A (eif2). hexA's within-genome paralog (eif1) hits In2's eif2 even
# harder (1e-70) than hexA hits hex1 — but only on a *different* target gene.
HEXA_HITS = (
    'hexA\thex1\t1e-69\t230\tIn1\tIn2\n'   # hexA -> hex1: the real ortholog hit
    'eif1\teif2\t1e-70\t233\tIn1\tIn2\n'   # eIF5A paralog -> In2's eIF5A: wins proteome-wide
    'eif1\thex1\t5e-12\t45\tIn1\tIn2\n'    # eIF5A paralog -> hex1: loses on this target
)
# hexA's within-genome paralog is eif1; eif1's paralog is hexA. Only filter 2 (paralog
# competition) is exercised here -- filter 1 is the flat default and both hits (1e-69,
# 1e-70) easily clear it.
HEXA_PARALOGS = 'hexA\teif1\t42\t4.2e-11\neif1\thexA\t42\t4.2e-11\n'


def test_competition_proteome_scope_drops_hexa_like_ortholog(run_dir):
    # 'proteome' scope: the eIF5A paralog out-scores hexA *anywhere* in In2, so hexA's
    # only qualifying hit is dropped. With no surviving cross-hit hexA gets no matrix row
    # at all -> never a candidate. The rescue floor is disabled here to isolate the scope
    # behaviour -- see test_default_floor_protects_the_hexa_ortholog_under_proteome_scope
    # for what the shipped default does to this same case.
    matrix, candidates = run(run_dir, HEXA_HITS, paralog_text=HEXA_PARALOGS,
                             competition_scope='proteome', rescue_evalue='0')
    assert (matrix['protein_id'] == 'hexA').sum() == 0
    assert 'In1::hexA' not in candidates


def test_default_floor_protects_the_hexa_ortholog_under_proteome_scope(run_dir):
    # Behaviour change from issue #128. hexA -> hex1 is a real ortholog hit at 1e-69,
    # which clears the 1e-20 default floor, so it is no longer discarded even under the
    # strict 'proteome' scope. This is the floor doing its job: 'proteome' scope's known
    # failure mode was exactly this -- discarding a real ortholog because a paralog won
    # somewhere else in the target genome.
    matrix, candidates = run(run_dir, HEXA_HITS, paralog_text=HEXA_PARALOGS,
                             competition_scope='proteome')
    row = matrix[matrix['protein_id'] == 'hexA'].iloc[0]
    assert row['In1'] == 1 and row['In2'] == 1
    assert 'In1::hexA' in candidates


def test_competition_target_scope_keeps_hexa_like_ortholog(run_dir):
    # 'target' scope: on the shared target gene hex1, hexA (1e-69) beats its
    # paralog eif1 (5e-12), so the call survives -> present in 2/2 ingroup -> candidate.
    matrix, candidates = run(run_dir, HEXA_HITS, paralog_text=HEXA_PARALOGS,
                             competition_scope='target', rescue_evalue='0')
    row = matrix[matrix['protein_id'] == 'hexA'].iloc[0]
    assert row['In1'] == 1 and row['In2'] == 1
    assert 'In1::hexA' in candidates


# --- Filter 2 rescue (--paralog-rescue-evalue) ------------------------------
#
# Mirrors the pezizo_set1 A7UWR3/NCU11312-vs-HEX1 case: A7UWR3 (In1) hits Out1's
# target at 1e-95 -- independently strong evidence of homology -- but its
# in-genome paralog Q7SEN2 (here `strongParalog`) hits the *same* Out1 target
# even harder (1e-131), so 'target'-scope filter 2 disqualifies A7UWR3's hit on
# its own. HEX1 (weakQuery), by contrast, only reaches Out1 at a marginal 1e-9,
# with its own paralog (IF5A) again winning head-to-head -- there's no
# independently strong signal to rescue.
RESCUE_HITS = (
    'a7uwr3\ttarget1\t1e-95\t310\tIn1\tOut1\n'   # strong on its own; paralog still wins here
    'q7sen2\ttarget1\t1e-131\t409\tIn1\tOut1\n'  # in-genome paralog: wins head-to-head
    'hex1\ttarget2\t1e-9\t53\tIn1\tOut1\n'       # marginal; paralog wins by a wide margin too
    'if5a\ttarget2\t1e-69\t204\tIn1\tOut1\n'
)
RESCUE_PARALOGS = (
    'a7uwr3\tq7sen2\t303\t1e-92\nq7sen2\ta7uwr3\t297\t1e-90\n'
    'hex1\tif5a\t44\t2.8e-6\nif5a\thex1\t39\t1.1e-4\n'
)


def test_target_scope_disqualifies_a7uwr3_like_hit_without_rescue(run_dir):
    # With the floor explicitly disabled (0), a7uwr3's only hit is disqualified outright
    # -> no qualifying hit anywhere -> no matrix row at all for it (same "fully
    # disqualified" shape as test_flat_default_evalue_rejects_a_weak_hit above). This is
    # the pre-2026-09-20 behaviour, kept reachable via --paralog-rescue-evalue 0.
    matrix, candidates = run(run_dir, RESCUE_HITS, paralog_text=RESCUE_PARALOGS,
                             competition_scope='target', rescue_evalue='0')
    assert (matrix['protein_id'] == 'a7uwr3').sum() == 0
    assert 'In1::a7uwr3' not in candidates


def test_paralog_rescue_evalue_keeps_a7uwr3_like_hit(run_dir):
    # 1e-20 floor: a7uwr3's own 1e-95 clears it, so its hit survives filter 2 despite
    # q7sen2 winning head-to-head -- a7uwr3 is no longer called a lineage-specific novelty.
    matrix, candidates = run(run_dir, RESCUE_HITS, paralog_text=RESCUE_PARALOGS,
                             competition_scope='target', rescue_evalue='1e-20')
    row = matrix[matrix['protein_id'] == 'a7uwr3'].iloc[0]
    assert row['Out1'] == 1
    assert 'In1::a7uwr3' not in candidates  # present in Out1 now -> no longer an IN-only novelty


def test_paralog_rescue_evalue_still_drops_hex1_like_marginal_hit(run_dir):
    # Same 1e-20 floor: hex1's own hit is only 1e-9, well above the floor, so it still
    # doesn't clear the rescue and stays disqualified -- unlike a7uwr3, there's no
    # independently strong signal here for the floor to protect, so hex1's Out1 hit is
    # dropped exactly as it is without --paralog-rescue-evalue at all.
    matrix, candidates = run(run_dir, RESCUE_HITS, paralog_text=RESCUE_PARALOGS,
                             competition_scope='target', rescue_evalue='1e-20')
    assert (matrix['protein_id'] == 'hex1').sum() == 0
    assert 'In1::hex1' not in candidates


# --- The floor arm is ON by default (issue #128) -----------------------------
#
# Measured on NovInvenio_Investigations/results/pezizo_set1 (3393 candidates), the
# 1e-20 floor removes 213 candidates, 57.3% of which have TBLASTN hits in >=4 of the
# 6 outgroup GENOMES -- against a 48.2% rate for removing every filter-2-suppressed
# candidate indiscriminately. It is the only rule measured that beat that baseline.


def test_floor_arm_is_enabled_by_default(run_dir):
    # No --paralog-rescue-evalue flag at all: the 1e-20 default applies, so the
    # a7uwr3-like 1e-95 hit is rescued and the protein is not called a novelty.
    matrix, candidates = run(run_dir, RESCUE_HITS, paralog_text=RESCUE_PARALOGS,
                             competition_scope='target')
    assert matrix[matrix['protein_id'] == 'a7uwr3'].iloc[0]['Out1'] == 1
    assert 'In1::a7uwr3' not in candidates


def test_default_floor_still_drops_the_hex1_like_marginal_hit(run_dir):
    # The positive control must survive the new default untouched: hex1's own hit is
    # 1e-9, well above the 1e-20 floor, so filter 2 still removes it and hex1 keeps
    # its clean outgroup-absent row.
    matrix, _ = run(run_dir, RESCUE_HITS, paralog_text=RESCUE_PARALOGS,
                    competition_scope='target')
    assert (matrix['protein_id'] == 'hex1').sum() == 0


def test_rescue_evalue_zero_disables_the_floor(run_dir):
    # 0 is the off switch -- no e-value is <= 0, so the rescue never fires. This is how
    # nextflow.config's `paralog_rescue_evalue = null` reaches the script.
    matrix, _ = run(run_dir, RESCUE_HITS, paralog_text=RESCUE_PARALOGS,
                    competition_scope='target', rescue_evalue='0')
    assert (matrix['protein_id'] == 'a7uwr3').sum() == 0


# --- Filter 2 rescue, delta arm (--paralog-rescue-delta) --------------------
#
# delta = log10(query_evalue) - log10(paralog_evalue): how many orders of magnitude
# the paralog beat the query by. Asks "did the paralog explain the hit away, or only
# narrowly win?" rather than "is the hit strong in absolute terms?".
#
# OFF BY DEFAULT, and deliberately so. Measured on pezizo_set1, a delta rule is NOT
# selective: because the rescue fires per cell and any one rescued cell kills novelty,
# the trigger is the MINIMUM delta across a candidate's suppressed cells, whose median
# is only 8.6 -- so `delta < 45` fires for 291 of the 334 suppressed candidates, barely
# narrower than removing all 334. Its breadth>=4 hit-rate (47.1%) is below the 48.2%
# indiscriminate baseline, and OR-ing it onto the floor drags the floor's 57.3% down to
# 49.0%. Kept as an opt-in flag so the rule can be re-swept against better evidence
# (see issue #129, alignment coverage), not because it is recommended.
#
# RESCUE_HITS supplies both controls: a7uwr3 delta 36, hex1 delta 60, so D=45 separates
# them -- in the same direction as the real data (A7UWR3 cells 2..38, HEX1 cells 56/59).


def test_delta_arm_is_disabled_by_default(run_dir):
    # Floor off, no delta flag: nothing rescues the a7uwr3-like hit.
    matrix, _ = run(run_dir, RESCUE_HITS, paralog_text=RESCUE_PARALOGS,
                    competition_scope='target', rescue_evalue='0')
    assert (matrix['protein_id'] == 'a7uwr3').sum() == 0


def test_paralog_rescue_delta_keeps_a7uwr3_like_hit(run_dir):
    # delta 36 < 45 rescues it even with the floor explicitly off.
    matrix, candidates = run(run_dir, RESCUE_HITS, paralog_text=RESCUE_PARALOGS,
                             competition_scope='target', rescue_evalue='0', rescue_delta='45')
    assert matrix[matrix['protein_id'] == 'a7uwr3'].iloc[0]['Out1'] == 1
    assert 'In1::a7uwr3' not in candidates


def test_paralog_rescue_delta_still_drops_hex1_like_marginal_hit(run_dir):
    # delta 60 >= 45: the paralog beat hex1 by a wide margin, which is what a genuine
    # "this hit belongs to the ancestral gene" case looks like. hex1 stays disqualified.
    matrix, _ = run(run_dir, RESCUE_HITS, paralog_text=RESCUE_PARALOGS,
                    competition_scope='target', rescue_evalue='0', rescue_delta='45')
    assert (matrix['protein_id'] == 'hex1').sum() == 0


# The two arms cover different cells, so when both are on they must be OR-ed, not AND-ed.
#   strongFloor  q=1e-95  p=1e-200  -> delta 105: clears a 1e-20 floor, fails delta 45.
#   nearDelta    q=1e-13  p=1e-15   -> delta 2:   fails a 1e-20 floor, clears delta 45.
# nearDelta is the shape of A7UWR3's real Ccin cell (q=4.6e-13, p=3.9e-15, delta=2).
OR_HITS = (
    'strongFloor\ttargetA\t1e-95\t310\tIn1\tOut1\n'
    'sfParalog\ttargetA\t1e-200\t600\tIn1\tOut1\n'
    'nearDelta\ttargetB\t1e-13\t60\tIn1\tOut1\n'
    'ndParalog\ttargetB\t1e-15\t70\tIn1\tOut1\n'
)
OR_PARALOGS = (
    'strongFloor\tsfParalog\t300\t1e-90\nsfParalog\tstrongFloor\t300\t1e-90\n'
    'nearDelta\tndParalog\t80\t1e-20\nndParalog\tnearDelta\t80\t1e-20\n'
)


def test_floor_arm_alone_rescues_only_the_strong_hit(run_dir):
    matrix, _ = run(run_dir, OR_HITS, paralog_text=OR_PARALOGS,
                    competition_scope='target', rescue_evalue='1e-20')
    assert matrix[matrix['protein_id'] == 'strongFloor'].iloc[0]['Out1'] == 1
    assert (matrix['protein_id'] == 'nearDelta').sum() == 0


def test_delta_arm_alone_rescues_only_the_narrow_margin_hit(run_dir):
    matrix, _ = run(run_dir, OR_HITS, paralog_text=OR_PARALOGS,
                    competition_scope='target', rescue_evalue='0', rescue_delta='45')
    assert matrix[matrix['protein_id'] == 'nearDelta'].iloc[0]['Out1'] == 1
    assert (matrix['protein_id'] == 'strongFloor').sum() == 0


def test_rescue_arms_are_ored_not_anded(run_dir):
    # Each hit is rescued by the one arm that covers it. If the arms were AND-ed,
    # neither would be rescued and both rows would vanish.
    matrix, _ = run(run_dir, OR_HITS, paralog_text=OR_PARALOGS,
                    competition_scope='target', rescue_evalue='1e-20', rescue_delta='45')
    assert matrix[matrix['protein_id'] == 'strongFloor'].iloc[0]['Out1'] == 1
    assert matrix[matrix['protein_id'] == 'nearDelta'].iloc[0]['Out1'] == 1


# diamond reports 0.0 for an overwhelming hit. log10(0) is -inf, so a zero-e-value
# paralog beat the query by an unbounded margin: delta is +inf and the hit must stay
# disqualified, with no nan and no crash.
ZERO_PARALOG_HITS = (
    'q\ttargetZ\t1e-95\t310\tIn1\tOut1\n'
    'p\ttargetZ\t0.0\t900\tIn1\tOut1\n'
)
ZERO_PARALOG_PARALOGS = 'q\tp\t300\t1e-90\np\tq\t300\t1e-90\n'


def test_delta_arm_handles_a_zero_paralog_evalue(run_dir):
    matrix, _ = run(run_dir, ZERO_PARALOG_HITS, paralog_text=ZERO_PARALOG_PARALOGS,
                    competition_scope='target', rescue_evalue='0', rescue_delta='45')
    assert (matrix['protein_id'] == 'q').sum() == 0


def test_output_evalues_sidecar_matches_presence_calls(run_dir):
    # issue #44: --output-evalues emits the qualifying hit's e-value alongside each
    # presence=1 cell, empty for absence and for the protein's own source proteome.
    hits_path = run_dir / 'hits.tsv'
    hits_path.write_text(HIT_HEADER + 'g1\tg1_in2\t1e-10\t100\tIn1\tIn2\n')
    matrix_out = run_dir / 'matrix.tsv'
    candidates_out = run_dir / 'candidates.txt'
    evalues_out = run_dir / 'evalues.tsv'
    subprocess.run([
        sys.executable, str(SCRIPT),
        '--hits', str(hits_path),
        '--config', str(run_dir / 'config.csv'),
        '--ingroup-min-frac', '1.0',
        '--query-group', 'IN',
        '--other-max-frac', '0.0',
        '--output-matrix', str(matrix_out),
        '--output-candidates', str(candidates_out),
        '--output-evalues', str(evalues_out),
    ], check=True, capture_output=True, text=True)

    evalues = pd.read_csv(evalues_out, sep='\t', dtype=str, keep_default_na=False)
    row = evalues[evalues['protein_id'] == 'g1'].iloc[0]
    assert row['In2'] == '1e-10'
    assert row['In1'] == ''   # source proteome: presence is definitional, not a hit
    assert row['Out1'] == ''  # absent: no qualifying hit


def test_output_targets_sidecar_names_the_winning_hits_target_id(run_dir):
    # --output-targets mirrors --output-evalues' shape, but records the *target_id*
    # of the best (lowest-evalue) qualifying hit instead of its e-value -- so the
    # report can resolve and show which target protein a presence call came from.
    hits_path = run_dir / 'hits.tsv'
    hits_path.write_text(
        HIT_HEADER +
        'g1\tg1_in2_weak\t1e-8\t80\tIn1\tIn2\n'
        'g1\tg1_in2_strong\t1e-30\t300\tIn1\tIn2\n'
    )
    matrix_out = run_dir / 'matrix.tsv'
    candidates_out = run_dir / 'candidates.txt'
    evalues_out = run_dir / 'evalues.tsv'
    targets_out = run_dir / 'targets.tsv'
    subprocess.run([
        sys.executable, str(SCRIPT),
        '--hits', str(hits_path),
        '--config', str(run_dir / 'config.csv'),
        '--ingroup-min-frac', '1.0',
        '--query-group', 'IN',
        '--other-max-frac', '0.0',
        '--output-matrix', str(matrix_out),
        '--output-candidates', str(candidates_out),
        '--output-evalues', str(evalues_out),
        '--output-targets', str(targets_out),
    ], check=True, capture_output=True, text=True)

    targets = pd.read_csv(targets_out, sep='\t', dtype=str, keep_default_na=False)
    row = targets[targets['protein_id'] == 'g1'].iloc[0]
    assert row['In2'] == 'g1_in2_strong'  # the lower-evalue hit's target, not the weaker one
    assert row['In1'] == ''
    assert row['Out1'] == ''
