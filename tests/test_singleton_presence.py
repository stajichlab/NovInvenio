"""Filter-2 rescue arms on the singleton pathway (issue #138).

lib/singleton_presence.py applies the same paralog-competition filter as
bin/build_presence_matrix.py, so it needs the same rescue arms -- otherwise the same
protein gets two different novelty calls depending on --cluster_tool (#136 vs here).

Fixtures mirror tests/test_build_presence_matrix.py's RESCUE_HITS so the two pathways are
checked against the same two curated controls:
  a7uwr3  q=1e-95  paralog=1e-131  -> delta 36   (real A7UWR3 cells span 2..38)
  hex1    q=1e-9   paralog=1e-69   -> delta 60   (real HEX1 cells are 56 and 59)
"""
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'lib'))
from singleton_presence import score_singleton_hits  # noqa: E402

# (query_id, target_id, evalue, proteome_short)
HITS = [
    ('a7uwr3', 'target1', 1e-95,  'Out1'),   # strong on its own; paralog still wins here
    ('q7sen2', 'target1', 1e-131, 'Out1'),   # in-genome paralog: wins head-to-head
    ('hex1',   'target2', 1e-9,   'Out1'),   # marginal; paralog wins by a wide margin
    ('if5a',   'target2', 1e-69,  'Out1'),
]
SINGLETONS = {'a7uwr3', 'hex1'}
PARALOG_OF = {'a7uwr3': 'q7sen2', 'q7sen2': 'a7uwr3', 'hex1': 'if5a', 'if5a': 'hex1'}


def score(**kw):
    kw.setdefault('competition_scope', 'target')
    presence, _ = score_singleton_hits(HITS, SINGLETONS, PARALOG_OF, 1e-5, **kw)
    return presence['Out1']


def test_floor_is_enabled_by_default():
    # Matches bin/build_presence_matrix.py's default: a7uwr3's own 1e-95 clears the
    # 1e-20 floor, so its hit survives filter 2 and it is not called outgroup-absent.
    assert 'a7uwr3' in score()


def test_default_floor_still_drops_the_hex1_like_marginal_hit():
    # The positive control must survive the default untouched: 1e-9 is well above the
    # floor, so filter 2 still removes it.
    assert 'hex1' not in score()


def test_rescue_evalue_zero_disables_the_floor():
    # 0 is the off switch -- no e-value is <= 0 -- restoring the pre-#138 behaviour.
    assert 'a7uwr3' not in score(rescue_evalue=0)


def test_delta_arm_is_disabled_by_default():
    assert 'a7uwr3' not in score(rescue_evalue=0)


def test_delta_arm_rescues_the_a7uwr3_like_hit():
    # delta 36 < 45, so the paralog did not explain the hit away.
    assert 'a7uwr3' in score(rescue_evalue=0, rescue_delta=45)


def test_delta_arm_still_drops_the_hex1_like_marginal_hit():
    # delta 60 >= 45: a genuine "this hit belongs to the ancestral gene" case.
    assert 'hex1' not in score(rescue_evalue=0, rescue_delta=45)


# Each arm covers a cell the other misses, so when both are on they must be OR-ed.
#   strongFloor  q=1e-95 p=1e-200 -> delta 105: clears a 1e-20 floor, fails delta 45
#   nearDelta    q=1e-13 p=1e-15  -> delta 2:   fails the floor, clears delta 45
OR_HITS = [
    ('strongFloor', 'targetA', 1e-95,  'Out1'),
    ('sfParalog',   'targetA', 1e-200, 'Out1'),
    ('nearDelta',   'targetB', 1e-13,  'Out1'),
    ('ndParalog',   'targetB', 1e-15,  'Out1'),
]
OR_SINGLETONS = {'strongFloor', 'nearDelta'}
OR_PARALOG_OF = {'strongFloor': 'sfParalog', 'sfParalog': 'strongFloor',
                 'nearDelta': 'ndParalog', 'ndParalog': 'nearDelta'}


def or_score(**kw):
    kw.setdefault('competition_scope', 'target')
    presence, _ = score_singleton_hits(OR_HITS, OR_SINGLETONS, OR_PARALOG_OF, 1e-5, **kw)
    return presence['Out1']


def test_rescue_arms_are_ored_not_anded():
    got = or_score(rescue_evalue=1e-20, rescue_delta=45)
    assert 'strongFloor' in got and 'nearDelta' in got


def test_floor_arm_alone_rescues_only_the_strong_hit():
    got = or_score(rescue_evalue=1e-20)
    assert 'strongFloor' in got and 'nearDelta' not in got


def test_delta_arm_alone_rescues_only_the_narrow_margin_hit():
    got = or_score(rescue_evalue=0, rescue_delta=45)
    assert 'nearDelta' in got and 'strongFloor' not in got


ZERO_HITS = [('q', 'targetZ', 1e-95, 'Out1'), ('p', 'targetZ', 0.0, 'Out1')]


def test_delta_arm_handles_a_zero_paralog_evalue():
    # log10(0) is -inf, so delta is +inf: the paralog beat the query by an unbounded
    # margin and the hit must stay disqualified, with no nan and no crash.
    presence, _ = score_singleton_hits(
        ZERO_HITS, {'q'}, {'q': 'p', 'p': 'q'}, 1e-5,
        competition_scope='target', rescue_evalue=0, rescue_delta=45)
    assert 'q' not in presence['Out1']
