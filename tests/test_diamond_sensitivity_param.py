"""--diamond_sensitivity wiring (todo/diamond-very-sensitive-main-search.md, issue #135).

Structural checks on the Nextflow side. There is no cheap way to unit-test a Groovy
script block from Python, so these assert the three properties that matter and that a
future edit could silently break:

  1. The param exists with a default that PRESERVES current behaviour (diamond's default
     fast mode). This matters more than usual: adopting a more sensitive mode as the
     default invalidates results for every existing study, so it must be an explicit,
     recorded decision -- never a side effect of adding the flag.
  2. DIAMOND_SEARCH actually emits the flag when set, and NOTHING when unset (an empty
     string interpolated into the command line must not leave a stray token).
  3. A typo fails at pipeline START via main.nf's existing validation block, not hours
     later inside a queued task after the search step has finally been scheduled.

Why the flag exists at all: on pezizo_set1, 841 of 1479 candidate x outgroup cells with
genome-level (TBLASTN) evidence but no protein-level evidence get ZERO diamond hits --
not even weak noise to E=0.01 -- at loci where the outgroup protein exists in the FASTA.
See issue #135.
"""
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CONFIG = (REPO / 'nextflow.config').read_text()
DIAMOND = (REPO / 'modules' / 'diamond.nf').read_text()
MAIN = (REPO / 'main.nf').read_text()

# diamond's own sensitivity modes, weakest to strongest. 'fast' is diamond's default
# and is spelled by leaving the param empty, so it is deliberately not in this list.
MODES = ['sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive']


def _search_script_block():
    """The DIAMOND_SEARCH process only -- DIAMOND_MAKEDB above it must not be touched."""
    return DIAMOND.split('process DIAMOND_SEARCH', 1)[1]


def test_param_is_declared_and_defaults_to_diamonds_own_fast_mode():
    m = re.search(r"diamond_sensitivity\s*=\s*(['\"])(.*?)\1", CONFIG)
    assert m, 'params.diamond_sensitivity is not declared in nextflow.config'
    assert m.group(2) == '', (
        'the default must be empty (diamond fast mode) so adding this flag changes no '
        'existing result; adopting a sensitive default is a separate, recorded decision')


def test_diamond_search_passes_the_flag_through_when_set():
    block = _search_script_block()
    assert 'params.diamond_sensitivity' in block


def test_diamond_search_emits_nothing_when_unset():
    # The interpolation must be conditional: a bare `--${...}` would emit a lone `--`
    # when the param is empty, which diamond rejects.
    block = _search_script_block()
    assert re.search(r"params\.diamond_sensitivity\s*\?", block), (
        'the flag must be emitted conditionally (ternary/elvis) so an empty param '
        'contributes no token to the diamond command line')


def test_makedb_is_not_affected():
    # Sensitivity is a search-time option; `diamond makedb` does not accept it, and
    # DIAMOND_MAKEDB is storeDir-cached across runs.
    makedb = DIAMOND.split('process DIAMOND_SEARCH', 1)[0]
    assert 'diamond_sensitivity' not in makedb


def test_a_typo_fails_at_pipeline_start():
    assert 'diamond_sensitivity' in MAIN, (
        'main.nf has no validation for --diamond_sensitivity; a typo would only surface '
        'inside a queued task, after the search step is finally scheduled')
    for mode in MODES:
        assert mode in MAIN, f"main.nf validation does not accept '{mode}'"


def test_self_search_stays_at_very_sensitive():
    # DIAMOND_SELF is hardcoded --very-sensitive on purpose (HEX-1/eIF-5A paralog
    # detection, see modules/self_search.nf). This param must never loosen it.
    self_search = (REPO / 'modules' / 'self_search.nf').read_text()
    assert '--very-sensitive' in self_search
    assert 'diamond_sensitivity' not in self_search
