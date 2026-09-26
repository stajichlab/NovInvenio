"""Issue #189: BUILD_PRESENCE_MATRIX only had `label 'low_cpu'`'s flat 4 GB memory
directive, with no scaling across retries, so a study large enough to OOM it hit
the identical wall on every attempt (exit 137) and failed the whole run instead
of recovering.

There is no cheap way to unit-test a `nextflow.config` `process {}` block from
Python without actually running Nextflow (which would need a real OOM to
trigger), so -- following the pattern already used for --diamond_sensitivity
(tests/test_diamond_sensitivity_param.py) and documented in the brief as
acceptable when nothing better is feasible -- these are structural checks on
the config text plus a confirmation that modules/build_presence_matrix.nf's
label was left untouched (the fix must not change it) and that the retry
policy already in effect for it (conf/ucr_hpcc_slurm.config) actually retries
on any non-zero exit, including OOM's 137.
"""
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CONFIG = (REPO / 'nextflow.config').read_text()
SLURM_CONFIG = (REPO / 'conf' / 'ucr_hpcc_slurm.config').read_text()
MODULE = (REPO / 'modules' / 'build_presence_matrix.nf').read_text()


def _withname_block(text, pattern):
    """Return the withName block's body, up to its OWN closing brace -- not the
    first '}' encountered, which may belong to a nested closure like
    `memory = { 16.GB * task.attempt }` inside the block."""
    start_m = re.search(r"withName:\s*'" + re.escape(pattern) + r"'\s*\{", text)
    assert start_m, f"no withName: '{pattern}' block found"
    depth = 1
    i = start_m.end()
    body_start = i
    while depth > 0:
        if text[i] == '{':
            depth += 1
        elif text[i] == '}':
            depth -= 1
        i += 1
    return text[body_start:i - 1]


def test_module_label_is_unchanged():
    """The fix must scale memory via nextflow.config's process selector, not by
    editing the process's own label directive."""
    m = re.search(r"process BUILD_PRESENCE_MATRIX\s*\{\s*label\s+'([^']+)'", MODULE)
    assert m, 'BUILD_PRESENCE_MATRIX should still declare a label directive'
    assert m.group(1) == 'low_cpu', (
        "modules/build_presence_matrix.nf's label must stay 'low_cpu' -- "
        "the memory override belongs in nextflow.config's withName block")


def test_withname_override_exists_and_scales_with_attempt():
    block = _withname_block(CONFIG, '.*BUILD_PRESENCE_MATRIX')
    m = re.search(r"memory\s*=\s*\{\s*(\d+)\.GB\s*\*\s*task\.attempt\s*\}", block)
    assert m, (
        "withName: '.*BUILD_PRESENCE_MATRIX' must set memory as a task.attempt-"
        "scaled closure (same pattern as ASSEMBLY_QUALITY_QC/HMMSEARCH_CHUNK), "
        "not a flat value -- a flat value hits the identical OOM wall on every retry")
    base_gb = int(m.group(1))
    assert base_gb >= 16, (
        f"base memory is {base_gb}GB, not enough headroom above low_cpu's flat "
        "4GB to expect a retry to actually succeed")


def test_withname_override_comes_after_low_cpu_label_default():
    """withName selectors override withLabel in Nextflow's directive precedence
    regardless of position, but keep the override textually after low_cpu's
    definition so the file continues to document the fixed-vs-scaled contrast
    for a reader scanning top-to-bottom (matches ASSEMBLY_QUALITY_QC's placement)."""
    low_cpu_idx = CONFIG.index("withLabel: 'low_cpu'")
    override_idx = CONFIG.index("withName: '.*BUILD_PRESENCE_MATRIX'")
    assert override_idx > low_cpu_idx


def test_slurm_retry_policy_is_unconditional_on_exit_status():
    """errorStrategy must not be gated on a specific exit code (e.g. only retrying
    137) -- confirm it retries on ANY non-zero exit, so OOM's 137 is covered
    without needing a bespoke clause for BUILD_PRESENCE_MATRIX."""
    assert re.search(r"errorStrategy\s*=\s*\{\s*task\.attempt\s*<\s*\d+\s*\?\s*'retry'",
                      SLURM_CONFIG), (
        "expected an unconditional (not exit-status-gated) retry errorStrategy "
        "in conf/ucr_hpcc_slurm.config")
    m = re.search(r"maxRetries\s*=\s*(\d+)", SLURM_CONFIG)
    assert m and int(m.group(1)) >= 2, (
        "maxRetries must allow at least 2 retries so a task.attempt-scaled "
        "memory override actually gets a chance to try a bigger allocation")
