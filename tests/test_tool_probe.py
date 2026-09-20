"""The runnable-tool probe behind the integration tests' skip (issue #142).

The probe exists because a suite that is reliably red for an unrelated reason is a suite
whose red stops being read. On this cluster's Abu Dhabi nodes both mmseqs and famsa
SIGILL (no AVX2), so the pipeline integration test fails for a reason that has nothing
to do with the code under test.

Probing beats grepping /proc/cpuinfo for 'avx2': it is self-correcting (working binaries
on PATH make the test run, with no flag to remember), it is not tied to one instruction
set, and it checks the binary that will actually be used rather than inferring from CPU
flags.

SIGILL is simulated here rather than depending on this node's hardware, so the test means
the same thing everywhere.
"""
import sys

from conftest import tool_runs

PY = sys.executable
SIGILL_CMD = [PY, '-c', 'import os, signal; os.kill(os.getpid(), signal.SIGILL)']
NONZERO_CMD = [PY, '-c', 'raise SystemExit(1)']
OK_CMD = [PY, '-c', 'pass']


def test_a_working_command_is_runnable():
    assert tool_runs(OK_CMD) is True


def test_a_command_killed_by_sigill_is_not_runnable():
    # The real signature on a no-AVX2 node: the binary exists and is executable, but the
    # CPU cannot run its instructions.
    assert tool_runs(SIGILL_CMD) is False


def test_a_merely_nonzero_exit_is_still_runnable():
    # Crucial distinction: `famsa` with no arguments exits non-zero on a HEALTHY node
    # (it prints usage). Treating any non-zero exit as "cannot run" would skip the
    # integration test everywhere, silently, which is worse than the false red.
    assert tool_runs(NONZERO_CMD) is True


def test_a_missing_binary_is_not_runnable():
    assert tool_runs(['definitely-not-a-real-binary-xyzzy']) is False
