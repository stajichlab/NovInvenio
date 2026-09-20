import shutil
import signal
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# Make the project's `lib/` importable in tests.
sys.path.insert(0, str(REPO / 'lib'))
# ...and this directory, so test modules can `from conftest import tool_runs`.
sys.path.insert(0, str(Path(__file__).resolve().parent))

# How the pipeline's own beforeScript resolves tools for a local (`-profile standard`)
# run. Probed here so a skip reflects the binary that would actually execute.
PIXI_BIN = REPO / '.pixi' / 'envs' / 'default' / 'bin'

# Enough to tell "this CPU cannot run the binary" from "the binary ran and complained".
_SIGILL_CODES = {132, -signal.SIGILL}


def tool_runs(cmd) -> bool:
    """True if `cmd` executes at all on this machine.

    A non-zero exit is still True -- several of these tools exit non-zero when given no
    arguments on a perfectly healthy node (famsa prints usage), and treating that as
    "cannot run" would skip the integration tests everywhere, silently. Only a missing
    binary or SIGILL counts as unrunnable.

    SIGILL is the signature of a binary compiled for instructions this CPU lacks --
    bioconda's mmseqs2 and famsa are built with AVX2 and die instantly on this cluster's
    Abu Dhabi nodes (see .living/learnings.md L-3).
    """
    try:
        p = subprocess.run(cmd, capture_output=True, timeout=60)
    except (FileNotFoundError, PermissionError, OSError):
        return False
    except subprocess.TimeoutExpired:
        return True      # it ran; being slow is not our concern here
    return p.returncode not in _SIGILL_CODES


def _resolve(tool: str) -> str | None:
    local = PIXI_BIN / tool
    if local.exists():
        return str(local)
    return shutil.which(tool)


def unrunnable_tools(*tools: str) -> list[str]:
    """Which of `tools` are missing or cannot execute here. Probe arguments are chosen
    to be harmless and fast; the exit code is ignored except for SIGILL."""
    probe_args = {'mmseqs': ['version'], 'famsa': [], 'mafft': ['--version'],
                  'hmmbuild': ['-h'], 'hmmsearch': ['-h'], 'diamond': ['version']}
    bad = []
    for t in tools:
        path = _resolve(t)
        if path is None or not tool_runs([path] + probe_args.get(t, [])):
            bad.append(t)
    return bad
