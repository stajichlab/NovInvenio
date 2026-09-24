"""Every bin/ script a Nextflow process calls must be executable.

Nextflow puts bin/ on PATH and runs scripts by bare name, so a script committed
without the exec bit fails only at run time, with "Permission denied" (exit 126).
Unit tests import these scripts, so they do not catch it. Found 2026-09-24:
pangenome_assembly_quality_qc.py and pangenome_diagnostics.py were mode 100644,
and a 529-strain Coccidioides run failed ASSEMBLY_QUALITY_QC after 7.8 h.
"""
import os
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
NF_FILES = [p for p in REPO.rglob("*.nf") if ".nextflow" not in p.parts and "work" not in p.parts]


def test_bin_scripts_called_from_nextflow_are_executable():
    bin_names = {p.name for p in (REPO / "bin").iterdir() if p.is_file()}
    called = set()
    for nf in NF_FILES:
        text = nf.read_text()
        for name in bin_names:
            if re.search(rf"(?<![\w./-]){re.escape(name)}(?![\w.-])", text):
                called.add(name)
    assert called, "no bin/ scripts found in any .nf file -- test is not looking in the right place"
    not_exec = sorted(n for n in called if not os.access(REPO / "bin" / n, os.X_OK))
    assert not not_exec, f"bin/ scripts called from .nf files but not executable: {not_exec}"
