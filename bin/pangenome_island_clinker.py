#!/usr/bin/env python3
"""ISLAND_CLINKER: one clinker figure per locus directory (spec section 8).

For each --locus_dirs entry (<key>/ holding <strain>.gbk files, an optional
clinker_order.txt and groups.csv from bin/pangenome_island_gbk_slice.py) runs

    clinker <key>/*.gbk -gf <key>/groups.csv -j <cpus> -p <tmp>.html -ufo

(-ufo/--use_file_order, M7: the .gbk files are given in clinker pick order --
exemplar first -- from clinker_order.txt, not clinker's own default
alignment-based ordering) and writes <out_dir>/<key>.html, slimmed by
lib/clinker_html.py unless --keep_sequences.

A missing clinker executable (I2a) fails the whole run loudly (exit 2)
instead of silently skipping every locus -- see main(). A single locus
whose clinker run fails, or whose output page is malformed (M1), is
skipped with a warning and no page instead, so one bad locus never costs
the others (plan Ruling R15); the synteny page then says "No synteny
figure for this locus". Exit status is otherwise 0.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from clinker_html import inject_ui_fixes, slim_clinker_html  # noqa: E402


def ordered_gbks(locus_dir: Path) -> list[str]:
    """This locus's .gbk file paths in clinker pick order (M7 -- exemplar
    first): bin/pangenome_island_gbk_slice.py writes clinker_order.txt with
    one file name per line in that order. Falls back to alphabetical when
    the order file is missing (an older ISLAND_GBK_SLICE run) or lists a
    file that is no longer there."""
    order_file = locus_dir / "clinker_order.txt"
    present = {p.name: str(p) for p in locus_dir.glob("*.gbk")}
    if order_file.is_file():
        names = order_file.read_text().split()
        ordered = [present[name] for name in names if name in present]
        ordered += [p for name, p in sorted(present.items()) if name not in set(names)]
        return ordered
    return [present[name] for name in sorted(present)]


def run_locus(locus_dir: Path, out_dir: Path, clinker: str, cpus: int, keep: bool) -> bool:
    key = locus_dir.name
    gbks = ordered_gbks(locus_dir)
    if not gbks:
        print(f"WARNING: {key}: no .gbk files; skipped", file=sys.stderr)
        return False
    raw = out_dir / f"{key}.raw.html"
    cmd = [clinker, *gbks, "-gf", str(locus_dir / "groups.csv"), "-j", str(cpus), "-p", str(raw),
          "-ufo"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except OSError as e:
        print(f"WARNING: {key}: cannot run {clinker}: {e}; skipped", file=sys.stderr)
        return False
    if proc.returncode != 0 or not raw.is_file():
        print(f"WARNING: {key}: clinker exited {proc.returncode}; skipped\n{proc.stderr[-2000:]}",
              file=sys.stderr)
        raw.unlink(missing_ok=True)
        return False
    html = raw.read_text()
    try:
        page = inject_ui_fixes(html) if keep else slim_clinker_html(html)
    except ValueError as e:
        # M1: a clinker page with no "const data=" block (a version whose
        # markup differs, or a truncated/corrupted -p output) is a
        # per-locus failure like any other (R15), not a fatal error.
        print(f"WARNING: {key}: malformed clinker page ({e}); skipped", file=sys.stderr)
        raw.unlink(missing_ok=True)
        return False
    (out_dir / f"{key}.html").write_text(page)
    raw.unlink()
    return True


def clinker_missing(clinker: str) -> bool:
    """True unless `clinker` resolves to a real executable, either by name
    on PATH (shutil.which()) or as an existing executable file path."""
    if shutil.which(clinker):
        return False
    path = Path(clinker)
    return not (path.is_file() and os.access(path, os.X_OK))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--locus_dirs", nargs="+", required=True)
    ap.add_argument("--out_dir", default=".")
    ap.add_argument("--cpus", type=int, default=1)
    ap.add_argument("--clinker", default="clinker", help="clinker executable")
    ap.add_argument("--keep_sequences", action="store_true",
                    help="write clinker's page unchanged (--pangenome_clinker_slim false)")
    args = ap.parse_args(argv)
    # I2a (final review): a missing clinker executable would otherwise make
    # every locus "fail" individually (R15's per-locus skip), which reads
    # like a real per-locus clinker problem instead of the actual cause --
    # the published container image (0.5.0) has no clinker in it. Fail the
    # whole run loudly instead.
    if clinker_missing(args.clinker):
        print(f"ERROR: clinker executable '{args.clinker}' not found (the published "
              "container image 0.5.0 has no clinker; rebuild the image or use the pixi "
              "environment, or run with --pangenome_clinker false)", file=sys.stderr)
        return 2
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    done = sum(run_locus(Path(d), out_dir, args.clinker, args.cpus, args.keep_sequences)
               for d in args.locus_dirs)
    print(f"pangenome_island_clinker: {done} of {len(args.locus_dirs)} loci drawn", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
