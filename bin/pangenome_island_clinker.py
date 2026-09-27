#!/usr/bin/env python3
"""ISLAND_CLINKER: one clinker figure per locus directory (spec section 8).

For each --locus_dirs entry (<key>/ holding <strain>.gbk files and
groups.csv from bin/pangenome_island_gbk_slice.py) runs

    clinker <key>/*.gbk -gf <key>/groups.csv -j <cpus> -p <tmp>.html

and writes <out_dir>/<key>.html, slimmed by lib/clinker_html.py unless
--keep_sequences. A locus whose clinker run fails (or a missing clinker
executable) is skipped with a warning and no page, so one bad locus never
costs the others (plan Ruling R15); the synteny page then says "No synteny
figure for this locus". Exit status is 0 unless an argument is wrong.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from clinker_html import inject_ui_fixes, slim_clinker_html  # noqa: E402


def run_locus(locus_dir: Path, out_dir: Path, clinker: str, cpus: int, keep: bool) -> bool:
    key = locus_dir.name
    gbks = sorted(str(p) for p in locus_dir.glob("*.gbk"))
    if not gbks:
        print(f"WARNING: {key}: no .gbk files; skipped", file=sys.stderr)
        return False
    raw = out_dir / f"{key}.raw.html"
    cmd = [clinker, *gbks, "-gf", str(locus_dir / "groups.csv"), "-j", str(cpus), "-p", str(raw)]
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
    (out_dir / f"{key}.html").write_text(
        inject_ui_fixes(html) if keep else slim_clinker_html(html))
    raw.unlink()
    return True


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--locus_dirs", nargs="+", required=True)
    ap.add_argument("--out_dir", default=".")
    ap.add_argument("--cpus", type=int, default=1)
    ap.add_argument("--clinker", default="clinker", help="clinker executable")
    ap.add_argument("--keep_sequences", action="store_true",
                    help="write clinker's page unchanged (--pangenome_clinker_slim false)")
    args = ap.parse_args(argv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    done = sum(run_locus(Path(d), out_dir, args.clinker, args.cpus, args.keep_sequences)
               for d in args.locus_dirs)
    print(f"pangenome_island_clinker: {done} of {len(args.locus_dirs)} loci drawn", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
