#!/usr/bin/env python3
"""Per-module Pfam domain summary for `pangenome_detect_trans_modules.py`'s
Leiden trans-module output: one row per module (with >= --min_module_size
members) listing how many member families have >=1 Pfam domain hit and
which domain names those hits are.

This is deliberately general -- it does not hardcode any specific domain
class (NACHT/HET, PKS/NRPS, etc.). A targeted search for a specific domain
family (e.g. "does this study have any HET-domain-carrying modules?") is a
one-line `grep`/`awk` over this script's output, not a pipeline parameter --
see CHANGES.md for a worked example from the coccidioides_pangenome study.

Usage:
  pangenome_module_domains.py --family_modules family_modules.tsv \\
      --domtblout pfam.domtblout --min_module_size 2 \\
      --output module_domains.tsv
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import open_maybe_compressed  # noqa: E402


def load_family_modules(path: str) -> dict[str, tuple[str, int]]:
    """{family: (module_id, module_size)}"""
    out = {}
    with open_maybe_compressed(path) as fh:
        next(fh)
        for line in fh:
            family, module_id, module_size = line.rstrip("\n").split("\t")
            out[family] = (module_id, int(module_size))
    return out


def parse_domtblout(path: str) -> dict[str, set[tuple[str, str]]]:
    """{query_protein: {(domain_name, accession), ...}} from a hmmscan
    --domtblout file (target_name=domain, accession, tlen, query_name=protein, ...)."""
    hits: dict[str, set[tuple[str, str]]] = defaultdict(set)
    with open_maybe_compressed(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            domain_name, accession, query = parts[0], parts[1], parts[3]
            hits[query].add((domain_name, accession))
    return hits


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family_modules", required=True)
    ap.add_argument("--domtblout", required=True)
    ap.add_argument("--min_module_size", type=int, default=2)
    ap.add_argument("--top_n_domains", type=int, default=15,
                     help="max domain names listed per module (most frequent first)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    family_module = load_family_modules(args.family_modules)
    domain_hits = parse_domtblout(args.domtblout) if Path(args.domtblout).stat().st_size else {}
    print(f"{len(family_module)} families in modules, {len(domain_hits)} proteins with >=1 Pfam hit",
          file=sys.stderr)

    module_families: dict[str, list[str]] = defaultdict(list)
    for family, (module_id, size) in family_module.items():
        if size >= args.min_module_size:
            module_families[module_id].append(family)

    rows = []
    for module_id, families in module_families.items():
        domains: Counter = Counter()
        n_with = 0
        for fam in families:
            fam_domains = {d for d, _ in domain_hits.get(fam, set())}
            if fam_domains:
                n_with += 1
            domains.update(fam_domains)
        top_domains = ",".join(f"{d}({c})" for d, c in domains.most_common(args.top_n_domains))
        rows.append((module_id, len(families), n_with, top_domains))

    with open(args.output, "w") as out:
        out.write("module_id\tmodule_size\tn_families_with_domain\tdomain_names\n")
        for module_id, size, n_with, top_domains in sorted(rows, key=lambda r: -r[1]):
            out.write(f"{module_id}\t{size}\t{n_with}\t{top_domains or '-'}\n")

    print(f"pangenome_module_domains: {len(rows)} modules (>= {args.min_module_size} members) summarized",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
