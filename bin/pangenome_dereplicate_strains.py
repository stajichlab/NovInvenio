#!/usr/bin/env python3
"""Strain inventory: per-strain assembly-quality proxy (N50/contig count)
and mash-based dereplication, feeding the frequency/co-occurrence steps'
"dereplicated strain count" denominator.

Ported from NovInvenio_Investigations' Afumigatus_pangenome study
(bin/dereplicate_strains.py). Two changes from that copy:

1. `config_parser` is imported directly from this repo's own lib/ (no
   cross-repo NOVINVENIO_ROOT lookup needed).
2. `--mash_dist_tsv` (optional): pass a precomputed `mash dist -t` matrix
   (e.g. from a single, shared MASH_SKETCH Nextflow process) instead of
   having this script run `mash sketch`/`mash dist` itself. The study's own
   NEXTFLOW_MIGRATION_NOTES.md flagged this and pangenome_assign_clades.py
   as independently re-sketching the exact same genome set -- real
   duplicated work; both scripts can now share one upstream sketch/dist
   run. Omitting `--mash_dist_tsv` preserves the original standalone
   behavior (this script sketches itself), so it still works as a
   plain CLI tool outside Nextflow.

Usage:
  pangenome_dereplicate_strains.py --config config.csv --data_dir data_dir \\
      --output strain_inventory.tsv
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def compute_assembly_stats(fasta_path: str | Path) -> dict:
    """Return {n_contigs, n50, total_length} for a DNA FASTA, without
    external tools -- a fast, dependency-free fragmentation proxy."""
    lengths = []
    current = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
        if current:
            lengths.append(current)
    lengths.sort(reverse=True)
    total = sum(lengths)
    half = total / 2
    running = 0
    n50 = 0
    for length in lengths:
        running += length
        if running >= half:
            n50 = length
            break
    return {"n_contigs": len(lengths), "n50": n50, "total_length": total}


def parse_mash_dist(lines: list[str], threshold: float) -> list[set[str]]:
    """Parse `mash dist -t` output (a query x reference distance matrix,
    tab-separated, first row is '#query\\tref1\\tref2\\t...') into
    dereplication groups: strains whose pairwise mash distance is below
    `threshold` are grouped together (union-find over the threshold graph)."""
    if not lines:
        return []
    header = lines[0].lstrip("#").split("\t")
    names = header[1:]
    parent = {name: name for name in names}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for row in lines[1:]:
        parts = row.split("\t")
        query = parts[0]
        for ref, dist_str in zip(names, parts[1:]):
            if query == ref:
                continue
            if float(dist_str) < threshold:
                union(query, ref)

    groups: dict[str, set[str]] = {}
    for name in names:
        root = find(name)
        groups.setdefault(root, set()).add(name)
    return list(groups.values())


def groups_to_shorts(
    groups: list[set[str]], path_to_short: dict[str, str]
) -> list[set[str]]:
    """Translate mash's path-string dedup groups into config Short IDs.

    mash reports whatever path string it was actually invoked with. Matched
    by BASENAME, not the full path string: when a shared upstream Nextflow
    MASH_SKETCH process does the sketching (see
    modules/pangenome/mash.nf -- this script's `--mash_dist_tsv` path), the
    genome FASTAs are Nextflow-staged files in that process's own work
    directory, so mash reports their bare filenames, not
    "<data_dir>/dna/<file>". Matching by basename works identically whether
    this script sketches itself (the standalone-CLI fallback path, where
    mash would report the full `data_dir/dna/...` path -- its basename is
    still unique per strain) or consumes a pre-sketched matrix, as long as
    every strain's DNA filename is unique (already required: the samplesheet
    keys everything on Short + a per-strain DNA filename).

    `path_to_short` is keyed by basename here (see `main()`).

    Raises:
        KeyError: with the offending path and an example of an expected one --
            a bare KeyError here means a strain's DNA filename in the
            samplesheet doesn't match what was actually sketched.
    """
    out: list[set[str]] = []
    for group in groups:
        shorts = set()
        for path in group:
            key = Path(path).name
            if key not in path_to_short:
                example = next(iter(path_to_short), "<none>")
                raise KeyError(
                    f"mash reported path {path!r} (basename {key!r}), which is not "
                    f"one of the {len(path_to_short)} samplesheet DNA filenames "
                    f"(e.g. {example!r}). Every strain's DNA filename must be unique "
                    "and must match what was actually sketched."
                )
            shorts.add(path_to_short[key])
        out.append(shorts)
    return out


def choose_representatives(
    dedup_groups: list[set[str]], assembly_stats: dict[str, dict]
) -> dict[str, str]:
    """For every strain in every group, map it to the group's chosen
    representative (the member with the highest N50 -- the least
    fragmented assembly)."""
    result: dict[str, str] = {}
    for group in dedup_groups:
        # Use deterministic tie-break on equal N50 (lexically last Short ID)
        rep = max(group, key=lambda s: (assembly_stats[s]["n50"], s))
        for member in group:
            result[member] = rep
    return result


def sketch_and_dist(dna_paths: dict[str, Path], sketch_prefix: str) -> str:
    """Run `mash sketch` + `mash dist -t` over `dna_paths` and return the
    dist matrix's stdout text. The fallback path used when no
    `--mash_dist_tsv` was supplied (i.e. this script is run standalone,
    outside a Nextflow MASH_SKETCH process)."""
    subprocess.run(
        ["mash", "sketch", "-o", sketch_prefix] + [str(p) for p in dna_paths.values()],
        check=True,
    )
    return subprocess.run(
        ["mash", "dist", "-t", f"{sketch_prefix}.msh", f"{sketch_prefix}.msh"],
        check=True, capture_output=True, text=True,
    ).stdout


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--data_dir", required=True)
    ap.add_argument("--mash_threshold", type=float, default=0.001)
    ap.add_argument(
        "--mash_dist_tsv",
        help="precomputed 'mash dist -t' output (e.g. from a shared MASH_SKETCH "
        "process); when omitted, this script runs mash sketch/dist itself",
    )
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
    from config_parser import parse_config  # noqa: E402

    samples = parse_config(args.config)
    data_dir = Path(args.data_dir)
    dna_paths = {s.short: data_dir / "dna" / s.dna for s in samples if s.dna}

    stats = {short: compute_assembly_stats(p) for short, p in dna_paths.items()}

    if args.mash_dist_tsv:
        dist_out = Path(args.mash_dist_tsv).read_text()
    else:
        dist_out = sketch_and_dist(dna_paths, "strain_sketches")
    groups = parse_mash_dist(dist_out.splitlines(), args.mash_threshold)

    # Translate groups from mash-reported paths to Short IDs, matching by
    # basename (see groups_to_shorts' docstring for why).
    path_to_short = {p.name: short for short, p in dna_paths.items()}
    if len(path_to_short) != len(dna_paths):
        # Two different strains whose DNA files share a basename (different
        # directories) collapse into one dict entry here, silently dropping
        # a strain from mash-distance grouping even though this function's
        # docstring assumes basenames are unique.
        by_basename: dict[str, list[str]] = {}
        for short, p in dna_paths.items():
            by_basename.setdefault(p.name, []).append(short)
        collisions = {name: shorts for name, shorts in by_basename.items() if len(shorts) > 1}
        raise ValueError(
            f"Duplicate DNA basenames across strains: {collisions} -- "
            "each strain's DNA filename must be unique (across its full "
            "samplesheet directory), since dereplication matches mash's "
            "reported paths back to Short IDs by basename alone."
        )
    groups = groups_to_shorts(groups, path_to_short)

    reps = choose_representatives(groups, stats)
    dedup_group_id = {}
    for i, group in enumerate(groups):
        for member in group:
            dedup_group_id[member] = i

    with open(args.output, "w") as fh:
        fh.write("Short\tn_contigs\tn50\ttotal_length\tdedup_group\tis_representative\n")
        for short in sorted(dna_paths):
            s = stats[short]
            fh.write(
                f"{short}\t{s['n_contigs']}\t{s['n50']}\t{s['total_length']}\t"
                f"{dedup_group_id[short]}\t{int(reps[short] == short)}\n"
            )


if __name__ == "__main__":
    main()
