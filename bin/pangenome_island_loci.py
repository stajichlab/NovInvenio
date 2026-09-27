#!/usr/bin/env python3
"""Compute the island locus view (spec
docs/superpowers/specs/2026-09-24-island-locus-view-design.md, sections 1-6)
and write it as island_loci.json for bin/pangenome_island_synteny.py.

family_positions.tsv[.zst] is streamed three times, each pass filtered:
  1. copies of the candidate loci's member families + every contig's
     (min, max) rank, to find carriers and choose each exemplar;
  2. the exemplars' own contigs, to read their gene order (the columns);
  3. copies of every column family, for the per-strain cell states.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import open_maybe_compressed  # noqa: E402
from config_parser import parse_config  # noqa: E402
from island_locus import (  # noqa: E402
    DEFAULT_CONTAINMENT, DEFAULT_EMPTY_FRAC, DEFAULT_FLANK, DEFAULT_FLANK_MIN, DEFAULT_K,
    candidate_loci, carrier_placements, choose_exemplar, compute_locus, group_loci,
    locus_columns, locus_payload, rank_key,
)
from island_locus import (  # noqa: E402
    DEFAULT_DNA_MIN_COV, DEFAULT_DNA_MIN_ID, apply_dna_calls, dna_checked_strains, dna_query,
    dna_target, rescue_hit_spans,
)
from pangenome_matrix import PresenceMatrix  # noqa: E402
from pfam_classes import dominant_class  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from pangenome_domain_enrichment import parse_domtblout  # noqa: E402
from pangenome_island_synteny import load_gene_locations  # noqa: E402
from island_synteny import column_locations  # noqa: E402


@dataclass
class PositionScan:
    positions: dict = field(default_factory=dict)  # {(strain, family): [(contig, rank)]}
    orders: dict = field(default_factory=dict)  # {(strain, contig): [(rank, family)]}
    spans: dict = field(default_factory=dict)  # {(strain, contig): (min_rank, max_rank)}


def scan_family_positions(path: str, families: set[str] | None = None,
                          contigs: set[tuple[str, str]] | None = None,
                          spans: bool = False) -> PositionScan:
    """One streaming pass over family_positions.tsv[.zst] (Short, family,
    contig, rank). Keeps every copy of a family in `families`, the full
    order of each (strain, contig) in `contigs`, and with `spans` the
    (min, max) rank of every (strain, contig)."""
    out = PositionScan()
    with open_maybe_compressed(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_s, i_f, i_c, i_r = (header.index(c) for c in ("Short", "family", "contig", "rank"))
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_s, i_f, i_c, i_r):
                continue
            try:
                rank = int(parts[i_r])
            except ValueError:
                continue
            strain, fam, contig = parts[i_s], parts[i_f], parts[i_c]
            if families is not None and fam in families:
                out.positions.setdefault((strain, fam), []).append((contig, rank))
            if contigs is not None and (strain, contig) in contigs:
                out.orders.setdefault((strain, contig), []).append((rank, fam))
            if spans:
                key = (strain, contig)
                cur = out.spans.get(key)
                if cur is None:
                    out.spans[key] = (rank, rank)
                elif rank < cur[0] or rank > cur[1]:
                    out.spans[key] = (min(cur[0], rank), max(cur[1], rank))
    return out


def read_islands(path: str) -> list[dict]:
    with open_maybe_compressed(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_matrix_strains(path: str) -> list[str]:
    with open_maybe_compressed(path) as fh:
        return fh.readline().rstrip("\n").split("\t")[1:]


def read_n50(path: str | None) -> dict[str, int]:
    """{Short: n50} from assembly_quality_vs_content.tsv. Missing or empty
    file gives {}."""
    if not path or not Path(path).is_file() or Path(path).stat().st_size == 0:
        return {}
    out = {}
    with open_maybe_compressed(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                out[row["Short"]] = int(row["n50"])
            except (KeyError, ValueError):
                continue
    return out


def read_bins(path: str | None) -> dict[str, str]:
    if not path:
        return {}
    with open_maybe_compressed(path) as fh:
        return {r["family"]: r.get("bin", "") for r in csv.DictReader(fh, delimiter="\t")}


def build(args) -> dict:
    islands = read_islands(args.islands_with_domains)
    strains = read_matrix_strains(args.presence_matrix)
    loci = group_loci(islands, args.containment)
    cands = candidate_loci(loci, len(strains), args.rank_by, args.candidates, args.min_strains)
    n50 = read_n50(args.assembly_quality)
    species_of = ({s.short: s.species for s in parse_config(args.config)} if args.config else {})

    members = {f for loc in cands for f in loc.members}
    scan1 = scan_family_positions(args.family_positions, families=members, spans=True)
    chosen = []
    n_unplaced = 0
    for loc in cands:
        places = carrier_placements(loc.members, scan1.positions, strains, scan1.spans, args.k)
        pick = choose_exemplar(places, n50, args.flank, args.flank_min)
        if pick is None:
            n_unplaced += 1
            continue
        chosen.append((loc, pick[0], pick[1]))

    scan2 = scan_family_positions(args.family_positions,
                                  contigs={(p.strain, p.contig) for _, p, _ in chosen})
    with_cols = []
    for loc, place, tier in chosen:
        flank = args.flank_min if tier == "short_flanks" else args.flank
        cols = locus_columns(scan2.orders.get((place.strain, place.contig), []), loc.members,
                             flank, block=(place.lo, place.hi))
        if cols is None:
            n_unplaced += 1
            continue
        with_cols.append((loc, place, tier, cols))

    col_fams = {f for _, _, _, cols in with_cols for f in cols.families}
    scan3 = scan_family_positions(args.family_positions, families=col_fams)
    matrix = PresenceMatrix.from_tsv(args.presence_matrix, families=col_fams)
    results = [compute_locus(loc, place, tier, cols, strains, scan3.positions, scan1.spans,
                             matrix, species_of, args.k, args.empty_frac)
               for loc, place, tier, cols in with_cols]
    results.sort(key=lambda r: rank_key(r, args.rank_by))
    drawn = results[:args.top_loci]

    # DNA presence check (spec section 4b, plan Ruling R17): pass 1 writes
    # the work list for ISLAND_DNA_CHECK; pass 2 (the same inputs, so the
    # same drawn loci) applies its calls and re-orders the drawn loci.
    if args.dna_targets_dir:
        if not (args.gene_positions and args.cluster_tsv):
            raise SystemExit("--dna_targets_dir needs --gene_positions and --cluster_tsv")
        exemplar_ranks = {
            loc.locus_id: [r for r, _ in sorted(scan2.orders.get((place.strain, place.contig), []))
                           if place.lo <= r <= place.hi]
            for loc, place, _, _ in with_cols}
        dna_locs = load_gene_locations(
            args.cluster_tsv, args.gene_positions, {f for r in drawn for f in r["families"]},
            {s for r in drawn for s in dna_checked_strains(r)} | {r["exemplar"] for r in drawn},
            args.id_sep)
        # Spec 4b, Ruling R26: an exemplar locus column that is a TBLASTN
        # rescue hit takes the hit's span (rescue_positions start + the end
        # of that HSP in the exemplar's own tblastn output).
        rescue_spans = load_rescue_spans(
            args.rescue_positions, args.rescue_tblastn,
            {(r["exemplar"], f) for r in drawn
             for f in r["families"][r["n_left"]:r["n_left"] + r["n_locus"]]})
        write_dna_targets(args.dna_targets_dir,
                          dna_target_rows(drawn, exemplar_ranks, scan3.positions, dna_locs, args.k,
                                          rescue_spans),
                          args.dna_batch)
    dna_on = args.dna_check == "true"
    if dna_on:
        calls = read_dna_calls(args.dna_calls)
        for r in drawn:
            apply_dna_calls(r, calls.get(r["locus_id"], {}), species_of, args.empty_frac)
        drawn.sort(key=lambda r: rank_key(r, args.rank_by))

    bins = read_bins(args.frequency_table)
    fam_domains = (parse_domtblout([args.domtblout], max_ievalue=args.domain_evalue)
                   if args.domtblout else {})
    gene_locs = None
    if args.gene_positions and args.cluster_tsv:
        gene_locs = load_gene_locations(args.cluster_tsv, args.gene_positions,
                                        {f for r in drawn for f in r["families"]},
                                        {r["exemplar"] for r in drawn}, args.id_sep)
    out_loci = []
    for i, r in enumerate(drawn):
        dom_sets = [sorted(fam_domains.get(f, set())) for f in r["families"]]
        locs = span = None
        if gene_locs is not None:
            locs = column_locations(r["families"], r["exemplar"], r["exemplar_contig"], -1, -1,
                                    gene_locs)
            block = [x for x in locs[r["n_left"]:r["n_left"] + r["n_locus"]]
                     if x and not x.get("rescued")]
            if block:
                span = (min(x["start"] for x in block), max(x["end"] for x in block))
        out_loci.append(locus_payload(
            r, "L%03d" % (i + 1), bins,
            [dominant_class(d) for d in dom_sets], [",".join(d) for d in dom_sets],
            dominant_class(sorted({d for ds in dom_sets for d in ds})), locs, span))
    return {
        "project": args.project,
        "locus_params": {"flank": args.flank, "flank_min": args.flank_min, "k": args.k,
                         "empty_frac": args.empty_frac, "containment": args.containment,
                         "rank_by": args.rank_by, "top_loci": args.top_loci,
                         "candidates": args.candidates, "min_strains": args.min_strains,
                         "dna_check": dna_on, "dna_min_id": args.dna_min_id,
                         "dna_min_cov": args.dna_min_cov},
        "n_loci_total": len(loci),
        "n_loci_candidates": len(cands),
        "n_loci_unplaced": n_unplaced,
        "n_strains_with_n50": sum(1 for s in strains if s in n50),
        "n_strains": len(strains),
        "loci": out_loci,
        "_results": drawn,
    }


DNA_TARGET_COLUMNS = ["locus_id", "role", "strain", "contig", "start", "end", "genes"]
DNA_CALL_COLUMNS = ["locus_id", "strain", "col", "status", "coverage"]


def load_rescue_spans(rescue_positions: str | None, tblastn_paths: list[str] | None,
                      wanted: set[tuple[str, str]]) -> dict[tuple[str, str], list[tuple]]:
    """Spans of the TBLASTN rescue hits of the `wanted` (strain, family)
    pairs (Ruling R26): rescue_positions.tsv gives (contig, start); the
    strain's own tblastn file (`<Short>.tblastn.tsv[.gz|.zst]`, matched by
    file name, so only the files of strains with a wanted row are read)
    gives the end. Missing or empty inputs give no spans."""
    rows = []
    if rescue_positions and Path(rescue_positions).is_file() \
            and Path(rescue_positions).stat().st_size > 0:
        with open_maybe_compressed(rescue_positions) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if (row.get("Short"), row.get("family")) not in wanted:
                    continue
                try:
                    rows.append((row["Short"], row["family"], row["contig"], int(row["start"])))
                except (KeyError, ValueError):
                    continue
    files: dict[str, list[str]] = {}
    for path in tblastn_paths or []:
        files.setdefault(Path(path).name.split(".tblastn")[0], []).append(path)
    spans: dict[tuple[str, str], list[tuple]] = {}
    for strain in sorted({x[0] for x in rows}):
        mine = [x for x in rows if x[0] == strain]
        for path in files.get(strain, []):
            with open_maybe_compressed(path) as fh:
                for key, val in rescue_hit_spans(mine, fh).items():
                    spans.setdefault(key, []).extend(val)
    n_spans = sum(len(v) for v in spans.values())
    print(f"pangenome_island_loci: {len(rows)} rescue-only exemplar locus copies, "
          f"{n_spans} with a TBLASTN span, {len(rows) - n_spans} without (unchecked)",
          file=sys.stderr)
    return spans


def dna_target_rows(drawn: list[dict], exemplar_ranks: dict[str, list[int]],
                    positions: dict, gene_locs: dict, k: int,
                    rescue_spans: dict | None = None) -> list[list[list]]:
    """The DNA check's work list, one block per drawn locus (spec 4b): a
    query row (the exemplar's locus DNA; genes = 'column:start-end;...') and
    a target row per checked strain. A strain without a target interval gets
    contig '-' and start = end = 0 (Ruling R20). A rescue-only exemplar
    column uses its `rescue_spans` span (Ruling R26). A locus with no
    checked strain, or no exemplar locus column with a span, has no block."""
    blocks = []
    for r in drawn:
        strains = dna_checked_strains(r)
        query = dna_query(r, exemplar_ranks.get(r["locus_id"], []), positions, gene_locs,
                          rescue_spans)
        if not strains or query is None:
            continue
        contig, start, end, genes = query
        rows = [[r["locus_id"], "query", r["exemplar"], contig, start, end,
                 ";".join(f"{c}:{a}-{b}" for c, a, b in genes)]]
        for s in strains:
            target = dna_target(r, s, positions, gene_locs, k)
            rows.append([r["locus_id"], "target", s, *(target or ("-", 0, 0)), ""])
        blocks.append(rows)
    return blocks


def write_dna_targets(out_dir: str, blocks: list[list[list]], batch: int) -> list[Path]:
    """batch_001.tsv, batch_002.tsv, ... in `out_dir`, `batch` loci each
    (one ISLAND_DNA_CHECK task per file). No blocks, no files."""
    d = Path(out_dir)
    d.mkdir(parents=True, exist_ok=True)
    paths = []
    step = max(1, batch)
    for i in range(0, len(blocks), step):
        path = d / ("batch_%03d.tsv" % (i // step + 1))
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(DNA_TARGET_COLUMNS)
            for rows in blocks[i:i + step]:
                w.writerows(rows)
        paths.append(path)
    return paths


def read_dna_calls(paths: list[str]) -> dict[str, dict[str, dict[int, str]]]:
    """{locus_id: {strain: {column: status}}} from ISLAND_DNA_CHECK's
    dna_calls_*.tsv files. Missing or empty files are skipped."""
    out: dict[str, dict[str, dict[int, str]]] = {}
    for path in paths or []:
        if not Path(path).is_file() or Path(path).stat().st_size == 0:
            continue
        with open_maybe_compressed(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                out.setdefault(row["locus_id"], {}).setdefault(row["strain"], {})[
                    int(row["col"])] = row["status"]
    return out


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--islands_with_domains", required=True)
    ap.add_argument("--presence_matrix", required=True)
    ap.add_argument("--family_positions", required=True)
    ap.add_argument("--project", required=True)
    ap.add_argument("--frequency_table", default=None)
    ap.add_argument("--assembly_quality", default=None,
                    help="assembly_quality_vs_content.tsv (Short, n50, ...) for exemplar ties")
    ap.add_argument("--config", default=None, help="samplesheet CSV, for species")
    ap.add_argument("--domtblout", default=None)
    ap.add_argument("--domain_evalue", type=float, default=1e-3)
    ap.add_argument("--gene_positions", default=None)
    ap.add_argument("--cluster_tsv", default=None)
    ap.add_argument("--id_sep", default="|")
    ap.add_argument("--flank", type=int, default=DEFAULT_FLANK)
    ap.add_argument("--flank_min", type=int, default=DEFAULT_FLANK_MIN)
    ap.add_argument("--k", type=int, default=DEFAULT_K)
    ap.add_argument("--empty_frac", type=float, default=DEFAULT_EMPTY_FRAC)
    ap.add_argument("--containment", type=float, default=DEFAULT_CONTAINMENT)
    ap.add_argument("--rank_by", choices=["informative", "strains", "size"], default="informative")
    ap.add_argument("--top_loci", type=int, default=50)
    ap.add_argument("--candidates", type=int, default=200)
    ap.add_argument("--dna_targets_dir", default=None,
                    help="pass 1 of the DNA presence check: write batch_NNN.tsv work lists here")
    ap.add_argument("--dna_batch", type=int, default=50, help="loci per DNA check batch file")
    ap.add_argument("--dna_check", choices=["true", "false"], default="false",
                    help="pass 2: apply --dna_calls (spec section 4b)")
    ap.add_argument("--dna_calls", nargs="*", default=[],
                    help="dna_calls_*.tsv from bin/pangenome_island_dna_check.py")
    ap.add_argument("--dna_min_id", type=float, default=DEFAULT_DNA_MIN_ID,
                    help="recorded in locus_params for the page note")
    ap.add_argument("--dna_min_cov", type=float, default=DEFAULT_DNA_MIN_COV,
                    help="recorded in locus_params for the page note")
    ap.add_argument("--rescue_positions", default=None,
                    help="pass 1: rescue_positions.tsv (Short, family, contig, start)")
    ap.add_argument("--rescue_tblastn", nargs="*", default=[],
                    help="pass 1: per-strain <Short>.tblastn.tsv[.zst] files, for the end of "
                    "a rescue-only exemplar column's hit (Ruling R26)")
    ap.add_argument("--min_strains", type=int, default=2)
    ap.add_argument("--output", required=True)
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    payload = build(args)
    payload.pop("_results")
    Path(args.output).write_text(json.dumps(payload, separators=(",", ":")))
    print(f"pangenome_island_loci: {payload['n_loci_total']} loci, "
          f"{payload['n_loci_candidates']} candidates, {payload['n_loci_unplaced']} without "
          f"a carrier, {len(payload['loci'])} drawn (--rank_by {args.rank_by}); wrote "
          f"{args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
