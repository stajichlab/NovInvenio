#!/usr/bin/env python3
"""ISLAND_GBK_SLICE: one GenBank file per (locus, strain) region listed in
island_regions.tsv (bin/pangenome_island_loci.py --regions_out), for clinker.

Writes, under --out_dir:
  <locus_key>/<strain>.gbk   gene + CDS features, /locus_tag, /translation,
                             /note="family=<tier-1 family>"
  <locus_key>/groups.csv     locus_tag,family (clinker -gf; no header)
  island_slices.tsv          locus_key, strain, contig, bp_start, bp_end,
                             n_genes, n_rescue, n_missing
Exits non-zero if a region's anchor families do not match the families at
those ranks after the rank rebuild (lib/genbank_slice.rank_entries), since
that means the regions and the slices describe different genes.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from Bio import SeqIO  # noqa: E402
from compressed_io import open_maybe_compressed  # noqa: E402
from config_parser import parse_config  # noqa: E402
from genbank_slice import build_record, cds_exons, rank_entries, safe_name  # noqa: E402

DNA_SUBDIRS = ["dna", "genome", "scaffolds"]
PEP_SUBDIRS = ["pep", "proteins"]


def resolve_under(data_dir: Path, name: str, subdirs: list[str]) -> Path | None:
    """Flat layout first, then each subdir -- pangenome.nf's resolve_fa()."""
    for cand in [data_dir / name] + [data_dir / sub / name for sub in subdirs]:
        if cand.exists():
            return cand
    return None


def read_regions(path: str) -> list[dict]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["rank_lo"], r["rank_hi"] = int(r["rank_lo"]), int(r["rank_hi"])
        r["anchor_list"] = [(int(a.split(":", 1)[0]), a.split(":", 1)[1])
                            for a in r["anchors"].split(";") if a]
    return rows


def iter_gene_rows(path: str):
    with open_maybe_compressed(path) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            try:
                yield p[0], p[1], p[2], int(p[3]), int(p[4])
            except ValueError:
                continue


def iter_rescue_rows(path: str | None):
    if not path or not Path(path).is_file() or Path(path).stat().st_size == 0:
        return
    with open_maybe_compressed(path) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            try:
                yield p[0], p[1], p[2], int(p[3])
            except ValueError:
                continue


def read_families(cluster_tsv: str, members: set[str]) -> dict[str, str]:
    """{member_id: family} for the wanted `Short<sep>protein` member IDs."""
    out = {}
    with open_maybe_compressed(cluster_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1] in members:
                out[parts[1]] = parts[0]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--data_dir", required=True)
    ap.add_argument("--gff3_dir", required=True)
    ap.add_argument("--gene_positions", required=True)
    ap.add_argument("--rescue_positions", default=None)
    ap.add_argument("--cluster_tsv", required=True)
    ap.add_argument("--id_sep", default="|")
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args(argv)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    regions = read_regions(args.regions)
    slices_path = out_dir / "island_slices.tsv"
    header = ["locus_key", "strain", "contig", "bp_start", "bp_end", "n_genes", "n_rescue",
              "n_missing"]
    if not regions:
        slices_path.write_text("\t".join(header) + "\n")
        print("pangenome_island_gbk_slice: no regions", file=sys.stderr)
        return 0

    samples = {s.short: s for s in parse_config(args.config)}
    wanted = {(r["strain"], r["contig"]) for r in regions}
    index = rank_entries(iter_gene_rows(args.gene_positions),
                         iter_rescue_rows(args.rescue_positions), wanted)
    region_entries = []
    for r in regions:
        entries = [e for e in index.get((r["strain"], r["contig"]), [])
                   if r["rank_lo"] <= e[0] <= r["rank_hi"]]
        region_entries.append(entries)
    members = {f"{r['strain']}{args.id_sep}{e[1]}"
               for r, ents in zip(regions, region_entries) for e in ents if e[1]}
    fam_of_member = read_families(args.cluster_tsv, members)

    errors = []
    for r, ents in zip(regions, region_entries):
        at = {e[0]: (fam_of_member.get(f"{r['strain']}{args.id_sep}{e[1]}") if e[1] else e[2])
              for e in ents}
        for rank, fam in r["anchor_list"]:
            if at.get(rank) != fam:
                errors.append(f"{r['locus_key']} {r['strain']} rank {rank}: region says "
                              f"{fam}, rank rebuild gives {at.get(rank)}")
    if errors:
        print("ERROR: rank rebuild does not match island_regions.tsv:\n  " +
              "\n  ".join(errors[:20]), file=sys.stderr)
        return 1

    by_strain: dict[str, list[int]] = {}
    for i, r in enumerate(regions):
        by_strain.setdefault(r["strain"], []).append(i)
    rows_out = []
    group_rows: dict[str, dict[str, str]] = {}
    pid_owner: dict[str, dict[str, set[str]]] = {}
    for i, r in enumerate(regions):
        for e in region_entries[i]:
            if e[1]:
                pid_owner.setdefault(r["locus_key"], {}).setdefault(e[1], set()).add(r["strain"])
    for strain, idxs in sorted(by_strain.items()):
        s = samples.get(strain)
        if s is None:
            print(f"WARNING: {strain} not in --config; skipped", file=sys.stderr)
            continue
        dna = resolve_under(Path(args.data_dir), s.dna, DNA_SUBDIRS)
        pep = resolve_under(Path(args.data_dir), s.protein, PEP_SUBDIRS)
        gff = Path(args.gff3_dir) / s.gff3 if s.gff3 else None
        if dna is None or pep is None or gff is None or not gff.exists():
            print(f"WARNING: {strain}: DNA, protein or GFF3 file not found; skipped",
                  file=sys.stderr)
            continue
        contigs = {regions[i]["contig"] for i in idxs}
        pids = {e[1] for i in idxs for e in region_entries[i] if e[1]}
        with open_maybe_compressed(str(dna)) as fh:
            seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(fh, "fasta") if rec.id in contigs}
        with open_maybe_compressed(str(pep)) as fh:
            prots = {rec.id: str(rec.seq) for rec in SeqIO.parse(fh, "fasta") if rec.id in pids}
        with open_maybe_compressed(str(gff)) as fh:
            exons = cds_exons(fh, contigs, pids)
        for i in idxs:
            r = regions[i]
            ents = region_entries[i]
            if not ents or r["contig"] not in seqs:
                print(f"WARNING: {r['locus_key']} {strain}: no genes or contig "
                      f"{r['contig']} not in {dna}; skipped", file=sys.stderr)
                continue
            clash = any(len(owners) > 1 for owners in pid_owner.get(r["locus_key"], {}).values())
            family_of = {e[1]: fam_of_member.get(f"{strain}{args.id_sep}{e[1]}", "")
                         for e in ents if e[1]}
            rec, summ = build_record(
                strain, r["contig"], seqs[r["contig"]], ents, exons, prots, family_of,
                (lambda pid, st=strain: f"{st}{args.id_sep}{pid}") if clash else (lambda pid: pid))
            locus_dir = out_dir / r["locus_key"]
            locus_dir.mkdir(exist_ok=True)
            SeqIO.write(rec, str(locus_dir / f"{safe_name(strain)}.gbk"), "genbank")
            group_rows.setdefault(r["locus_key"], {}).update(summ["labels"])
            rows_out.append([r["locus_key"], strain, r["contig"], summ["bp_start"],
                             summ["bp_end"], summ["n_genes"], summ["n_rescue"], summ["n_missing"]])
    for key, labels in group_rows.items():
        with open(out_dir / key / "groups.csv", "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            for label, fam in sorted(labels.items()):
                w.writerow([label, fam])
    with open(slices_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(sorted(rows_out))
    print(f"pangenome_island_gbk_slice: {len(rows_out)} slices for {len(group_rows)} loci",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
