#!/usr/bin/env python3
"""DNA presence check for the island locus view (spec
docs/superpowers/specs/2026-09-24-island-locus-view-design.md, section 4b).

Input: one batch_NNN.tsv work list from bin/pangenome_island_loci.py
--dna_targets_dir (columns locus_id, role, strain, contig, start, end,
genes). Per locus, the query row is the exemplar's locus DNA and each target
row is one checked strain's DNA from the start of its innermost left-flank
gene to the end of its innermost right-flank gene (flank genes included).

For every (locus, strain) this runs `blastn -task megablast` with the query
against the target (subject mode). A locus column is DNA present when HSPs
with identity >= --min_id cover >= --min_cov % of that exemplar gene's span
(overlapping HSPs are merged). A target shorter than 50 bp is DNA absent in
every column without an alignment. A strain without a target interval, or
whose genome FASTA is missing, is "unchecked".

Each strain's genome FASTA (data_dir/dna/<DNA column of the samplesheet>,
plain, .gz or .zst) is read once per batch.

Output: TSV locus_id, strain, col, status (present | absent | unchecked),
coverage (fraction of the gene span, 3 decimals).
"""
from __future__ import annotations

import argparse
import collections
import csv
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import open_maybe_compressed  # noqa: E402
from config_parser import parse_config  # noqa: E402
from island_locus import (  # noqa: E402
    DEFAULT_DNA_MIN_COV, DEFAULT_DNA_MIN_ID, DNA_MIN_TARGET_BP, gene_coverage,
)

CALL_COLUMNS = ["locus_id", "strain", "col", "status", "coverage"]


def read_targets(path: str) -> dict[str, dict]:
    """{locus_id: {"query": row, "targets": [row, ...]}} in file order.
    Rows are dicts with int start/end; the query row gains "genes" =
    [(col, start, end)]."""
    loci: dict[str, dict] = {}
    with open_maybe_compressed(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row["start"], row["end"] = int(row["start"]), int(row["end"])
            entry = loci.setdefault(row["locus_id"], {"query": None, "targets": []})
            if row["role"] == "query":
                genes = []
                for part in filter(None, row["genes"].split(";")):
                    col, span = part.split(":")
                    a, b = span.split("-")
                    genes.append((int(col), int(a), int(b)))
                row["genes"] = genes
                entry["query"] = row
            else:
                entry["targets"].append(row)
    return loci


def read_fasta_slices(path: str, wanted: dict[str, set[tuple[int, int]]]) -> dict:
    """{(contig, start, end): sequence} for 1-based closed intervals on the
    contigs named in `wanted`; the contig name is the header's first word.
    The file is read once, whole (a fungal genome is tens of MB)."""
    out = {}
    with open_maybe_compressed(path) as fh:
        text = fh.read()
    for record in text.split("\n>"):
        record = record.lstrip(">")
        head, _, body = record.partition("\n")
        name = head.split(None, 1)[0] if head.strip() else ""
        if name not in wanted:
            continue
        seq = body.replace("\n", "").replace("\r", "")
        for start, end in wanted[name]:
            out[(name, start, end)] = seq[start - 1:end]
    return out


def blastn_hsps(query_fa: Path, subject_fa: Path, blastn: str) -> list[tuple[int, int, float]]:
    """[(qstart, qend, pident)] of blastn -task megablast, query vs subject."""
    proc = subprocess.run([blastn, "-task", "megablast", "-query", str(query_fa),
                           "-subject", str(subject_fa), "-outfmt", "6 qstart qend pident"],
                          capture_output=True, text=True, check=True)
    hsps = []
    for line in proc.stdout.splitlines():
        q0, q1, pid = line.split("\t")
        hsps.append((int(q0), int(q1), float(pid)))
    return hsps


def check_batch(loci: dict[str, dict], genome_of: dict[str, str], blastn: str = "blastn",
                min_id: float = DEFAULT_DNA_MIN_ID, min_cov: float = DEFAULT_DNA_MIN_COV,
                cpus: int = 1) -> list[list]:
    """Call rows (CALL_COLUMNS order) for every (locus, target strain, query
    gene column). `genome_of` maps strain -> genome FASTA path."""
    wanted: dict[str, dict[str, set]] = collections.defaultdict(lambda: collections.defaultdict(set))
    for entry in loci.values():
        q = entry["query"]
        wanted[q["strain"]][q["contig"]].add((q["start"], q["end"]))
        for t in entry["targets"]:
            if t["contig"] != "-" and t["end"] - t["start"] + 1 >= DNA_MIN_TARGET_BP:
                wanted[t["strain"]][t["contig"]].add((t["start"], t["end"]))
    seqs: dict[tuple[str, str, int, int], str] = {}
    for strain, contigs in sorted(wanted.items()):
        path = genome_of.get(strain)
        if not path or not Path(path).is_file():
            print(f"WARNING: no genome FASTA for {strain}; its targets are unchecked",
                  file=sys.stderr)
            continue
        for (contig, start, end), seq in read_fasta_slices(path, contigs).items():
            seqs[(strain, contig, start, end)] = seq

    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        jobs = []
        for n, (locus_id, entry) in enumerate(loci.items()):
            q = entry["query"]
            qseq = seqs.get((q["strain"], q["contig"], q["start"], q["end"]))
            genes = [(col, a - q["start"] + 1, b - q["start"] + 1) for col, a, b in q["genes"]]
            qfa = Path(tmp, f"q{n}.fa")
            if qseq:
                qfa.write_text(f">query\n{qseq}\n")
            for m, t in enumerate(entry["targets"]):
                length = t["end"] - t["start"] + 1
                key = (t["strain"], t["contig"], t["start"], t["end"])
                if not qseq or t["contig"] == "-" or (length >= DNA_MIN_TARGET_BP and key not in seqs):
                    rows += [[locus_id, t["strain"], col, "unchecked", ""] for col, _, _ in genes]
                elif length < DNA_MIN_TARGET_BP:
                    rows += [[locus_id, t["strain"], col, "absent", "0.000"] for col, _, _ in genes]
                else:
                    tfa = Path(tmp, f"t{n}_{m}.fa")
                    tfa.write_text(f">target\n{seqs[key]}\n")
                    jobs.append((locus_id, t["strain"], genes, qfa, tfa))

        def run(job):
            locus_id, strain, genes, qfa, tfa = job
            hsps = blastn_hsps(qfa, tfa, blastn)
            out = []
            for col, g0, g1 in genes:
                cov = gene_coverage(hsps, (g0, g1), min_id)
                status = "present" if cov * 100 >= min_cov else "absent"
                out.append([locus_id, strain, col, status, f"{cov:.3f}"])
            return out

        with ThreadPoolExecutor(max_workers=max(1, cpus)) as pool:
            for out in pool.map(run, jobs):
                rows += out
    return rows


def genome_paths(config: str, data_dir: str) -> dict[str, str]:
    """{Short: data_dir/dna/<DNA>} from the samplesheet."""
    return {s.short: str(Path(data_dir) / "dna" / s.dna) for s in parse_config(config)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--targets", required=True, help="batch_NNN.tsv from pangenome_island_loci.py")
    ap.add_argument("--config", required=True, help="samplesheet CSV (Short, DNA columns)")
    ap.add_argument("--data_dir", required=True, help="study data_dir holding dna/")
    ap.add_argument("--min_id", type=float, default=DEFAULT_DNA_MIN_ID,
                    help="minimum HSP identity, percent (default 90)")
    ap.add_argument("--min_cov", type=float, default=DEFAULT_DNA_MIN_COV,
                    help="minimum covered share of an exemplar gene, percent (default 80)")
    ap.add_argument("--cpus", type=int, default=1, help="parallel blastn calls")
    ap.add_argument("--blastn", default="blastn")
    ap.add_argument("--output", required=True)
    args = ap.parse_args(argv)
    loci = read_targets(args.targets)
    rows = check_batch(loci, genome_paths(args.config, args.data_dir), args.blastn,
                       args.min_id, args.min_cov, args.cpus)
    with open(args.output, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(CALL_COLUMNS)
        w.writerows(rows)
    status = collections.Counter(r[3] for r in rows)
    print(f"pangenome_island_dna_check: {len(loci)} loci, "
          f"{sum(len(e['targets']) for e in loci.values())} strains checked; cells "
          f"{dict(sorted(status.items()))}; wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
