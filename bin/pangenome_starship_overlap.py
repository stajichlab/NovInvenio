#!/usr/bin/env python3
"""Overlap of Starship-annotated pangenome families with pangenome accessory islands.

Ported from NovInvenio_Investigations (NII) bin/starship_island_overlap.py (Jason
Stajich, NII commits acd0961 and 034f37b, 2026-10-17). The statistics, the output file
names, columns and metric names are unchanged, so NII's bin/sync_pangenome_report.py
(stage_starship_overlap) reads this script's output directory as it reads the NII one.
Changes from the NII script: inputs can be given as single files (the Nextflow process
does this) instead of only --run-dir; the element table follows lib/starship_tables.py
(column 'strain', alias 'Short'; label/confidence optional; a gene-level form); the
member ID separator is --id_sep; 0-byte stub files count as "not given"; a gene/loss
table is optional and its absence is reported in input_checks.tsv.

Inputs (pangenome.nf outputs; each file flag overrides the --run_dir lookup):
  --frequency_table     frequency_table.tsv
  --significant_islands significant_islands.tsv
  --cluster_tsv         cluster/tier1_cluster.tsv (rep <tab> member, member = Short|protein_id)
  --gene_positions      gene_positions.tsv.zst
  --family_positions    family_positions.tsv.zst
  --rescue_positions    rescue_positions.tsv (0-byte when the rescue pass was off)
  --starship_annotation element table (lib/starship_tables.py), coordinates in >= 1 strain
  --captain_tblout      HMMER --tblout of a captain-gene search (e.g. DUF3435) against the
                        run's proteins; target names are cluster member IDs. Hits with
                        full-sequence E-value <= --captain_evalue count.
  --gain_loss           optional per-family table with columns family, direction (e.g.
                        NII analysis/gain_loss/<run>/family_gain_loss.tsv.gz). Without it
                        the direction column of class_breakdown.tsv is 'NA'.

Definitions (all written to summary.tsv as well):

  family            a row of frequency_table.tsv (tier-1 cluster representative).
  island            a row of significant_islands.tsv (one distinct member-family set
                    per row, deduplicated across strains by the pipeline).
  island family     a family listed in member_families of >= 1 island.
  Starship family   a family with >= 1 member gene, in a strain that has called
                    elements, whose midpoint (start+end)//2 lies inside an element
                    (gene-level table: a member gene listed for the element).
  captain family    a family with >= 1 member hit by the captain HMM.
  island instance   for a strain with called elements: a maximal run of consecutive
                    non-core (bin not core/soft_core) entries of
                    family_positions.tsv on one contig whose member-family set equals
                    an island's member-family set (the pipeline's own rule). Its span
                    is taken from gene_positions (genes) and rescue_positions
                    (rescued hits, start only) in rank order.
  universes         'cooc_eligible': bin shell or cloud with strain_count >=
                    --min_carriers (pangenome_cooccurrence.py: shell+cloud, >= 5
                    representative carriers); 'noncore': any bin except core,
                    soft_core and '-'; for the Starship set also each universe
                    restricted to families with >= 1 gene in an annotated strain
                    ('..._in_annotated_strain'), because elements are only known there.

Outputs in --out_dir: summary.tsv, enrichment.tsv, class_breakdown.tsv,
per_starship.tsv, per_island.tsv.gz, annotated_families.tsv.gz,
coord_permutation.tsv, input_checks.tsv, PROVENANCE.yaml (inputs with sha256).
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from compressed_io import open_maybe_compressed  # noqa: E402
import starship_tables  # noqa: E402

from scipy.stats import fisher_exact, hypergeom  # noqa: E402

CORE_BINS = {"core", "soft_core"}


# ---------------------------------------------------------------- io helpers
def open_text(path: Path):
    return open_maybe_compressed(path)


def given(path: Path | None) -> Path | None:
    """None for an unset flag or a 0-byte stub file (EMPTY_EVALUES_STUB)."""
    if path is None:
        return None
    path = Path(path)
    if not path.exists():
        raise SystemExit(f"input not found: {path}")
    return path if path.stat().st_size > 0 else None


def first_existing(run_dir: Path | None, *names: str) -> Path | None:
    if run_dir is None:
        return None
    for n in names:
        p = run_dir / n
        if p.exists():
            return p
    return None


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------- loaders
def load_frequency(path: Path) -> dict[str, dict]:
    with open_text(path) as fh:
        return {r["family"]: r for r in csv.DictReader(fh, delimiter="\t")}


def load_islands(path: Path) -> list[dict]:
    out = []
    with open_text(path) as fh:
        for i, r in enumerate(csv.DictReader(fh, delimiter="\t")):
            fams = r["member_families"].split(",")
            out.append({"island_idx": i, "n_strains": int(r["n_strains"]),
                        "example_strain": r["example_strain"],
                        "island_size": int(r["island_size"]),
                        "families": fams, "family_set": frozenset(fams),
                        "classifications": r.get("classifications", "")})
    return out


def load_members(cluster_tsv: Path, strains: set[str], sep: str = "|") -> tuple[dict[str, str], dict[str, str]]:
    """member->rep for all members (needed for captain hits) and for the annotated strains."""
    all_m2r: dict[str, str] = {}
    with open_text(cluster_tsv) as fh:
        for line in fh:
            rep, mem = line.rstrip("\n").split("\t")[:2]
            all_m2r[mem] = rep
    strain_m2r = {m: r for m, r in all_m2r.items() if m.split(sep, 1)[0] in strains}
    return all_m2r, strain_m2r


def load_strain_genes(gene_positions: Path, strains: set[str]) -> list[tuple]:
    """[(Short, protein_id, contig, start, end)] for the annotated strains."""
    out = []
    with open_text(gene_positions) as fh:
        if next(fh, None) is None:
            return out
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[0] in strains:
                out.append((p[0], p[1], p[2], int(p[3]), int(p[4])))
    return out


def load_captain_families(tblout: Path, all_m2r: dict[str, str], evalue: float) -> tuple[set[str], int, int]:
    fams, n_hits, n_unmapped = set(), 0, 0
    with open_text(tblout) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if float(p[4]) > evalue:
                continue
            n_hits += 1
            rep = all_m2r.get(p[0])
            if rep is None:
                n_unmapped += 1
            else:
                fams.add(rep)
    return fams, n_hits, n_unmapped


def load_gain_loss(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    with open_text(path) as fh:
        return {r["family"]: r["direction"] for r in csv.DictReader(fh, delimiter="\t")}


def derive_gene_level_spans(els: list[dict], genes: list[tuple]) -> tuple[list[dict], list[str]]:
    """Gene-level elements -> span on the contig that holds most of their placed genes.

    Returns (elements with contig/start/end set, element_ids with no placed gene). An
    element whose genes lie on several contigs keeps the contig with the most genes.
    """
    pos = {(g[0], g[1]): g for g in genes}
    out, unplaced = [], []
    for e in els:
        if e.get("genes") is None:
            out.append(e)
            continue
        placed = [pos[(e["Short"], gid)] for gid in e["genes"] if (e["Short"], gid) in pos]
        if not placed:
            unplaced.append(e["element_id"])
            continue
        contig = Counter(g[2] for g in placed).most_common(1)[0][0]
        on = [g for g in placed if g[2] == contig]
        out.append({**e, "contig": contig, "start": min(g[3] for g in on), "end": max(g[4] for g in on)})
    return out, unplaced


# ---------------------------------------------------------------- core logic
def element_of(contig: str, mid: int, els_by_contig: dict[str, list[dict]]) -> str | None:
    for e in els_by_contig.get(contig, ()):
        if e["start"] <= mid <= e["end"]:
            return e["element_id"]
    return None


def starship_families(genes, strain_m2r, els, sep: str = "|") -> tuple[dict[str, set[str]], dict[str, list[tuple]]]:
    """element_id -> families; element_id -> genes [(Short, pid, contig, start, end, family)].

    Coordinate elements: a gene belongs when its midpoint lies inside the element.
    Gene-level elements (e['genes'] is a list): a gene belongs when it is listed.
    """
    by_contig: dict[tuple, list[dict]] = defaultdict(list)
    listed: dict[tuple, list[str]] = defaultdict(list)
    for e in els:
        if e.get("genes") is not None:
            for gid in e["genes"]:
                listed[(e["Short"], gid)].append(e["element_id"])
        else:
            by_contig[(e["Short"], e["contig"])].append(e)
    fams: dict[str, set[str]] = {e["element_id"]: set() for e in els}
    egenes: dict[str, list[tuple]] = {e["element_id"]: [] for e in els}
    for short, pid, contig, s, t in genes:
        hits = list(listed.get((short, pid), ()))
        cand = by_contig.get((short, contig))
        if cand:
            mid = (s + t) // 2
            hits += [e["element_id"] for e in cand if e["start"] <= mid <= e["end"]]
        for eid in hits:
            fam = strain_m2r.get(f"{short}{sep}{pid}")
            egenes[eid].append((short, pid, contig, s, t, fam))
            if fam is not None:
                fams[eid].add(fam)
    return fams, egenes


def strain_order(family_positions: Path, rescue_positions: Path | None, genes, strain_m2r, strains):
    """Short -> list of (contig, family, start, end) in the run's rank order.

    family_positions gives (family, contig, rank); coordinates come from the
    gene and rescue position tables sorted by (contig, start), which is the
    order the pipeline ranks them in (ties may swap; the coordinate is the same).
    """
    fp: dict[str, list[tuple]] = defaultdict(list)
    with open_text(family_positions) as fh:
        next(fh, None)
        for line in fh:
            s, fam, contig, rank = line.rstrip("\n").split("\t")[:4]
            if s in strains:
                fp[s].append((int(rank), contig, fam))
    coords: dict[str, list[tuple]] = defaultdict(list)
    for short, pid, contig, a, b in genes:
        coords[short].append((contig, a, b))
    if rescue_positions is not None:
        with open_text(rescue_positions) as fh:
            next(fh, None)
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if p[0] in strains:
                    coords[p[0]].append((p[2], int(p[3]), int(p[3])))
    out = {}
    for s in strains:
        rows = sorted(fp.get(s, []))
        cs = sorted(coords.get(s, []))
        if len(rows) != len(cs):
            raise SystemExit(f"{s}: family_positions has {len(rows)} rows but genes+rescues give {len(cs)}")
        out[s] = [(c, fam, cs[i][1], cs[i][2]) for i, (rank, c, fam) in enumerate(rows)]
    return out


def island_instances(order, freq, island_by_set) -> list[dict]:
    """Maximal non-core runs whose member set equals an island's member set."""
    inst = []
    for short, rows in order.items():
        run: list[tuple] = []
        prev = None

        def flush():
            if len(run) >= 2:
                key = frozenset(r[1] for r in run)
                idx = island_by_set.get(key)
                if idx is not None:
                    inst.append({"Short": short, "contig": run[0][0], "start": min(r[2] for r in run),
                                 "end": max(r[3] for r in run), "island_idx": idx, "n_genes": len(run)})

        for row in rows:
            contig, fam = row[0], row[1]
            core = freq.get(fam, {}).get("bin", "core") in CORE_BINS or fam not in freq
            if contig != prev and run:
                flush()
                run = []
            if not core:
                run.append(row)
            elif run:
                flush()
                run = []
            prev = contig
        if run:
            flush()
    return inst


def hyper_row(name, universe_name, universe, annotated, island_fams):
    N = len(universe)
    K = len(annotated & universe)
    n = len(island_fams & universe)
    k = len(annotated & island_fams & universe)
    exp = K * n / N if N else float("nan")
    p = float(hypergeom.sf(k - 1, N, K, n)) if N else float("nan")
    table = [[k, K - k], [n - k, N - K - n + k]]
    odds = float(fisher_exact(table, alternative="greater")[0]) if N else float("nan")
    return {"annotation": name, "universe": universe_name, "N_universe": N, "K_annotated": K,
            "n_island": n, "k_both": k, "expected_k": round(exp, 2),
            "fold": round(k / exp, 3) if exp else "", "odds_ratio": round(odds, 3),
            "p_hypergeom_greater": f"{p:.3g}",
            "frac_annotated_in_island": round(k / K, 4) if K else "",
            "frac_universe_in_island": round(n / N, 4) if N else ""}


def coord_permutation(genes_ann, els, island_fams, strain_m2r, freq, n_perm, seed, sep: str = "|"):
    """Shift each element to a random position on its own contig (same length).
    Statistic: annotated-strain genes in island families whose midpoint lies in an element.
    Contig length = the largest gene end on that contig (no genome FASTA is read)."""
    contig_len: dict[tuple, int] = defaultdict(int)
    for short, pid, contig, a, b in genes_ann:
        contig_len[(short, contig)] = max(contig_len[(short, contig)], b)
    pts: dict[tuple, list[tuple[int, bool, bool]]] = defaultdict(list)
    for short, pid, contig, a, b in genes_ann:
        fam = strain_m2r.get(f"{short}{sep}{pid}")
        noncore = fam is not None and freq.get(fam, {}).get("bin", "core") not in CORE_BINS | {"-"}
        pts[(short, contig)].append(((a + b) // 2, fam in island_fams, noncore))

    def stat(placed):
        isl = nonc = 0
        for (key, s, t) in placed:
            for mid, is_isl, is_nc in pts.get(key, ()):
                if s <= mid <= t:
                    isl += is_isl
                    nonc += is_nc
        return isl, nonc

    obs = stat([((e["Short"], e["contig"]), e["start"], e["end"]) for e in els])
    rng = random.Random(seed)
    null = []
    for _ in range(n_perm):
        placed = []
        for e in els:
            L = e["end"] - e["start"]
            clen = contig_len[(e["Short"], e["contig"])]
            s = rng.randint(1, max(1, clen - L))
            placed.append(((e["Short"], e["contig"]), s, s + L))
        null.append(stat(placed)[0])
    null.sort()
    ge = sum(1 for x in null if x >= obs[0])
    return {"statistic": "annotated-strain genes in island families inside elements",
            "observed": obs[0], "observed_noncore_genes_inside": obs[1],
            "n_perm": n_perm, "seed": seed,
            "null_mean": round(sum(null) / len(null), 2) if null else "",
            "null_q95": null[int(0.95 * (len(null) - 1))] if null else "",
            "null_max": null[-1] if null else "",
            "p_perm_greater": round((ge + 1) / (n_perm + 1), 5)}


# ---------------------------------------------------------------- writers
def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as fh:
        if not rows:
            fh.write("\t".join(fields or []) + "\n")
            return
        w = csv.DictWriter(fh, fieldnames=fields or list(rows[0]), delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    # Underscore flags (pipeline convention) with the NII dash spellings as aliases.
    ap.add_argument("--run_dir", "--run-dir", type=Path, help="pangenome.nf output dir (.../pangenome)")
    ap.add_argument("--frequency_table", type=Path)
    ap.add_argument("--significant_islands", type=Path)
    ap.add_argument("--cluster_tsv", type=Path)
    ap.add_argument("--gene_positions", type=Path)
    ap.add_argument("--family_positions", type=Path)
    ap.add_argument("--rescue_positions", type=Path)
    ap.add_argument("--samplesheet", type=Path,
                    help="optional: check that element strains are Short IDs of these rows")
    ap.add_argument("--starship_annotation", "--starship-annotation", type=Path,
                    help="element table (lib/starship_tables.py)")
    ap.add_argument("--captain_tblout", "--captain-tblout", type=Path,
                    help="HMMER tblout of a captain HMM vs the run's proteins")
    ap.add_argument("--captain_evalue", "--captain-evalue", type=float, default=1e-5)
    ap.add_argument("--gain_loss", "--gain-loss", type=Path, help="per-family table (columns family, direction)")
    ap.add_argument("--extra_gene_positions", "--extra-gene-positions", type=Path,
                    help="TSV Short, protein_id, contig, start, end for annotated-strain members missing from "
                         "gene_positions (used for Starship families, gene metrics and the permutation only)")
    ap.add_argument("--id_sep", default="|", help="Short/protein_id separator of cluster member IDs")
    ap.add_argument("--out_dir", "--out-dir", required=True, type=Path)
    ap.add_argument("--run_label", "--run-label", default="")
    ap.add_argument("--min_carriers", "--min-carriers", type=int, default=5)
    ap.add_argument("--n_perm", "--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261017)
    return ap


def main(argv=None) -> None:
    ap = build_parser()
    a = ap.parse_args(argv)
    a.starship_annotation = given(a.starship_annotation)
    a.captain_tblout = given(a.captain_tblout)
    a.gain_loss = given(a.gain_loss)
    a.extra_gene_positions = given(a.extra_gene_positions)
    if not a.starship_annotation and not a.captain_tblout:
        ap.error("give --starship_annotation and/or --captain_tblout (a 0-byte file counts as not given)")
    a.out_dir.mkdir(parents=True, exist_ok=True)
    rd = a.run_dir
    sep = a.id_sep
    paths = {
        "frequency_table": a.frequency_table or first_existing(rd, "frequency_table.tsv"),
        "significant_islands": a.significant_islands or first_existing(rd, "significant_islands.tsv"),
        "tier1_cluster": a.cluster_tsv or first_existing(rd, "cluster/tier1_cluster.tsv",
                                                          "cluster/tier1_cluster.tsv.zst"),
        "gene_positions": a.gene_positions or first_existing(rd, "gene_positions.tsv.zst", "gene_positions.tsv"),
        "family_positions": a.family_positions or first_existing(rd, "family_positions.tsv.zst",
                                                                 "family_positions.tsv"),
        "rescue_positions": (given(a.rescue_positions) if a.rescue_positions is not None
                             else first_existing(rd, "rescue_positions.tsv", "rescue_positions.tsv.zst")),
    }
    for k in ("frequency_table", "significant_islands", "tier1_cluster"):
        if paths[k] is None:
            raise SystemExit(f"missing {k} (give --{k if k != 'tier1_cluster' else 'cluster_tsv'} or --run_dir)")

    checks: list[dict] = []

    def C(check, value, note=""):
        checks.append({"check": check, "value": value, "note": note})

    freq = load_frequency(paths["frequency_table"])
    islands = load_islands(paths["significant_islands"])
    island_by_set = {}
    for isl in islands:
        island_by_set.setdefault(isl["family_set"], isl["island_idx"])
    island_fams = set().union(*(i["family_set"] for i in islands)) & set(freq) if islands else set()
    island_fams_multi = (set().union(*(i["family_set"] for i in islands if i["n_strains"] >= 2)) & set(freq)
                         if islands else set())
    gl = load_gain_loss(a.gain_loss)
    C("gain_loss", "given" if a.gain_loss else "not_computed",
      "per-family direction table" if a.gain_loss else
      "no per-family gain/loss table; class_breakdown.tsv direction = NA")

    known = starship_tables.samplesheet_strains(a.samplesheet) if a.samplesheet else None
    try:
        els = starship_tables.load_elements(a.starship_annotation, known) if a.starship_annotation else []
    except starship_tables.TableError as exc:
        raise SystemExit(f"ERROR: --starship_annotation: {exc}")
    strains = {e["Short"] for e in els}
    all_m2r, strain_m2r = load_members(paths["tier1_cluster"], strains, sep)

    sets: dict[str, set[str]] = {}
    summary: list[dict] = []

    def S(metric, value, note=""):
        summary.append({"run": a.run_label, "metric": metric, "value": value, "definition": note})

    S("n_families", len(freq), "rows of frequency_table.tsv")
    S("n_islands", len(islands), "rows of significant_islands.tsv")
    S("n_islands_multi_strain", sum(i["n_strains"] >= 2 for i in islands), "islands with n_strains >= 2")
    S("n_island_families", len(island_fams), "families in >= 1 island (and in frequency_table)")
    S("n_island_families_multi_strain", len(island_fams_multi), "families in >= 1 island with n_strains >= 2")

    fam_el: dict[str, set[str]] = defaultdict(set)
    egenes = {}
    genes_ann = []
    genes_run = []
    if els:
        if paths["gene_positions"] is None:
            raise SystemExit("gene_positions table needed for --starship_annotation")
        genes_run = load_strain_genes(paths["gene_positions"], strains)
        genes_ann = list(genes_run)
        extra = []
        if a.extra_gene_positions:
            have = {(g[0], g[1]) for g in genes_run}
            extra = [g for g in load_strain_genes(a.extra_gene_positions, strains) if (g[0], g[1]) not in have]
            genes_ann += extra
        C("element_rows", len(els), "elements in --starship_annotation")
        C("element_strains", ",".join(sorted(strains)), "")
        if any(e.get("genes") is not None for e in els):
            els, unplaced = derive_gene_level_spans(els, genes_ann)
            C("elements_gene_level_unplaced", len(unplaced),
              "gene-level elements with no gene in gene_positions; dropped: " + ",".join(unplaced))
            strains = {e["Short"] for e in els}
        contigs = {(g[0], g[2]) for g in genes_ann}
        unmatched = [e["element_id"] for e in els if (e["Short"], e["contig"]) not in contigs]
        C("elements_contig_not_in_gene_positions", len(unmatched),
          "element (strain, contig) with no gene in gene_positions: " + ",".join(unmatched))
        if unmatched:
            print(f"WARNING: {len(unmatched)} element(s) on a contig with no gene in gene_positions "
                  f"(contig names must match the strain's GFF3): {','.join(unmatched[:10])}", file=sys.stderr)
        if a.extra_gene_positions:
            S("n_extra_gene_positions_used", len(extra),
              "annotated-strain genes added from --extra-gene-positions (not used for island instances)")
        efams, egenes = starship_families(genes_ann, strain_m2r, els, sep)
        for eid, fs in efams.items():
            for f in fs:
                fam_el[f].add(eid)
        sets["starship_all"] = set(fam_el) & set(freq)
        high = {e["element_id"] for e in els if e["confidence"] == "high"}
        sets["starship_high"] = {f for f, es in fam_el.items() if es & high} & set(freq)
        S("annotated_strains", ",".join(sorted(strains)), "strains with called elements")
        S("n_elements", len(els), "rows of --starship-annotation")
        S("n_elements_high", len(high), "confidence == high")
        S("n_genes_in_elements", sum(len(v) for v in egenes.values()), "annotated-strain genes with midpoint in an element")
    if a.captain_tblout:
        cap, n_hits, n_unmapped = load_captain_families(a.captain_tblout, all_m2r, a.captain_evalue)
        sets["captain"] = cap & set(freq)
        S("captain_hits", n_hits, f"tblout targets with full-sequence E <= {a.captain_evalue}")
        S("captain_hits_unmapped", n_unmapped, "hits whose target is not a cluster member")
    for name, fs in sets.items():
        S(f"n_{name}_families", len(fs), "families in the annotation set")
        S(f"n_{name}_families_in_island", len(fs & island_fams), "of those, in >= 1 island")
        S(f"n_{name}_families_in_multi_strain_island", len(fs & island_fams_multi), "of those, in >= 1 island with n_strains >= 2")
        n_isl = sum(1 for i in islands if i["family_set"] & fs)
        S(f"n_islands_with_{name}_family", n_isl, "islands with >= 1 member family in the set")
        S(f"frac_islands_with_{name}_family", round(n_isl / len(islands), 4) if islands else "", "")
        S(f"frac_island_families_{name}", round(len(fs & island_fams) / len(island_fams), 4) if island_fams else "", "")

    # ---- enrichment
    cooc = {f for f, r in freq.items() if r["bin"] in ("shell", "cloud") and int(r["strain_count"]) >= a.min_carriers}
    noncore = {f for f, r in freq.items() if r["bin"] not in CORE_BINS | {"-"}}
    in_strain = set(strain_m2r.values())
    enrich = []
    for name, fs in sets.items():
        for uname, U in (("cooc_eligible", cooc), ("noncore", noncore)):
            enrich.append(hyper_row(name, uname, U, fs, island_fams))
            if name.startswith("starship"):
                enrich.append(hyper_row(name, uname + "_in_annotated_strain", U & in_strain, fs, island_fams))
    write_tsv(a.out_dir / "enrichment.tsv", enrich)

    # ---- class / direction breakdown
    rows = []
    for name, fs in sets.items():
        cnt = Counter()
        for f in fs:
            cnt[(freq[f]["bin"], gl.get(f, "NA"), f in island_fams)] += 1
        for (b, d, inn), n in sorted(cnt.items()):
            rows.append({"annotation": name, "bin": b, "direction": d, "in_island": "yes" if inn else "no", "n_families": n})
    write_tsv(a.out_dir / "class_breakdown.tsv", rows, ["annotation", "bin", "direction", "in_island", "n_families"])

    # ---- coordinate level (annotated strains)
    inst = []
    if els and paths["family_positions"] is not None:
        order = strain_order(paths["family_positions"], paths["rescue_positions"], genes_run, strain_m2r, strains)
        inst = island_instances(order, freq, island_by_set)
        S("n_island_instances_in_annotated_strains", len(inst), "non-core runs in annotated strains equal to an island member set")
    per_ship = []
    island_hits: dict[int, set[str]] = defaultdict(set)
    island_coord_hits: dict[int, set[str]] = defaultdict(set)
    for isl in islands:
        for f in isl["family_set"]:
            for eid in fam_el.get(f, ()):
                island_hits[isl["island_idx"]].add(eid)
    for e in els:
        eid = e["element_id"]
        fs = {g[5] for g in egenes[eid] if g[5] is not None} & set(freq)
        ov = [i for i in inst if i["Short"] == e["Short"] and i["contig"] == e["contig"]
              and i["start"] <= e["end"] and i["end"] >= e["start"]]
        for i in ov:
            island_coord_hits[i["island_idx"]].add(eid)
        genes_in_inst = sum(1 for g in egenes[eid]
                            if any(i["start"] <= (g[3] + g[4]) // 2 <= i["end"] and i["contig"] == g[2] for i in ov))
        bins = Counter(freq[f]["bin"] for f in fs)
        per_ship.append({
            "element_id": eid, "label": e["label"], "confidence": e["confidence"], "Short": e["Short"],
            "contig": e["contig"], "start": e["start"], "end": e["end"], "length_bp": e["end"] - e["start"] + 1,
            "n_genes": len(egenes[eid]), "n_families": len(fs),
            "n_noncore_families": sum(v for b, v in bins.items() if b not in CORE_BINS | {"-"}),
            "n_families_in_island": len(fs & island_fams),
            "frac_families_in_island": round(len(fs & island_fams) / len(fs), 4) if fs else "",
            "n_islands_hit_family_level": sum(1 for v in island_hits.values() if eid in v),
            "n_island_instances_overlapping": len(ov),
            "n_genes_inside_island_instances": genes_in_inst,
            "frac_genes_inside_island_instances": round(genes_in_inst / len(egenes[eid]), 4) if egenes[eid] else "",
            "bins": ";".join(f"{b}:{v}" for b, v in sorted(bins.items())),
        })
    if els:
        write_tsv(a.out_dir / "per_starship.tsv", per_ship)
        S("n_islands_hit_family_level", len(island_hits), "islands with >= 1 member family that is a Starship family")
        S("n_islands_overlapping_element_in_annotated_strain", len(island_coord_hits),
          "islands with an instance in an annotated strain overlapping an element by >= 1 bp")
        hi = {e["element_id"] for e in els if e["confidence"] == "high"}
        S("n_high_elements_hit_family_level", len(hi & set().union(*island_hits.values())) if island_hits else 0, "")
        S("n_high_elements_overlapped_by_instance", len(hi & set().union(*island_coord_hits.values())) if island_coord_hits else 0, "")
        # gene metric of the NII island genome-count sweep (score_starships.py)
        fam_genes: dict[str, list] = defaultdict(list)
        for short, pid, contig, s, t in genes_ann:
            fam = strain_m2r.get(f"{short}{sep}{pid}")
            if fam:
                fam_genes[fam].append((short, contig, (s + t) // 2))
        by_c = defaultdict(list)
        for e in els:
            by_c[(e["Short"], e["contig"])].append(e)

        def inside(g):
            return any(e["start"] <= g[2] <= e["end"] for e in by_c.get((g[0], g[1]), ()))

        placed = placed_in = 0
        for isl in islands:
            for f in isl["families"]:
                for g in fam_genes.get(f, ()):
                    placed += 1
                    placed_in += inside(g)
        nc = [g for f, gs in fam_genes.items() if freq.get(f, {}).get("bin") not in (CORE_BINS | {"-", None}) for g in gs]
        S("island_gene_frac_in_starship_sweep_metric", round(placed_in / placed, 4) if placed else "",
          f"sum over islands of annotated-strain genes of member families inside elements / all such genes ({placed_in}/{placed})")
        S("noncore_gene_frac_in_starship_sweep_metric", round(sum(map(inside, nc)) / len(nc), 4) if nc else "",
          f"annotated-strain genes of non-core families inside elements ({sum(map(inside, nc))}/{len(nc)})")
        perm = coord_permutation(genes_ann, els, island_fams, strain_m2r, freq, a.n_perm, a.seed, sep)
        write_tsv(a.out_dir / "coord_permutation.tsv", [perm])

    # ---- per island
    cap = sets.get("captain", set())
    ss = sets.get("starship_all", set())
    pi = []
    inst_by_isl = defaultdict(list)
    for i in inst:
        inst_by_isl[i["island_idx"]].append(f"{i['Short']}:{i['contig']}:{i['start']}-{i['end']}")
    for isl in islands:
        fs = isl["family_set"]
        pi.append({"island_idx": isl["island_idx"], "n_strains": isl["n_strains"], "example_strain": isl["example_strain"],
                   "island_size": isl["island_size"], "n_member_families": len(fs),
                   "n_starship_families": len(fs & ss), "starship_elements_family_level": ",".join(sorted(island_hits.get(isl["island_idx"], ()))),
                   "starship_elements_coordinate": ",".join(sorted(island_coord_hits.get(isl["island_idx"], ()))),
                   "n_captain_families": len(fs & cap),
                   "instances_in_annotated_strains": ";".join(inst_by_isl.get(isl["island_idx"], ()))})
    write_tsv(a.out_dir / "per_island.tsv.gz", pi)

    af = []
    for f in sorted(set().union(*sets.values()) if sets else ()):
        af.append({"family": f, "bin": freq[f]["bin"], "strain_count": freq[f]["strain_count"],
                   "direction": gl.get(f, "NA"), "starship_elements": ",".join(sorted(fam_el.get(f, ()))),
                   "captain": "yes" if f in cap else "no", "in_island": "yes" if f in island_fams else "no",
                   "in_multi_strain_island": "yes" if f in island_fams_multi else "no"})
    write_tsv(a.out_dir / "annotated_families.tsv.gz", af)
    write_tsv(a.out_dir / "summary.tsv", summary)
    write_tsv(a.out_dir / "input_checks.tsv", checks, ["check", "value", "note"])

    # ---- provenance
    inputs = {k: v for k, v in paths.items() if v is not None}
    for k in ("starship_annotation", "captain_tblout", "gain_loss", "extra_gene_positions"):
        if getattr(a, k):
            inputs[k] = getattr(a, k)
    with open(a.out_dir / "PROVENANCE.yaml", "w") as fh:
        fh.write("derived_by: bin/pangenome_starship_overlap.py (port of NII bin/starship_island_overlap.py)\n")
        fh.write(f"date: {dt.date.today().isoformat()}\n")
        fh.write(f"command: {' '.join(sys.argv)}\n")
        fh.write(f"parameters: {{min_carriers: {a.min_carriers}, captain_evalue: {a.captain_evalue}, "
                 f"n_perm: {a.n_perm}, seed: {a.seed}, id_sep: '{sep}'}}\n")
        fh.write("inputs:\n")
        for k, p in inputs.items():
            fh.write(f"  {k}:\n    path: {Path(p).resolve()}\n    sha256: {sha256(Path(p))}\n")
    print(f"wrote {a.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
