#!/usr/bin/env python3
"""Secondary-metabolite (SM) clusters vs Starship elements and pangenome islands.

Generic, table-driven port of NovInvenio_Investigations (NII)
studies/fungi/Afumigatus_pangenome/analysis/starship_island_overlap/full_v072_antismash/
scripts/antismash_sm_starship.py (Jason Stajich, NII commit 034f37b, 2026-10-17). The NII
script parsed one antiSMASH JSON for one strain (Af293); this one reads a precomputed SM
cluster table for any number of strains (schema: lib/starship_tables.py) and maps cluster
genes to the run through gene_positions. Output file names, columns and metric names
follow the NII script where the meaning is the same, so NII's
bin/sync_pangenome_report.py reads this output directory as an SM dir
(--starship_sm_dir). Dropped: the antiSMASH JSON parser, the locus-tag mapping, the
PKS/NRPS proxy comparison (proxy_* metrics) and the Af293 constant.

Inputs
  --sm_clusters          SM cluster table (lib/starship_tables.py).
  --starship_annotation  optional element table; without it no Starship overlap is
                         computed (overlap columns are 'NA', no sm_starship_overlap.tsv,
                         no sm_permutation.tsv).
  --frequency_table, --significant_islands, --cluster_tsv, --gene_positions
                         pangenome.nf outputs (or --run_dir).
  --extra_gene_positions optional TSV Short, protein_id, contig, start, end for genes
                         missing from gene_positions.
  --captain_tblout       optional captain-gene hmmsearch tblout (targets Short|protein_id).
  --gain_loss            optional per-family table (family, direction).
  --contig_lengths       optional TSV strain, contig, length for the permutation; default
                         contig length = the largest gene end or cluster end on the contig.

Definitions
  SM-annotated strain   a strain with >= 1 row in the SM table. Strains without rows are
                        treated as not annotated (not as "no SM cluster").
  strain genes          the strain's rows of gene_positions (+ --extra_gene_positions).
  cluster genes         the table's 'genes' list when given (an ID that is not a strain
                        gene is kept, counted in n_genes and as unplaced, and is never
                        inside an element); else the strain genes with start >= cluster
                        start and end <= cluster end, in start order.
  core gene             a gene listed in 'core_genes' (empty when the column is absent).
  gene inside element   gene midpoint (start+end)//2 inside an element (1-based inclusive),
                        the rule of bin/pangenome_starship_overlap.py.
  overlap rule A        >= 1 cluster gene inside an element.
  overlap rule B        >= 50% of the cluster genes inside elements.
  Overlap is tested only for clusters in strains that also have elements; other
  clusters get 'NA'.
  captain gene          a cluster gene whose member has a captain hit, E <= --captain_evalue.

Outputs (--out_dir): sm_clusters.tsv.gz, sm_by_type.tsv, sm_starship_overlap.tsv,
sm_permutation.tsv, sm_core_families.tsv.gz, sm_family_class_direction.tsv, summary.tsv
(metric, value, source), input_checks.tsv, PROVENANCE.yaml.
"""
from __future__ import annotations

import argparse
import datetime as dt
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
import pangenome_starship_overlap as sio  # noqa: E402
import starship_tables  # noqa: E402

from scipy.stats import hypergeom  # noqa: E402

NA = "NA"


def hyper(N, K, n, k):
    exp = K * n / N if N else float("nan")
    return {"N": N, "K": K, "n": n, "k": k, "expected": round(exp, 2),
            "fold": round(k / exp, 2) if exp else "",
            "p_greater": f"{hypergeom.sf(k - 1, N, K, n):.3g}" if N else ""}


def load_contig_lengths(path: Path | None) -> dict[tuple, int]:
    if path is None:
        return {}
    _, rows = starship_tables.read_table(path)
    return {(r["strain"], r["contig"]): int(r["length"]) for r in rows}


def load_captain_members(tblout: Path | None, evalue: float) -> set[str]:
    out: set[str] = set()
    if tblout is None:
        return out
    with sio.open_text(tblout) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if float(p[4]) <= evalue:
                out.add(p[0])
    return out


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sm_clusters", required=True, type=Path)
    ap.add_argument("--starship_annotation", type=Path)
    ap.add_argument("--run_dir", type=Path)
    ap.add_argument("--frequency_table", type=Path)
    ap.add_argument("--significant_islands", type=Path)
    ap.add_argument("--cluster_tsv", type=Path)
    ap.add_argument("--gene_positions", type=Path)
    ap.add_argument("--extra_gene_positions", type=Path)
    ap.add_argument("--captain_tblout", type=Path)
    ap.add_argument("--captain_evalue", type=float, default=1e-5)
    ap.add_argument("--gain_loss", type=Path)
    ap.add_argument("--contig_lengths", type=Path)
    ap.add_argument("--samplesheet", type=Path)
    ap.add_argument("--id_sep", default="|")
    ap.add_argument("--n_perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261017)
    ap.add_argument("--out_dir", required=True, type=Path)
    return ap


def main(argv=None) -> None:
    a = build_parser().parse_args(argv)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    rd = a.run_dir
    sep = a.id_sep
    for k in ("starship_annotation", "extra_gene_positions", "captain_tblout", "gain_loss", "contig_lengths"):
        setattr(a, k, sio.given(getattr(a, k)))
    paths = {
        "frequency_table": a.frequency_table or sio.first_existing(rd, "frequency_table.tsv"),
        "significant_islands": a.significant_islands or sio.first_existing(rd, "significant_islands.tsv"),
        "tier1_cluster": a.cluster_tsv or sio.first_existing(rd, "cluster/tier1_cluster.tsv",
                                                              "cluster/tier1_cluster.tsv.zst"),
        "gene_positions": a.gene_positions or sio.first_existing(rd, "gene_positions.tsv.zst", "gene_positions.tsv"),
    }
    for k, p in paths.items():
        if p is None:
            raise SystemExit(f"missing {k} (give the file flag or --run_dir)")

    known = starship_tables.samplesheet_strains(a.samplesheet) if a.samplesheet else None
    try:
        clusters = starship_tables.load_sm_clusters(a.sm_clusters, known)
        els = starship_tables.load_elements(a.starship_annotation, known) if a.starship_annotation else []
    except starship_tables.TableError as exc:
        raise SystemExit(f"ERROR: {exc}")

    checks: list[dict] = []

    def C(check, value, note=""):
        checks.append({"check": check, "value": value, "note": note})

    # ---------------- pangenome tables
    freq = sio.load_frequency(paths["frequency_table"])
    islands = sio.load_islands(paths["significant_islands"])
    island_fams = set().union(*(i["family_set"] for i in islands)) & set(freq) if islands else set()
    gl = sio.load_gain_loss(a.gain_loss)
    C("gain_loss", "given" if a.gain_loss else "not_computed",
      "per-family direction table" if a.gain_loss else
      "no per-family gain/loss table; sm_family_class_direction.tsv direction = NA")
    sm_strains = {c["Short"] for c in clusters}
    all_m2r, strain_m2r = sio.load_members(paths["tier1_cluster"], sm_strains, sep)
    nonc_f = {f for f, r in freq.items() if r["bin"] not in sio.CORE_BINS | {"-"}}

    # ---------------- strain genes
    genes = sio.load_strain_genes(paths["gene_positions"], sm_strains)
    if a.extra_gene_positions:
        have = {(g[0], g[1]) for g in genes}
        extra = [g for g in sio.load_strain_genes(a.extra_gene_positions, sm_strains) if (g[0], g[1]) not in have]
        genes += extra
        C("extra_gene_positions_used", len(extra), "SM-strain genes added from --extra_gene_positions")
    if els and any(e.get("genes") is not None for e in els):
        els, unplaced = sio.derive_gene_level_spans(els, genes)
        C("elements_gene_level_unplaced", len(unplaced), ",".join(unplaced))
    gene = {}      # (Short, pid) -> dict
    by_sc = defaultdict(list)
    for short, pid, contig, s, t in genes:
        g = {"Short": short, "pid": pid, "contig": contig, "start": s, "end": t, "mid": (s + t) // 2}
        gene[(short, pid)] = g
        by_sc[(short, contig)].append(g)
    for lst in by_sc.values():
        lst.sort(key=lambda g: (g["start"], g["end"], g["pid"]))
    contigs = set(by_sc)
    unmatched = [c["cluster_id"] for c in clusters if (c["Short"], c["contig"]) not in contigs]
    C("sm_clusters_rows", len(clusters), "rows of --sm_clusters")
    C("sm_clusters_contig_not_in_gene_positions", len(unmatched), ",".join(unmatched[:50]))
    if unmatched:
        print(f"WARNING: {len(unmatched)} SM cluster(s) on a contig with no gene in gene_positions "
              f"(contig names must match the strain's GFF3): {','.join(unmatched[:10])}", file=sys.stderr)

    cap_members = load_captain_members(a.captain_tblout, a.captain_evalue)
    C("captain_tblout", "given" if a.captain_tblout else "not_given", "")

    # ---------------- elements
    el_strains = {e["Short"] for e in els}
    els_by_c = defaultdict(list)
    for e in els:
        els_by_c[(e["Short"], e["contig"])].append(e)
    if els:
        no_sm = sorted(el_strains - sm_strains)
        C("element_strains_without_sm_rows", len(no_sm), ",".join(no_sm))

    def el_of(key, mid, placed=None):
        for e in (placed if placed is not None else els_by_c).get(key, ()):
            if e["start"] <= mid <= e["end"]:
                return e["element_id"]
        return None

    for g in gene.values():
        g["element"] = el_of((g["Short"], g["contig"]), g["mid"]) or ""
        g["family"] = strain_m2r.get(f"{g['Short']}{sep}{g['pid']}", "")
        g["captain"] = f"{g['Short']}{sep}{g['pid']}" in cap_members

    # ---------------- cluster genes
    n_unplaced = n_core_not_in_genes = 0
    for c in clusters:
        if c["genes"]:
            c["gene_ids"] = list(c["genes"])
        else:
            c["gene_ids"] = [g["pid"] for g in by_sc.get((c["Short"], c["contig"]), ())
                             if c["start"] <= g["start"] and g["end"] <= c["end"]]
        c["placed"] = [gene[(c["Short"], x)] for x in c["gene_ids"] if (c["Short"], x) in gene]
        n_unplaced += len(c["gene_ids"]) - len(c["placed"])
        n_core_not_in_genes += sum(1 for x in c["core_genes"] if x not in c["gene_ids"])
        c["testable"] = c["Short"] in el_strains
    C("sm_cluster_genes_unplaced", n_unplaced, "listed cluster genes not in gene_positions (never inside an element)")
    C("sm_core_genes_not_in_cluster_genes", n_core_not_in_genes, "core_genes IDs not in the cluster's gene list")

    def gfield(c, x, key):
        g = gene.get((c["Short"], x))
        return g[key] if g else ""

    # ---------------- per cluster
    rows = []
    clusters.sort(key=lambda c: (c["Short"], c["contig"], c["start"], c["cluster_id"]))
    for c in clusters:
        g = c["gene_ids"]
        ins = [x for x in g if gfield(c, x, "element")]
        frac = len(ins) / len(g) if g else 0.0
        c.update(n_in=len(ins), frac=frac, ruleA=c["testable"] and len(ins) >= 1,
                 ruleB=c["testable"] and frac >= 0.5,
                 elements=sorted({gfield(c, x, "element") for x in ins}),
                 captain=[x for x in g if gfield(c, x, "captain")])

        def yn(flag):
            return ("yes" if flag else "no") if c["testable"] else NA
        rows.append({"cluster_id": c["cluster_id"], "strain": c["Short"], "contig": c["contig"],
                     "start": c["start"], "end": c["end"],
                     "length_bp": c["end"] - c["start"] + 1, "type": c["type"], "contig_edge": c["contig_edge"],
                     "n_genes": len(g), "n_core_genes": len(c["core_genes"]),
                     "core_genes": ",".join(c["core_genes"]),
                     "core_members": ",".join(c["core_genes"]),
                     "core_families": ",".join(gfield(c, x, "family") for x in c["core_genes"]),
                     "n_genes_in_starship": len(ins) if c["testable"] else NA,
                     "frac_genes_in_starship": round(frac, 4) if c["testable"] else NA,
                     "starship_elements": ",".join(c["elements"]),
                     "overlap_ge1gene": yn(c["ruleA"]),
                     "overlap_ge50pct": yn(c["ruleB"]),
                     "captain_genes": ",".join(c["captain"]),
                     "knowncluster_top_hit": c["knowncluster_top_hit"],
                     "knowncluster_description": c["knowncluster_description"],
                     "knowncluster_similarity": c["knowncluster_similarity"],
                     "genes": ",".join(g), "members": ",".join(g)})
    sio.write_tsv(a.out_dir / "sm_clusters.tsv.gz", rows)

    # ---------------- per type
    by_type = []
    types = sorted({c["type"] for c in clusters}, key=lambda t: (-sum(1 for c in clusters if c["type"] == t), t))
    for t in types + ["ALL"]:
        cs = [c for c in clusters if t == "ALL" or c["type"] == t]
        ov = bool(els)
        by_type.append({"type": t, "n_clusters": len(cs),
                        "n_overlap_ge1gene": sum(c["ruleA"] for c in cs) if ov else NA,
                        "n_overlap_ge50pct": sum(c["ruleB"] for c in cs) if ov else NA,
                        "n_with_captain": sum(1 for c in cs if c["captain"]) if a.captain_tblout else NA})
    sio.write_tsv(a.out_dir / "sm_by_type.tsv", by_type)

    # ---------------- per element (elements in SM-annotated strains)
    per_el = []
    sm_gene_keys = {(c["Short"], x) for c in clusters for x in c["gene_ids"]}
    for e in els:
        if e["Short"] not in sm_strains:
            continue
        eg = [g for g in gene.values() if g["Short"] == e["Short"] and g["element"] == e["element_id"]]
        eg.sort(key=lambda g: (g["contig"], g["start"], g["pid"]))
        cs = [c for c in clusters if c["Short"] == e["Short"] and e["element_id"] in c["elements"]]
        per_el.append({"element_id": e["element_id"], "label": e["label"], "confidence": e["confidence"],
                       "contig": e["contig"], "start": e["start"], "end": e["end"],
                       "n_genes_in_element": len(eg),
                       "n_sm_clusters_ge1gene": len(cs), "n_sm_clusters_ge50pct": sum(c["ruleB"] for c in cs),
                       "sm_clusters": ";".join(f"{c['cluster_id']}({c['type']};{c['n_in']}/{len(c['gene_ids'])} genes)"
                                               for c in cs),
                       "n_sm_cluster_genes_in_element": sum(1 for g in eg if (g["Short"], g["pid"]) in sm_gene_keys),
                       "captain_genes_in_element": ",".join(g["pid"] for g in eg if g["captain"])})
    if els:
        sio.write_tsv(a.out_dir / "sm_starship_overlap.tsv", per_el)

    # ---------------- permutation (shift each element on its own contig)
    perm_rows = []
    testable = [c for c in clusters if c["testable"]]
    if els and testable:
        clen = load_contig_lengths(a.contig_lengths)
        proxy_len: dict[tuple, int] = defaultdict(int)
        for (s, ctg), lst in by_sc.items():
            proxy_len[(s, ctg)] = max(g["end"] for g in lst)
        for c in clusters:
            proxy_len[(c["Short"], c["contig"])] = max(proxy_len[(c["Short"], c["contig"])], c["end"])
        for e in els:
            proxy_len[(e["Short"], e["contig"])] = max(proxy_len[(e["Short"], e["contig"])], e["end"])
        n_proxy = sum(1 for e in els if (e["Short"], e["contig"]) not in clen)
        C("permutation_contig_length", "contig_lengths" if a.contig_lengths and not n_proxy else "gene_extent",
          f"{n_proxy} element contig(s) use the largest gene/cluster/element end as contig length")

        def cluster_stats(placed):
            A = B = gin = 0
            for c in testable:
                k = sum(1 for g in c["placed"] if el_of((g["Short"], g["contig"]), g["mid"], placed))
                A += k >= 1
                B += bool(c["gene_ids"]) and k / len(c["gene_ids"]) >= 0.5
                gin += k
            return A, B, gin

        obs = cluster_stats(els_by_c)
        rng = random.Random(a.seed)
        null = []
        for _ in range(a.n_perm):
            placed = defaultdict(list)
            for e in els:
                key = (e["Short"], e["contig"])
                L = e["end"] - e["start"]
                s = rng.randint(1, max(1, clen.get(key, proxy_len[key]) - L))
                placed[key].append({"element_id": e["element_id"], "start": s, "end": s + L})
            null.append(cluster_stats(placed))
        for i, name in enumerate(("clusters_overlap_ge1gene", "clusters_overlap_ge50pct",
                                  "sm_cluster_genes_in_elements")):
            vals = sorted(x[i] for x in null)
            perm_rows.append({"statistic": name, "observed": obs[i], "n_perm": a.n_perm, "seed": a.seed,
                              "null_mean": round(sum(vals) / len(vals), 3) if vals else "",
                              "null_q95": vals[int(0.95 * (len(vals) - 1))] if vals else "",
                              "null_max": vals[-1] if vals else "",
                              "p_perm_greater": round((sum(1 for v in vals if v >= obs[i]) + 1) / (a.n_perm + 1), 5),
                              "p_perm_less": round((sum(1 for v in vals if v <= obs[i]) + 1) / (a.n_perm + 1), 5)})
        sio.write_tsv(a.out_dir / "sm_permutation.tsv", perm_rows)

    # ---------------- gene-level hypergeometric (strains with elements and SM rows)
    test_strains = sm_strains & el_strains
    univ = {k for k, g in gene.items() if k[0] in test_strains}
    in_el = {k for k in univ if gene[k]["element"]}
    sm_g = {(c["Short"], g["pid"]) for c in testable for g in c["placed"]}
    core_g = {(c["Short"], x) for c in testable for x in c["core_genes"] if (c["Short"], x) in gene}
    h_all = hyper(len(univ), len(in_el), len(sm_g), len(sm_g & in_el))
    h_core = hyper(len(univ), len(in_el), len(core_g), len(core_g & in_el))

    # ---------------- family view (core genes -> families), all SM-annotated strains
    fam_rows = []
    core_fams = defaultdict(set)
    gene_fams = set()
    for c in clusters:
        for g in c["placed"]:
            if g["family"]:
                gene_fams.add(g["family"])
        for x in c["core_genes"]:
            f = strain_m2r.get(f"{c['Short']}{sep}{x}", "")
            if f:
                core_fams[f].add(c["type"])
    for f in sorted(core_fams):
        fam_rows.append({"family": f, "sm_types": "+".join(sorted(core_fams[f])), "bin": freq.get(f, {}).get("bin", "NA"),
                         "strain_count": freq.get(f, {}).get("strain_count", ""), "direction": gl.get(f, "NA"),
                         "in_island": "yes" if f in island_fams else "no"})
    sio.write_tsv(a.out_dir / "sm_core_families.tsv.gz", fam_rows,
                  ["family", "sm_types", "bin", "strain_count", "direction", "in_island"])
    brk = []
    for label, fs in (("core_gene_families", set(core_fams)), ("all_cluster_gene_families", gene_fams)):
        cnt = Counter((freq.get(f, {}).get("bin", "NA"), gl.get(f, "NA"), f in island_fams) for f in fs)
        for (b, d, i), n in sorted(cnt.items()):
            brk.append({"gene_set": label, "bin": b, "direction": d, "in_island": "yes" if i else "no", "n_families": n})
    sio.write_tsv(a.out_dir / "sm_family_class_direction.tsv", brk,
                  ["gene_set", "bin", "direction", "in_island", "n_families"])
    cf_nc = set(core_fams) & nonc_f
    h_fam = hyper(len(nonc_f), len(cf_nc), len(island_fams & nonc_f), len(cf_nc & island_fams))
    gf_nc = gene_fams & nonc_f
    h_fam_all = hyper(len(nonc_f), len(gf_nc), len(island_fams & nonc_f), len(gf_nc & island_fams))

    # ---------------- summary
    S = []

    def M(metric, value, source):
        S.append({"metric": metric, "value": value, "source": source})

    M("sm_strains", ",".join(sorted(sm_strains)) if len(sm_strains) <= 20 else len(sm_strains),
      "strains with >= 1 row in the SM table")
    M("strains_with_sm_clusters", len(sm_strains), "")
    M("sm_strain_genes", len(gene), "genes of SM-annotated strains in gene_positions (+ extra)")
    M("sm_strain_genes_with_family", sum(1 for g in gene.values() if g["family"]), "member in tier1_cluster.tsv")
    M("sm_clusters", len(clusters), "rows of the SM table (sm_clusters.tsv.gz)")
    M("sm_cluster_genes", len({(c["Short"], x) for c in clusters for x in c["gene_ids"]}), "cluster genes")
    M("sm_cluster_genes_unplaced", n_unplaced, "listed cluster genes not in gene_positions")
    M("sm_core_genes", len({(c["Short"], x) for c in clusters for x in c["core_genes"]}), "core_genes column")
    M("sm_clusters_contig_edge", sum(1 for c in clusters if c["contig_edge"] == "True"), "contig_edge column")
    if els:
        M("starship_elements", len(els), "rows of the element table")
        M("sm_clusters_testable", len(testable), "clusters in strains with elements")
        M("sm_strain_cds_in_starships", len(in_el), "genes (strains with elements and SM rows) with midpoint inside an element")
        M("sm_clusters_overlap_ge1gene", sum(c["ruleA"] for c in clusters), "rule A (sm_clusters.tsv.gz)")
        M("sm_clusters_overlap_ge50pct", sum(c["ruleB"] for c in clusters), "rule B (sm_clusters.tsv.gz)")
    if a.captain_tblout:
        M("sm_clusters_with_captain", sum(1 for c in clusters if c["captain"]),
          f"captain hit E<={a.captain_evalue} in a cluster gene")
        M("sm_strain_captain_cds", sum(1 for g in gene.values() if g["captain"]), "SM-strain genes with a captain hit")
    if els:
        M("elements_with_sm_cluster_ge1gene", sum(1 for p in per_el if p["n_sm_clusters_ge1gene"]), "sm_starship_overlap.tsv")
        M("elements_with_sm_cluster_ge50pct", sum(1 for p in per_el if p["n_sm_clusters_ge50pct"]), "sm_starship_overlap.tsv")
    for pr in perm_rows:
        M(f"perm_{pr['statistic']}_observed", pr["observed"], "sm_permutation.tsv")
        M(f"perm_{pr['statistic']}_null_mean", pr["null_mean"], "sm_permutation.tsv")
        M(f"perm_{pr['statistic']}_p_greater", pr["p_perm_greater"], "sm_permutation.tsv")
        M(f"perm_{pr['statistic']}_p_less", pr["p_perm_less"], "sm_permutation.tsv")
    hyp = [("hypergeom_core_families_in_islands_noncore", h_fam),
           ("hypergeom_cluster_gene_families_in_islands_noncore", h_fam_all)]
    if els:
        hyp = [("hypergeom_sm_cluster_genes_in_starships", h_all),
               ("hypergeom_sm_core_genes_in_starships", h_core)] + hyp
    for nm, h in hyp:
        for k, v in h.items():
            M(f"{nm}_{k}", v, "universe: genes of strains with elements and SM rows" if "starships" in nm
              else "universe: non-core families; K = SM families (SM-strain genes); n = island families")
    M("sm_core_families", len(core_fams), "families of core genes (sm_core_families.tsv.gz)")
    M("sm_core_families_noncore", len(cf_nc), "")
    M("sm_core_families_noncore_in_island", len(cf_nc & island_fams), "")
    sio.write_tsv(a.out_dir / "summary.tsv", S)
    sio.write_tsv(a.out_dir / "input_checks.tsv", checks, ["check", "value", "note"])

    with open(a.out_dir / "PROVENANCE.yaml", "w") as fh:
        fh.write("derived_by: bin/pangenome_sm_overlap.py (port of NII full_v072_antismash/scripts/"
                 "antismash_sm_starship.py)\n"
                 f"date: {dt.date.today().isoformat()}\ncommand: {' '.join(sys.argv)}\n"
                 f"parameters: {{captain_evalue: {a.captain_evalue}, n_perm: {a.n_perm}, seed: {a.seed}, "
                 f"id_sep: '{sep}'}}\ninputs:\n")
        inputs = dict(paths)
        for k in ("sm_clusters", "starship_annotation", "extra_gene_positions", "captain_tblout", "gain_loss",
                  "contig_lengths"):
            if getattr(a, k):
                inputs[k] = getattr(a, k)
        for k, p in inputs.items():
            fh.write(f"  {k}:\n    path: {Path(p).resolve()}\n    sha256: {sio.sha256(Path(p))}\n")
    print(f"wrote {a.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
