#!/usr/bin/env python3
"""Pfam-domain enrichment test for accessory-island member families vs. the
CORRECT background: all families actually eligible for co-occurrence testing
(shell+cloud bins, matching pangenome_cooccurrence.py's own selection) --
never the whole genome, since core genes were never eligible for testing in
the first place and a whole-genome background would spuriously enrich for
"accessory-genome-typical" domains regardless of which specific island is
being tested.

Generalized from studies/fungi/Afumigatus_pangenome/bin/
summarize_island_functions.py's domain_enrichment() -- ported verbatim
(same Fisher-exact + BH-FDR statistics), split out of that study's combined
island-annotation-plus-enrichment script into a standalone enrichment-only
step (island annotation now lives in pangenome_report_tables.py, Task 7).

Usage:
  pangenome_domain_enrichment.py --significant_islands significant_islands.tsv \\
      --domtblout background_reps_vs_pfam.domtblout \\
      --frequency_table frequency_table.tsv \\
      --output island_pfam_enrichment.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys

import numpy as np
from scipy.stats import fisher_exact, false_discovery_control


def parse_domtblout(paths: list[str], max_ievalue: float = 1e-3) -> dict[str, set[str]]:
    """{family_id: {pfam_domain_name, ...}} from one or more hmmscan
    --domtblout files. Domain-level i-Evalue re-checked explicitly (can be
    looser than the sequence-level cutoff hmmscan's own -E already applied)."""
    hits: dict[str, set[str]] = {}
    for path in paths:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 13:
                    continue
                target_name, query_name = parts[0], parts[3]
                try:
                    i_evalue = float(parts[12])
                except ValueError:
                    continue
                if i_evalue > max_ievalue:
                    continue
                hits.setdefault(query_name, set()).add(target_name)
    return hits


def domain_enrichment(
    island_member_families: set[str],
    background_families: set[str],
    family_domains: dict[str, set[str]],
) -> list[dict]:
    """One-sided Fisher's exact test per Pfam domain, BH-FDR corrected across
    every domain tested. `island_member_families` is intersected with
    `background_families` first: singleton island members were never eligible
    for co-occurrence testing (shell/cloud only), so they're excluded from
    enrichment testing entirely -- a real, expected occurrence, not an
    error."""
    excluded = island_member_families - background_families
    island_member_families = island_member_families & background_families
    if excluded:
        print(
            f"domain_enrichment: {len(excluded)} island-member families are outside "
            "the eligible (shell+cloud) background -- most likely singleton families "
            "an island happened to include; excluded from the enrichment test "
            f"({len(island_member_families)} remain)", file=sys.stderr,
        )

    all_domains: set[str] = set()
    for family in island_member_families:
        all_domains |= family_domains.get(family, set())

    n_background = len(background_families)
    n_island = len(island_member_families)
    rows = []
    pvalues = []
    for domain in sorted(all_domains):
        families_with_domain = {f for f in background_families if domain in family_domains.get(f, ())}
        both = len(families_with_domain & island_member_families)
        island_only = n_island - both
        domain_only = len(families_with_domain) - both
        neither = n_background - n_island - domain_only
        _, p = fisher_exact([[both, island_only], [domain_only, neither]], alternative="greater")
        rows.append({
            "domain": domain,
            "n_with_domain_in_islands": both,
            "n_with_domain_in_background": len(families_with_domain),
            "n_island_families": n_island,
            "n_background_families": n_background,
            "fisher_p": p,
        })
        pvalues.append(p)

    if pvalues:
        qvalues = false_discovery_control(np.asarray(pvalues), method="bh")
        for row, q in zip(rows, qvalues):
            row["fdr_q"] = float(q)
    rows.sort(key=lambda r: r["fisher_p"])
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--significant_islands", required=True)
    ap.add_argument("--domtblout", required=True, action="append")
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    from pangenome_select_background_reps import select_background_families

    family_domains = parse_domtblout(args.domtblout)
    print(f"pangenome_domain_enrichment: {len(family_domains)} families with >=1 Pfam domain hit",
          file=sys.stderr)

    background = select_background_families(args.frequency_table)
    print(f"pangenome_domain_enrichment: {len(background)} eligible (shell+cloud) background families",
          file=sys.stderr)

    with open(args.significant_islands, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    island_member_families: set[str] = set()
    for row in rows:
        island_member_families.update(row["member_families"].split(","))

    enrichment = domain_enrichment(island_member_families, background, family_domains)
    with open(args.output, "w") as out:
        out.write(
            "domain\tn_with_domain_in_islands\tn_with_domain_in_background\t"
            "n_island_families\tn_background_families\tfisher_p\tfdr_q\n"
        )
        for r in enrichment:
            out.write(
                f"{r['domain']}\t{r['n_with_domain_in_islands']}\t{r['n_with_domain_in_background']}\t"
                f"{r['n_island_families']}\t{r['n_background_families']}\t"
                f"{r['fisher_p']:.3e}\t{r['fdr_q']:.3e}\n"
            )

    print(f"pangenome_domain_enrichment: {len(enrichment)} domains tested, "
          f"{sum(1 for r in enrichment if r['fdr_q'] < 0.05)} significant at FDR<0.05",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
