#!/usr/bin/env python3
"""Render the island synteny / presence-absence view (issue #116).

One self-contained HTML page: for each selected accessory island, a
presence/absence grid with rows = strains collapsed into distinct haplotypes and
columns = member families in locus order.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import open_maybe_compressed  # noqa: E402
from config_parser import parse_config  # noqa: E402
from island_synteny import build_payload, select_islands  # noqa: E402
from island_synteny_template import ISLAND_SYNTENY_TEMPLATE  # noqa: E402
from pangenome_matrix import PresenceMatrix  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from pangenome_domain_enrichment import parse_domtblout  # noqa: E402


def load_positions(path: str,
                   families: "set[str] | None" = None) -> dict[tuple[str, str], list[tuple[str, int]]]:
    """(strain, family) -> [(contig, gene rank), ...], from
    family_positions.tsv[.zst].

    One row per (strain, family, COPY) (see
    bin/pangenome_build_family_positions.py's own stderr: "(strain, family,
    copy) positions written") -- a family with a paralog elsewhere in the
    genome has more than one entry here, so every copy is kept rather than
    letting a later row silently overwrite an earlier one.
    order_families_by_locus() picks the right copy using the island's own
    locus contig.

    `families`, when given, skips rows whose family is not in the set WHILE
    STREAMING (issue #118 follow-up: this file has one row per (strain,
    family, copy) -- 4.5M rows / 54,421 families on the real 529-strain
    genus_vs_ureesii study -- and building the whole dict regardless of how
    few islands are drawn was still the dominant cost after the presence-
    matrix fix). The filter is applied by FAMILY ONLY, never by contig or
    copy -- a family kept because it's in `families` keeps EVERY one of its
    rows, on every contig. Filtering copies would silently break
    order_families_by_locus()'s off-contig fallback (it falls back to a
    family's global minimum rank when the island's own locus contig has no
    copy), which needs every copy of a selected family present to work.
    `families=None` (the default) is byte-identical to today's unfiltered
    behaviour -- no other caller passes this argument.
    """
    positions: dict[tuple[str, str], list[tuple[str, int]]] = {}
    with open_maybe_compressed(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                fam = row["family"]
                if families is not None and fam not in families:
                    continue
                key = (row["Short"], fam)
                entry = (row["contig"], int(row["rank"]))
            except (KeyError, ValueError):
                continue
            positions.setdefault(key, []).append(entry)
    return positions


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--islands_with_domains", required=True)
    ap.add_argument("--presence_matrix", required=True)
    ap.add_argument("--family_positions", required=True)
    ap.add_argument("--project", required=True)
    ap.add_argument("--min_strains", type=int, default=2)
    ap.add_argument("--top_islands", type=int, default=50)
    ap.add_argument("--domtblout", default=None)
    ap.add_argument("--domain_evalue", type=float, default=1e-3,
                    help="Domain-level i-Evalue ceiling for --domtblout, "
                    "matching REPORT_TABLES/DOMAIN_ENRICHMENT's "
                    "--pangenome_pfam_domain_evalue so the glyph strip and "
                    "the sidebar chip agree on the same page.")
    ap.add_argument("--config", default=None,
                    help="Optional analysis samplesheet CSV (GROUP,Species,...,Short,...). "
                    "Threads a {Short: Species} map into the payload for the page's "
                    "'by species' row sort (issue #119). Omitting it is never an "
                    "error -- the sort option is simply not offered.")
    ap.add_argument(
        "--diagnostics_banner", default=None,
        help="Optional pangenome_diagnostics.py diagnostics_banner.html "
        "snippet (issue #134) -- inserted at the top of the page body, "
        "before any results.",
    )
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    with open_maybe_compressed(args.islands_with_domains) as fh:
        island_rows = list(csv.DictReader(fh, delimiter="\t"))

    # issue #118: load only the families the islands we're actually going to
    # draw reference, not the whole matrix -- select_islands() applies the
    # same min-strains filter and top-N truncation build_payload() applies
    # internally (that duplication is deliberate and cheap; see this file's
    # docstring notes and lib/island_synteny.py's build_payload docstring),
    # so the family set computed here always matches what build_payload()
    # would draw. On the real 529-strain genus_vs_ureesii study this took
    # peak RSS from 5.54 GB (54,421 families materialised for 50 islands) to
    # a small fraction of that -- 298 distinct families actually referenced.
    selected = select_islands(island_rows, min_strains=args.min_strains,
                              top_n=args.top_islands)
    needed_families = {
        f for row in selected for f in row.get("member_families", "").split(",") if f
    }
    matrix = PresenceMatrix.from_tsv(args.presence_matrix, families=needed_families)
    # Same set, same reasoning, for family_positions.tsv -- see load_positions()'s
    # docstring for why this is safe (filtered by family only, never by copy).
    positions = load_positions(args.family_positions, families=needed_families)

    family_domains = (parse_domtblout([args.domtblout], max_ievalue=args.domain_evalue)
                      if args.domtblout else None)

    species_of = ({s.short: s.species for s in parse_config(args.config)}
                  if args.config else None)

    payload = build_payload(
        island_rows, matrix, positions, project=args.project,
        min_strains=args.min_strains, top_n=args.top_islands,
        family_domains=family_domains, species_of=species_of,
    )

    # Escape `</` so a Pfam description or family ID cannot close the
    # <script> block early -- these strings come from HMM output and FASTA
    # headers, which this pipeline does not control.
    blob = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
    diagnostics_banner_html = (
        Path(args.diagnostics_banner).read_text() if args.diagnostics_banner else ""
    )
    page = (ISLAND_SYNTENY_TEMPLATE
            .replace("__PROJECT_TITLE__", html.escape(args.project))
            .replace("<!--__DIAGNOSTICS_BANNER__-->", diagnostics_banner_html)
            .replace("/*__PAYLOAD__*/", blob))
    Path(args.output).write_text(page)

    print(f"pangenome_island_synteny: {len(payload['islands'])} islands drawn, "
          f"{payload['n_islands_excluded']} excluded (< {args.min_strains} "
          f"strains), {payload['n_islands_truncated']} truncated "
          f"(qualified but past --top_islands {args.top_islands}), "
          f"wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
