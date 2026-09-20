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
from island_synteny import build_payload, select_islands  # noqa: E402
from island_synteny_template import ISLAND_SYNTENY_TEMPLATE  # noqa: E402
from pangenome_matrix import PresenceMatrix  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from pangenome_domain_enrichment import parse_domtblout  # noqa: E402


def load_positions(path: str) -> dict[tuple[str, str], list[tuple[str, int]]]:
    """(strain, family) -> [(contig, gene rank), ...], from
    family_positions.tsv[.zst].

    One row per (strain, family, COPY) (see
    bin/pangenome_build_family_positions.py's own stderr: "(strain, family,
    copy) positions written") -- a family with a paralog elsewhere in the
    genome has more than one entry here, so every copy is kept rather than
    letting a later row silently overwrite an earlier one.
    order_families_by_locus() picks the right copy using the island's own
    locus contig.
    """
    positions: dict[tuple[str, str], list[tuple[str, int]]] = {}
    with open_maybe_compressed(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                key = (row["Short"], row["family"])
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
    positions = load_positions(args.family_positions)

    family_domains = (parse_domtblout([args.domtblout], max_ievalue=args.domain_evalue)
                      if args.domtblout else None)

    payload = build_payload(
        island_rows, matrix, positions, project=args.project,
        min_strains=args.min_strains, top_n=args.top_islands,
        family_domains=family_domains,
    )

    # Escape `</` so a Pfam description or family ID cannot close the
    # <script> block early -- these strings come from HMM output and FASTA
    # headers, which this pipeline does not control.
    blob = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
    page = (ISLAND_SYNTENY_TEMPLATE
            .replace("__PROJECT_TITLE__", html.escape(args.project))
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
