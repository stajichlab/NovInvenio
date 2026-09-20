"""Build the island synteny view's JSON payload.

For one accessory island: rows = strains, columns = member families in locus
order, filled cell = family present in that strain. Deletion breakpoints show
up as vertical edges.

No I/O at import time; every function here takes already-parsed rows and returns
deterministic values except for `build_payload`, which embeds a wall-clock
timestamp. The whole payload is unit-testable without touching a pipeline run.
"""
from __future__ import annotations

from datetime import datetime, timezone

from pfam_classes import CLASS_LABELS, dominant_class


def _as_int(value, default: int = 0) -> int:
    """Parse a TSV cell to int, tolerating '' and the '-' missing sentinel."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _located(rows: list[dict]) -> list[dict]:
    """Islands with resolved locus coordinates.

    An unplaced island (`locus_id == '-'`) has no column order, so the view
    cannot draw it. This filter is the same in both select_islands and
    build_payload to prevent the two from drifting.
    """
    return [r for r in rows if r.get("locus_id", "-") not in ("-", "", None)]


def select_islands(rows: list[dict], min_strains: int = 2,
                   top_n: int = 50) -> list[dict]:
    """Islands worth drawing, largest first.

    Two filters, in this order:

    1. A locus must be resolved (`locus_id != '-'`), since the view is
       locus-ordered and an unplaced island has no column order.
    2. The island must be carried by at least `min_strains` strains. On the
       real genus_vs_ureesii run 17,218 of 27,836 located islands (62%) are
       single-strain, and all 20 of the largest were -- an unfiltered size
       ranking is therefore entirely strain-private content, and a
       single-strain island's grid is one filled row with no breakpoint in
       it. This is the same filter PR #110 added to report.md's top-islands
       table, for the same reason.

    Truncation to `top_n` happens AFTER filtering, so the filter never costs
    slots in the output.
    """
    located = _located(rows)
    kept = [r for r in located if _as_int(r.get("n_strains")) >= min_strains]
    kept.sort(key=lambda r: -_as_int(r.get("island_size")))
    return kept[:top_n]


def collapse_haplotypes(presence_rows: dict[str, list[bool]]) -> list[dict]:
    """Collapse identical strain rows into distinct presence patterns.

    Returns one entry per distinct pattern: {"pattern", "count", "strains"},
    sorted lexicographically on the 0/1 pattern string. That sort is
    equivalent to a Hamming-distance leaf ordering for this data and is what
    makes deletion breakpoints align into visible vertical edges.

    Row collapsing is the load-bearing decision of this view: without it, a
    530-strain study draws 530 rows for a 20-56 column island and nothing is
    legible. With it, the row count is the number of haplotypes, typically
    tens.
    """
    by_pattern: dict[str, list[str]] = {}
    for strain, vector in presence_rows.items():
        pattern = "".join("1" if v else "0" for v in vector)
        by_pattern.setdefault(pattern, []).append(strain)
    return [
        {"pattern": p, "count": len(s), "strains": sorted(s)}
        for p, s in sorted(by_pattern.items())
    ]


def order_families_by_locus(families: list[str],
                            positions: dict[tuple[str, str], object],
                            strain: str,
                            contig: str | None = None) -> list[str]:
    """Member families in locus order within `strain`, on `contig`.

    Gene ORDER is what makes this view meaningful -- the question it answers
    is "where does the block break?" -- so the columns must follow the
    island's own layout, taken from the example strain's gene ranks in
    family_positions.tsv.

    `positions[(strain, family)]` may be either a bare rank (int, the legacy
    single-copy shape) or a list of `(contig, rank)` pairs -- family_positions.tsv
    has one row per (strain, family, COPY) (see
    bin/pangenome_build_family_positions.py), so a family with a paralog
    elsewhere in the genome has more than one entry. When a list is given, a
    family's locus-order key is the MINIMUM rank among its copies ON `contig`
    (the island's own locus contig); if it has no copy on that contig, its
    minimum rank across ALL contigs is used as a fallback rather than
    dropping the column. Using an off-contig rank verbatim would previously
    let a stray paralog on an unrelated, alphabetically-later contig drag
    the whole column to the end of the grid -- confirmed on a real 529-strain
    study where 2,623 of 4,509,470 (strain, family) pairs are multi-copy.

    A family with no position in this strain at all keeps a column, sorted
    last by name. Dropping it would silently shrink the island and change
    the breakpoint the reader sees.
    """
    def rank_key(family: str):
        entries = positions.get((strain, family))
        if entries is None:
            return None
        if isinstance(entries, int):
            return entries
        if not entries:
            return None
        on_contig = [r for c, r in entries if contig is not None and c == contig]
        if on_contig:
            return min(on_contig)
        return min(r for _, r in entries)

    ranked = []
    unranked = []
    for f in families:
        key = rank_key(f)
        if key is None:
            unranked.append(f)
        else:
            ranked.append((key, f))
    unranked.sort()
    return [f for _, f in sorted(ranked)] + unranked


def build_payload(island_rows: list[dict], matrix, positions: dict,
                  project: str, min_strains: int = 2,
                  top_n: int = 50,
                  family_domains: dict[str, set[str]] | None = None,
                  species_of: dict[str, str] | None = None) -> dict:
    """The JSON payload embedded in the island synteny page.

    `matrix` is a lib.pangenome_matrix.PresenceMatrix; `positions` maps
    (strain, family) -> gene rank. `family_domains` is
    bin/pangenome_domain_enrichment.py's `parse_domtblout()` result --
    {family_id: {domain_names}} from the pfam.domtblout FAMILY_PFAM_SCAN
    produces, keyed by family representative ID (the same IDs that appear in
    `member_families`).

    The payload records `project` because family IDs are NOT stable across
    pipeline re-runs -- a family's ID is its mmseqs cluster representative,
    and mmseqs does not pick the same one twice (measured: 73% ID overlap
    between two runs of the same study). The page is only ever valid for the
    run that produced it, and must never be joined to another run's output by
    family ID.

    Each island also carries `family_classes`/`family_domains` -- one entry
    per column, same order as `families` -- so the page's glyph strip can
    colour each column by ITS OWN Pfam class rather than a single island-wide
    average. `family_domains` is optional and a family may legitimately be
    absent from it (the Pfam scan covers the shell/cloud background set, not
    every family in an island): both cases resolve to the "unannotated"
    class and an empty domain string rather than raising, so a run with no
    domain-enrichment data at all still produces a valid page.

    `species_of` (issue #119) is an OPTIONAL trailing {Short: Species} map --
    the analysis config CSV's own Species column, keyed by Short, as
    bin/pangenome_island_synteny.py's optional --config resolves it via
    lib/config_parser.py::parse_config(). It powers the page's "by species"
    row sort (spec: "immitis/posadasii band" -- grouping strains of the same
    species together, NOT a phylogeny-ordered sort; this pipeline has no
    strain tree, see todo/pangenome-phylogeny-aware-gain-loss.md). The
    payload always carries a top-level `species` key so the page can check
    it unconditionally: when `species_of` is omitted, or a strain in the
    presence matrix has no entry in it (e.g. added to the run after the
    config CSV was last touched), that strain is simply absent from the map
    rather than raising -- the species sort option is hidden client-side
    when the map is empty, exactly mirroring how the novelty report hides
    its category filter when the payload carries no category data.
    """
    fam_domain_map = family_domains or {}
    species_map = species_of or {}
    # Separate locus-filter exclusions from top_n truncation to count them honestly.
    located = _located(island_rows)
    qualifying = [r for r in located if _as_int(r.get("n_strains")) >= min_strains]
    qualifying.sort(key=lambda r: -_as_int(r.get("island_size")))
    selected = qualifying[:top_n]
    strains = sorted(matrix.strains)

    islands = []
    for row in selected:
        families = [f for f in row.get("member_families", "").split(",") if f]
        example = row.get("example_strain", "")
        locus_contig = row.get("locus_contig", "")
        ordered = order_families_by_locus(families, positions, example,
                                          contig=locus_contig or None)
        presence_rows = {
            s: [matrix.is_present(f, s) for f in ordered] for s in strains
        }
        domains = [d for d in row.get("pfam_domains", "").split(",") if d]
        family_domain_sets = [sorted(fam_domain_map.get(f, set())) for f in ordered]
        islands.append({
            "locus_id": row.get("locus_id", "-"),
            "locus_contig": row.get("locus_contig", ""),
            "locus_start": _as_int(row.get("locus_start"), -1),
            "locus_end": _as_int(row.get("locus_end"), -1),
            "example_strain": example,
            "size": len(ordered),
            "n_strains": _as_int(row.get("n_strains")),
            "families": ordered,
            "family_classes": [dominant_class(doms) for doms in family_domain_sets],
            "family_domains": [",".join(doms) for doms in family_domain_sets],
            "domains": domains,
            "dominant_class": dominant_class(domains),
            "haplotypes": collapse_haplotypes(presence_rows),
        })

    return {
        "project": project,
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "strains": strains,
        "classes": CLASS_LABELS,
        "islands": islands,
        "n_islands_total": len(located),
        "n_islands_excluded": len(located) - len(qualifying),
        "n_islands_truncated": len(qualifying) - len(selected),
        "species": {s: species_map[s] for s in strains if s in species_map},
    }
