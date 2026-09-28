"""Regression: lib/island_locus.py reproduces the feasibility prototype.

Spec validation plan item 2. The prototype is
NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/analysis/
island_locus_view_feasibility.py; its numbers on the 529-strain run
`rescue_freqpol_immitis_in_posadasii_out` (NovInvenio 91e3157) are in the
spec's Feasibility table and are pinned below (re-measured 2026-09-26 with
the prototype, --flank 5 --k 10 --top 50).

The prototype's definitions differ from the pipeline's in three stated ways,
so this test drives the library in "prototype mode": islands (not loci) are
the units, the exemplar is the island's example_strain, the locus block is
every member copy on the exemplar contig (locus_columns(block=None)), and
flanks intact accepts any copy pair (flank_pair(in_place=None)). Cell states
use StrainCells.base, which has no rescue or contig-break state.

Opt-in: set NOVINVENIO_LOCUS_REGRESSION_RUN to that run's output/pangenome
directory. About 3 minutes (two sort orders, two passes each).
"""
import collections
import csv
import os
import statistics
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
from island_locus import (  # noqa: E402
    ABSENT, IN_PLACE, flank_pair, group_loci, locus_columns, row_class, strain_cells,
)
from pangenome_island_loci import read_matrix_strains, scan_family_positions  # noqa: E402

RUN = os.environ.get("NOVINVENIO_LOCUS_REGRESSION_RUN")
pytestmark = pytest.mark.skipif(not RUN, reason="set NOVINVENIO_LOCUS_REGRESSION_RUN")

EXPECTED = {
    "size": {
        "both_sides_ge5": 21, "touch_end": 21, "median_left": 6.0, "median_right": 7.0,
        "core_softcore_pct": 66, "in_place": 131158, "elsewhere": 3965, "absent": 304476,
        "median_intact": 61.0, "median_empty_080": 0.0, "median_full": 2.0,
        "informative_080": 17, "informative_100": 0,
    },
    "strains": {
        "both_sides_ge5": 30, "touch_end": 7, "median_left": 15.0, "median_right": 16.0,
        "core_softcore_pct": 76, "in_place": 50604, "elsewhere": 297, "absent": 22101,
        "median_intact": 500.5, "median_empty_080": 52.5, "median_full": 245.5,
        "informative_080": 39, "informative_100": 39,
    },
}


def _islands():
    with open(Path(RUN) / "report_tables" / "islands_with_domains.tsv") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["locus_id"] not in ("-", "")]
    for r in rows:
        r["members"] = [f for f in r["member_families"].split(",") if f]
    return rows


def _measure(sort: str) -> dict:
    run = Path(RUN)
    fp = str(run / "family_positions.tsv.zst")
    qual = [r for r in _islands() if int(r["n_strains"]) >= 2]
    if sort == "size":
        qual.sort(key=lambda r: -int(r["island_size"]))
    else:
        qual.sort(key=lambda r: (-int(r["n_strains"]), -int(r["island_size"])))
    top = qual[:50]
    with open(run / "frequency_table.tsv") as fh:
        bins = {r["family"]: r["bin"] for r in csv.DictReader(fh, delimiter="\t")}
    scan = scan_family_positions(fp, contigs={(r["example_strain"], r["locus_contig"]) for r in top})
    cols = []
    for r in top:
        c = locus_columns(scan.orders.get((r["example_strain"], r["locus_contig"]), []),
                          r["members"], 5)
        if c is not None:
            cols.append(c)
    need = {f for c in cols for f in c.families}
    pos = scan_family_positions(fp, families=need).positions
    strains = read_matrix_strains(str(run / "presence_matrix.rescued.tsv"))
    lefts = [c.left_avail for c in cols]
    rights = [c.right_avail for c in cols]
    flank_bins = collections.Counter(bins.get(f, "?") for c in cols for f in c.left + c.right)
    states = collections.Counter()
    intact_n, empty80, empty100, full = [], [], [], []
    for c in cols:
        fams = c.families
        nl, nb = len(c.left), len(c.locus)
        n_int = e80 = e100 = nf = 0
        for s in strains:
            base = strain_cells(s, fams, pos, {}).base
            states.update(base[nl:nl + nb])
            pair = flank_pair(pos, s, list(c.left), list(c.right), nb) if c.left and c.right else None
            if pair is None:
                continue
            n_int += 1
            if row_class(base, nl, nb, True, 0.8) == "empty":
                e80 += 1
            if row_class(base, nl, nb, True, 1.0) == "empty":
                e100 += 1
            if row_class(base, nl, nb, True, 0.8) == "full":
                nf += 1
        intact_n.append(n_int)
        empty80.append(e80)
        empty100.append(e100)
        full.append(nf)
    tot_bins = sum(flank_bins.values())
    return {
        "both_sides_ge5": sum(1 for a, b in zip(lefts, rights) if a >= 5 and b >= 5),
        "touch_end": sum(1 for a, b in zip(lefts, rights) if a == 0 or b == 0),
        "median_left": statistics.median(lefts), "median_right": statistics.median(rights),
        "core_softcore_pct": round((flank_bins["core"] + flank_bins["soft_core"]) / tot_bins * 100),
        "in_place": states[IN_PLACE], "elsewhere": states["2"], "absent": states[ABSENT],
        "median_intact": statistics.median(intact_n),
        "median_empty_080": statistics.median(empty80), "median_full": statistics.median(full),
        "informative_080": sum(1 for e, f in zip(empty80, full) if e >= 10 and f >= 2),
        "informative_100": sum(1 for e, f in zip(empty100, full) if e >= 10 and f >= 2),
    }


@pytest.mark.parametrize("sort", ["size", "strains"])
def test_matches_feasibility_prototype(sort):
    assert _measure(sort) == EXPECTED[sort]


def test_locus_grouping_matches_m4():
    loci = group_loci(_islands())
    joined = sum(len(loc.variants) - 1 for loc in loci)
    assert (joined, sum(len(loc.variants) for loc in loci)) == (10579, 11880)
