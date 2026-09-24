#!/usr/bin/env python3
"""Co-occurrence (co-loss/co-gain) statistics between shell+cloud gene
families: frequency floor + Fisher/BH-FDR, outgroup-based gain/loss
polarization, an exact within-clade-stratified permutation-null replacement
as a cheap stand-in for a full phylogenetic correction.

Ported from NovInvenio_Investigations' Afumigatus_pangenome study
(bin/cooccurrence.py). Differences from that copy:

1. `config_parser` is imported directly from this repo's own lib/ (no
   cross-repo NOVINVENIO_ROOT lookup needed).
2. `--ingroup_label`/`--outgroup_label` replace hardcoded literal
   `"IN"`/`"OUT"` samplesheet GROUP values (NEXTFLOW_MIGRATION_NOTES.md
   section E).
3. `--n_perms` is DROPPED (not just deprecated): the within-clade
   permutation null is computed EXACTLY via `exact_stratified_pvalue`
   (a sum of independent per-clade hypergeometric draws, since the
   within-clade shuffle keeps every margin fixed) -- there is no sampling
   count left to honor, and the study's own notes recommended dropping the
   dead CLI flag in a fresh Nextflow port rather than carrying it forward
   for back-compat.

Usage:
  pangenome_cooccurrence.py --matrix presence_matrix.rescued.tsv \\
      --frequency_table frequency_table.tsv --config config.csv \\
      --output cooccurring_pairs.tsv
"""
from __future__ import annotations

import argparse
import itertools
import math
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact, false_discovery_control, hypergeom

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from ancestral_states import dollo_polarize, tree_direction  # noqa: E402
from compressed_io import open_maybe_compressed_write  # noqa: E402
from pangenome_matrix import PresenceMatrix  # noqa: E402
from pangenome_strain_inventory import read_representative_shorts  # noqa: E402


def presence_vector_within(matrix: PresenceMatrix, family: str, strains: list[str]) -> list[bool]:
    """Like PresenceMatrix.presence_vector, but restricted to a strain
    subset (e.g. ingroup-only) instead of always using matrix.strains."""
    return [matrix.is_present(family, s) for s in strains]


def strain_count_within(matrix: PresenceMatrix, family: str, strains: list[str]) -> int:
    return sum(presence_vector_within(matrix, family, strains))


def presence_bitset(matrix: PresenceMatrix, family: str, strains: list[str]) -> int:
    """Pack a family's presence/absence over `strains` into a single Python
    int (bit i = strains[i]) -- every family's bitset is built exactly once
    and reused for every pair via O(1)-ish bitwise ops
    (`fisher_counts_from_bits`), avoiding an O(pairs x strains) cost from
    recomputing a fresh presence vector per pair."""
    bits = 0
    for i, s in enumerate(strains):
        if matrix.is_present(family, s):
            bits |= 1 << i
    return bits


def fisher_counts_from_bits(a_bits: int, b_bits: int, n: int) -> tuple[int, int, int, int]:
    """(both, a_only, b_only, neither) 2x2-table counts for two presence
    bitsets over `n` strains, via bitwise AND/OR/XOR + `.bit_count()`."""
    full = (1 << n) - 1
    both = (a_bits & b_bits).bit_count()
    a_only = (a_bits & (full ^ b_bits)).bit_count()
    b_only = (b_bits & (full ^ a_bits)).bit_count()
    neither = (full ^ (a_bits | b_bits)).bit_count()
    return both, a_only, b_only, neither


def fisher_pvalue_from_counts(both: int, a_only: int, b_only: int, neither: int) -> float:
    _, p = fisher_exact([[both, a_only], [b_only, neither]], alternative="greater")
    return p


def quick_screen_pvalue(both: int, a_only: int, b_only: int, neither: int) -> float:
    """Fast, deliberately GENEROUS continuity-corrected normal approximation
    to the one-sided Fisher exact test -- used ONLY to decide whether a pair
    is worth an actual `scipy.fisher_exact` call, since an all-pairs exact-
    Fisher pass over tens of thousands of families is the real O(n^2)
    scaling cliff. The continuity correction (-0.5) biases the result toward
    LARGER p-values than the true one-sided exact test would give -- i.e.
    toward keeping a borderline pair for the exact test rather than
    discarding it -- so this must NEVER be used as the reported statistic,
    only as a pre-filter with a threshold looser than the real significance
    level (see `screen_alpha` in `find_cooccurring_pairs`, well above
    `fdr_alpha`). Returns 1.0 (never significant, always screened out) for a
    degenerate table with a zero margin, matching what the exact test would
    also report."""
    n = both + a_only + b_only + neither
    row1 = both + a_only
    col1 = both + b_only
    if n == 0 or row1 == 0 or col1 == 0 or row1 == n or col1 == n:
        return 1.0
    expected_both = row1 * col1 / n
    var = row1 * col1 * (n - row1) * (n - col1) / (n * n * (n - 1))
    if var <= 0:
        return 1.0
    z = (both - expected_both - 0.5) / math.sqrt(var)
    if z <= 0:
        return 1.0
    return 0.5 * math.erfc(z / math.sqrt(2))


def jaccard(a: list[bool], b: list[bool]) -> float:
    intersection = sum(1 for x, y in zip(a, b) if x and y)
    union = sum(1 for x, y in zip(a, b) if x or y)
    return intersection / union if union else 0.0


def fisher_pvalue(a: list[bool], b: list[bool]) -> float:
    both = sum(1 for x, y in zip(a, b) if x and y)
    a_only = sum(1 for x, y in zip(a, b) if x and not y)
    b_only = sum(1 for x, y in zip(a, b) if y and not x)
    neither = sum(1 for x, y in zip(a, b) if not x and not y)
    return fisher_pvalue_from_counts(both, a_only, b_only, neither)


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    if not pvalues:
        return []
    return list(false_discovery_control(pvalues, method="bh"))


def benjamini_hochberg_sparse(pvalues: list[float], total_m: int) -> list[float]:
    """BH-FDR q-values for a SPARSE subset of hypotheses that already
    survived `quick_screen_pvalue` -- a "streaming two-pass BH". Every
    hypothesis NOT in `pvalues` is treated as having p == 1.0 (the maximum
    possible value): the prefilter is built to never underestimate a
    survivor's true significance, so nothing it screens out could ever
    outrank a kept p-value in the full sort BH's procedure depends on.
    `total_m` is the TRUE total hypothesis count (e.g. `len(eligible) choose
    2`, not `len(pvalues)`) -- using `len(pvalues)` as m would under-correct,
    silently inflating the false-discovery rate."""
    k = len(pvalues)
    if k == 0:
        return []
    if total_m < k:
        raise ValueError(f"total_m ({total_m}) must be >= number of tested pairs ({k})")
    padded = np.concatenate([np.asarray(pvalues, dtype=np.float64), np.ones(total_m - k)])
    q_padded = false_discovery_control(padded, method="bh")
    return [float(x) for x in q_padded[:k]]


def polarize_direction(outgroup_present_count: int, outgroup_total: int) -> str:
    """A family present in EVERY outgroup strain implies the ingroup
    strains lacking it lost it ("loss"); a family absent from every
    outgroup strain implies the ingroup strains carrying it gained it
    ("gain"). With no outgroup strains to check, or with outgroup strains
    disagreeing (present in some but not all), the polarity can't be
    called -- "ambiguous"."""
    if outgroup_total == 0 or 0 < outgroup_present_count < outgroup_total:
        return "ambiguous"
    if outgroup_present_count == outgroup_total:
        return "loss"
    return "gain"


# Outgroup-frequency polarity (design decision Q4, 2026-09-20). Defaults calibrated on
# the 529-strain Coccidioides rescued matrix: of families present in >=99% of one
# species' strains, 89.6% (immitis -> posadasii, n=4,069) and 94.2% (posadasii ->
# immitis, n=4,213) reach >=0.90 in the other species, while only 10.6% / 32.6% are
# present in ALL of its strains. See docs/pangenome-assumptions.md.
DEFAULT_LOSS_MIN_FRAC = 0.9
DEFAULT_GAIN_MAX_FRAC = 0.1


def outgroup_fraction(outgroup_present_count: int, outgroup_total: int) -> float | None:
    """`asymmetry_a`: fraction of outgroup strains carrying the family. None when
    there are no outgroup strains (nothing to measure)."""
    if outgroup_total == 0:
        return None
    return outgroup_present_count / outgroup_total


def validate_polarity_thresholds(loss_min: float, gain_max: float) -> None:
    if not (0.0 <= gain_max < loss_min <= 1.0):
        raise ValueError(
            f"polarity thresholds need 0 <= gain_max < loss_min <= 1 "
            f"(got gain_max={gain_max}, loss_min={loss_min})")


def polarize_by_frequency(outgroup_present_count: int, outgroup_total: int,
                          loss_min: float = DEFAULT_LOSS_MIN_FRAC,
                          gain_max: float = DEFAULT_GAIN_MAX_FRAC) -> str:
    """Frequency-threshold version of polarize_direction(). "loss" when at least
    `loss_min` of outgroup strains carry the family, "gain" when at most
    `gain_max` do, "ambiguous" in between or with no outgroup. Unlike the strict
    rule, one missed call in a many-strain outgroup does not make a family
    unpolarisable. With a single outgroup strain the fraction is 0 or 1 and the
    two rules agree."""
    frac = outgroup_fraction(outgroup_present_count, outgroup_total)
    if frac is None:
        return "ambiguous"
    if frac >= loss_min:
        return "loss"
    if frac <= gain_max:
        return "gain"
    return "ambiguous"


def clade_composition(strains_present: list[str], clade_of_strain: dict[str, str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for s in strains_present:
        clade = clade_of_strain.get(s, "unknown")
        counts[clade] = counts.get(clade, 0) + 1
    return counts


@lru_cache(maxsize=200_000)
def _hypergeom_pmf_vector(n_c: int, k_c: int, m_c: int) -> tuple[float, ...]:
    """PMF of Hypergeom(N=n_c, K=k_c, n=m_c) over its support 0..min(k_c, m_c),
    as a hashable tuple so `exact_stratified_pvalue` below can cache it --
    n_c/k_c/m_c (a clade's strain count and family A/B counts within it)
    repeat across many family pairs in a real run, and the PMF itself is
    cheap to look up but not free to compute (scipy.stats.hypergeom.pmf).
    `maxsize` is bounded: distinct (n_c, k_c, m_c) keys scale roughly as
    sum_c (n_c+1)^2 for an uneven clade split, which at real-study scale can
    be a real, not theoretical, memory risk if left unbounded."""
    xs = np.arange(0, min(k_c, m_c) + 1)
    return tuple(hypergeom.pmf(xs, n_c, k_c, m_c))


def exact_stratified_pvalue(a: list[bool], b: list[bool], clades: list[str]) -> float:
    """Exact within-clade-stratified association test. Shuffling `b` WITHIN
    each clade group keeps every margin fixed: |a|, |b|, n, and each clade's
    own |a_c|, |b_c|, n_c. With every margin fixed, the one-sided
    (`greater`) Fisher/CMH statistic reduces to the overlap count ("both"),
    which under that shuffle is EXACTLY a sum of independent
    Hypergeom(n_c, |a_c|, |b_c|) draws, one per clade -- no simulation
    needed. p = P(sum_c X_c >= both_observed), the convolution of the
    per-clade PMFs. A clade with |a_c| == 0 or |b_c| == 0 contributes X_c ==
    0 with probability 1 (that clade can supply no overlap either way) and
    is skipped -- convolving with a point mass at 0 is a no-op."""
    both_observed = sum(1 for x, y in zip(a, b) if x and y)
    indices_by_clade: dict[str, list[int]] = {}
    for i, clade in enumerate(clades):
        indices_by_clade.setdefault(clade, []).append(i)

    dist = np.array([1.0])
    for indices in indices_by_clade.values():
        n_c = len(indices)
        k_c = sum(1 for i in indices if a[i])
        m_c = sum(1 for i in indices if b[i])
        if k_c == 0 or m_c == 0:
            continue
        dist = np.convolve(dist, np.array(_hypergeom_pmf_vector(n_c, k_c, m_c)))

    support = np.arange(len(dist))
    # min(1.0, ...): float64 summation over the convolved PMF's tail can
    # round fractionally above 1.0 when both_observed sits at or below the
    # support floor -- a p-value must never be reported > 1.0.
    return min(1.0, float(dist[support >= both_observed].sum()))


def find_cooccurring_pairs(
    matrix: PresenceMatrix,
    frequency_table: list[dict],
    clade_of_strain: dict[str, str],
    outgroup_presence: dict[str, tuple[int, int]],
    min_strain_count: int = 5,
    fdr_alpha: float = 0.05,
    strains: list[str] | None = None,
    screen_alpha: float = 0.2,
    species_tree=None,
    polarity_loss_min: float = DEFAULT_LOSS_MIN_FRAC,
    polarity_gain_max: float = DEFAULT_GAIN_MAX_FRAC,
) -> list[dict]:
    """`strains` restricts every ingroup statistic (the strain-count floor,
    Fisher presence vectors, clade composition) to that subset -- callers
    doing a real ingroup-vs-outgroup study MUST pass the ingroup-only strain
    list here, or outgroup strains silently inflate strain counts and
    contaminate clade composition. `outgroup_presence` maps family ->
    (outgroup_present_count, outgroup_total), always computed over the FULL
    matrix (it needs to see the outgroup strains) -- see polarize_direction
    for how the counts are used.

    Every family's presence is packed into a bitset exactly once
    (`presence_bitset`); each pair's 2x2 table comes from O(1)-ish bitwise
    ops (`fisher_counts_from_bits`); a fast analytic prefilter
    (`quick_screen_pvalue`, deliberately looser than `fdr_alpha`) skips the
    expensive exact `scipy.fisher_exact` call for the vast majority of
    clearly-non-significant pairs; BH-FDR correction
    (`benjamini_hochberg_sparse`) is applied against the TRUE total pair
    count, not just the survivor count.

    Invariant: the prefilter must NEVER be tighter than the caller's actual
    significance threshold -- `screen_alpha` is clamped up to at least
    `fdr_alpha`, so a permissive `fdr_alpha` always implies an equally
    permissive screen (the prefilter is purely a performance optimization
    and must never change what a given `fdr_alpha` would otherwise have
    kept)."""
    screen_alpha = max(screen_alpha, fdr_alpha)
    strains = matrix.strains if strains is None else strains
    n = len(strains)

    shell_cloud_families = [row["family"] for row in frequency_table if row["bin"] in ("shell", "cloud")]
    family_bits = {fam: presence_bitset(matrix, fam, strains) for fam in shell_cloud_families}
    eligible = [fam for fam in shell_cloud_families if family_bits[fam].bit_count() >= min_strain_count]
    clades = [clade_of_strain.get(s, "unknown") for s in strains]

    total_m = len(eligible) * (len(eligible) - 1) // 2
    print(
        f"cooccurrence: {len(eligible)} eligible families, {total_m} candidate "
        f"pairs to screen", file=sys.stderr,
    )
    # Deliberately only (fam_a, fam_b, p) here, NOT the four 2x2 counts too --
    # peak memory during screening matters at real scale, and the counts are
    # trivially recomputed from `family_bits` for the much smaller
    # FDR-survivor set below.
    survivors = []  # (fam_a, fam_b, p)
    progress_every = max(1, total_m // 20)  # ~20 progress lines regardless of scale
    for i, (fam_a, fam_b) in enumerate(itertools.combinations(eligible, 2), start=1):
        both, a_only, b_only, neither = fisher_counts_from_bits(
            family_bits[fam_a], family_bits[fam_b], n
        )
        if quick_screen_pvalue(both, a_only, b_only, neither) > screen_alpha:
            continue
        p = fisher_pvalue_from_counts(both, a_only, b_only, neither)
        survivors.append((fam_a, fam_b, p))
        if i % progress_every == 0:
            print(
                f"cooccurrence: screened {i}/{total_m} pairs "
                f"({len(survivors)} survivors so far)", file=sys.stderr,
            )

    if not survivors:
        return []
    qvalues = benjamini_hochberg_sparse([s[2] for s in survivors], total_m)

    fdr_survivors = [
        (fam_a, fam_b, p, q) for (fam_a, fam_b, p), q in zip(survivors, qvalues) if q < fdr_alpha
    ]
    print(
        f"cooccurrence: {len(survivors)} pairs cleared the prefilter, "
        f"{len(fdr_survivors)} clear FDR ({fdr_alpha}) -- running the exact "
        "within-clade-stratified test", file=sys.stderr,
    )
    perm_progress_every = max(1, len(fdr_survivors) // 20)

    # Dollo reconstruction is per FAMILY but consumed per PAIR, and a family appears in
    # many pairs -- memoise or this becomes O(pairs) instead of O(families).
    ingroup_set = set(strains or matrix.strains)
    dollo_cache: dict[str, tuple] = {}

    def dollo_for(fam):
        if fam not in dollo_cache:
            # The FULL matrix on purpose: the tree spans ingroup + outgroup, while
            # `strains` is deliberately ingroup-only for every other statistic here.
            carriers = {t for t in matrix.strains if matrix.is_present(fam, t)}
            res = dollo_polarize(species_tree, carriers)
            dollo_cache[fam] = (res, tree_direction(res, carriers, ingroup_set))
        return dollo_cache[fam]

    results = []
    for j, (fam_a, fam_b, p, q) in enumerate(fdr_survivors, start=1):
        a_bits, b_bits = family_bits[fam_a], family_bits[fam_b]
        both, a_only, b_only, neither = fisher_counts_from_bits(a_bits, b_bits, n)
        vec_a = [bool((a_bits >> i) & 1) for i in range(n)]
        vec_b = [bool((b_bits >> i) & 1) for i in range(n)]
        strains_present_a = [s for s, present in zip(strains, vec_a) if present]
        out_count_a, out_total_a = outgroup_presence.get(fam_a, (0, 0))
        union = both + a_only + b_only
        results.append({
            "family_a": fam_a,
            "family_b": fam_b,
            "jaccard": both / union if union else 0.0,
            "fisher_p": p,
            "fdr_q": q,
            "permutation_p": exact_stratified_pvalue(vec_a, vec_b, clades),
            "direction_a": polarize_direction(out_count_a, out_total_a),
            "clade_composition": clade_composition(strains_present_a, clade_of_strain),
            # Reported ALONGSIDE direction_a (design decision Q4), like the tree columns.
            "asymmetry_a": outgroup_fraction(out_count_a, out_total_a),
            "direction_a_freq": polarize_by_frequency(
                out_count_a, out_total_a, polarity_loss_min, polarity_gain_max),
        })
        if species_tree is not None:
            res, direction = dollo_for(fam_a)
            results[-1].update({
                # Reported ALONGSIDE direction_a, never replacing it (issue #140), so a
                # run can be compared old-vs-new before anything is retired.
                "direction_a_tree": direction,
                "gain_node": res.gain_node or "",
                "n_loss_events": res.n_loss_events,
                "loss_clades": ";".join(",".join(c) for c in res.loss_clades),
                "asr_method": res.method,
            })
        if j % perm_progress_every == 0:
            print(
                f"cooccurrence: exact stratified test on {j}/{len(fdr_survivors)} "
                "FDR-significant pairs", file=sys.stderr,
            )
    return results


def unlabelled_clade_fraction(strains: list[str], clade_of_strain: dict[str, str]) -> float:
    """Fraction of `strains` with no usable TaxonGroup label."""
    if not strains:
        return 0.0
    unlabelled = sum(1 for s in strains if not (clade_of_strain.get(s) or "").strip())
    return unlabelled / len(strains)


def warn_if_unstratified(
    strains: list[str], clade_of_strain: dict[str, str], threshold: float = 0.5
) -> bool:
    """The permutation null shuffles WITHIN clade, so strains with no
    TaxonGroup all land in one "unknown" pool and are shuffled freely -- i.e.
    unstratified. Returns True when the warning fired."""
    frac = unlabelled_clade_fraction(strains, clade_of_strain)
    if frac <= threshold:
        return False
    print(
        f"WARNING: {frac:.0%} of the {len(strains)} analysed strains have an empty "
        "TaxonGroup. Those strains all fall into a single 'unknown' group, so the "
        "clade-stratified permutation null is effectively an UNSTRATIFIED shuffle "
        "for them -- permutation_p does not control for clonal-clade structure. "
        "Populate TaxonGroup (e.g. DAPC/Mash clade assignments) before "
        "interpreting permutation_p as a phylogenetic control.",
        file=sys.stderr,
    )
    return True


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--ingroup_label", default="IN")
    ap.add_argument("--outgroup_label", default="OUT")
    ap.add_argument("--species_tree", default=None,
                    help="Rooted Newick strain tree spanning ingroup + outgroup. When "
                         "given, adds phylogeny-aware gain/loss columns "
                         "(direction_a_tree, gain_node, n_loss_events, loss_clades, "
                         "asr_method) ALONGSIDE the count-based direction_a, which is "
                         "unchanged. Tip labels must match the matrix strain IDs "
                         "exactly. Build it per study from markers (nf_phyling); a "
                         "taxonomy or mash tree is not a substitute. See "
                         "lib/ancestral_states.py for the Dollo assumption and its "
                         "known bias against horizontally-transferred families.")
    ap.add_argument("--inventory",
                    help="strain_inventory.tsv from pangenome_dereplicate_strains.py; when "
                         "given, the ingroup strain set is further restricted to "
                         "is_representative == 1 strains")
    ap.add_argument("--min_strain_count", type=int, default=5)
    ap.add_argument("--fdr_alpha", type=float, default=0.05)
    ap.add_argument(
        "--screen_alpha", type=float, default=0.2,
        help="analytic-prefilter threshold (must stay looser than --fdr_alpha) "
        "-- a pair only gets the expensive exact Fisher test if the fast "
        "normal-approximation screen clears this; see quick_screen_pvalue",
    )
    ap.add_argument(
        "--polarity_loss_min_frac", type=float, default=DEFAULT_LOSS_MIN_FRAC,
        help="direction_a_freq = loss when at least this fraction of outgroup strains "
        "carry the family (default %(default)s; see DEFAULT_LOSS_MIN_FRAC)",
    )
    ap.add_argument(
        "--polarity_gain_max_frac", type=float, default=DEFAULT_GAIN_MAX_FRAC,
        help="direction_a_freq = gain when at most this fraction of outgroup strains "
        "carry the family (default %(default)s; symmetric choice, not calibrated)",
    )
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    try:
        validate_polarity_thresholds(args.polarity_loss_min_frac, args.polarity_gain_max_frac)
    except ValueError as e:
        ap.error(str(e))

    sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
    from config_parser import parse_config  # noqa: E402

    samples = parse_config(args.config)
    clade_of_strain = {s.short: s.taxon_group for s in samples}
    outgroup_shorts = [s.short for s in samples if s.group == args.outgroup_label]
    ingroup_shorts = [s.short for s in samples if s.group == args.ingroup_label]
    if args.inventory:
        reps = set(read_representative_shorts(args.inventory))
        ingroup_shorts = [s for s in ingroup_shorts if s in reps]
        print(f"Dereplication: {len(ingroup_shorts)} representative ingroup strains",
              file=sys.stderr)
    if not ingroup_shorts:
        print("ERROR: no ingroup strains left after --config/--inventory filtering",
              file=sys.stderr)
        sys.exit(1)
    warn_if_unstratified(ingroup_shorts, clade_of_strain)

    matrix = PresenceMatrix.from_tsv(args.matrix)
    with open(args.frequency_table) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        frequency_table = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]

    # Intersect with the matrix's own columns rather than trusting the
    # samplesheet's full outgroup list: a matrix built with --groups IN (or
    # any other ingroup-only build) has no outgroup columns at all, and
    # computing outgroup_total from the samplesheet alone would silently
    # treat "not present in a column that doesn't exist" as "absent from
    # the outgroup" -- producing a confident, wrong "gain" call for every
    # family instead of the "ambiguous" that a missing outgroup should
    # produce.
    outgroup_in_matrix = [s for s in outgroup_shorts if s in matrix.strains]
    if outgroup_shorts and not outgroup_in_matrix:
        print(
            "WARNING: none of the samplesheet's outgroup strains are columns in "
            f"{args.matrix} (built with an ingroup-only --groups?) -- "
            "gain/loss polarization will be 'ambiguous' for every family.",
            file=sys.stderr,
        )
    outgroup_total = len(outgroup_in_matrix)
    outgroup_presence = {
        fam: (
            sum(1 for s in outgroup_in_matrix if matrix.is_present(fam, s)),
            outgroup_total,
        )
        for fam in matrix.families
    }

    species_tree = None
    if args.species_tree:
        from Bio import Phylo
        species_tree = Phylo.read(args.species_tree, "newick")
        tips = {t.name for t in species_tree.get_terminals()}
        missing = [t for t in matrix.strains if t not in tips]
        if missing:
            # Fail loudly: silently dropping strains would quietly change every
            # reconstruction without changing any column that shows it.
            raise SystemExit(
                f"--species_tree is missing {len(missing)} matrix strain(s): "
                f"{sorted(missing)[:10]} -- tip labels must match strain IDs exactly"
            )

    pairs = find_cooccurring_pairs(
        matrix, frequency_table, clade_of_strain, outgroup_presence,
        args.min_strain_count, args.fdr_alpha,
        strains=ingroup_shorts, screen_alpha=args.screen_alpha,
        species_tree=species_tree,
        polarity_loss_min=args.polarity_loss_min_frac,
        polarity_gain_max=args.polarity_gain_max_frac,
    )
    # asymmetry_a / direction_a_freq are always appended after clade_composition
    # (consumers read columns by name). New columns are APPENDED, and only when --species_tree was given, so a run
    # without a tree produces a byte-identical header to before (issue #140).
    tree_cols = ["direction_a_tree", "gain_node", "n_loss_events", "loss_clades",
                 "asr_method"] if species_tree is not None else []
    with open_maybe_compressed_write(args.output) as fh:
        fh.write("\t".join(
            ["family_a", "family_b", "jaccard", "fisher_p", "fdr_q", "permutation_p",
             "direction_a", "clade_composition", "asymmetry_a", "direction_a_freq"]
            + tree_cols) + "\n")
        for row in pairs:
            line = (
                f"{row['family_a']}\t{row['family_b']}\t{row['jaccard']:.4f}\t"
                f"{row['fisher_p']:.2e}\t{row['fdr_q']:.2e}\t{row['permutation_p']:.4f}\t"
                f"{row['direction_a']}\t{row['clade_composition']}\t"
                f"{'' if row['asymmetry_a'] is None else format(row['asymmetry_a'], '.4f')}\t"
                f"{row['direction_a_freq']}"
            )
            for c in tree_cols:
                line += f"\t{row[c]}"
            fh.write(line + "\n")


if __name__ == "__main__":
    main()
