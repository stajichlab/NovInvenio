"""Shared family-HMM presence-calling helpers for the novelty_discovery /
novelty_screen two-phase pathway (todo/novelty-discovery-screen.md).

Used by bin/novelty_presence_matrix.py (discovery phase, TARGET vs DISC_OUT) and
bin/novelty_screen.py (screen phase, calibrated families vs NEAR_IN/BROAD_OUT) so both
phases call "present" the same way: a family HMM hit with full-sequence E-value below its
(optionally calibrated) per-family threshold AND profile coverage >= a minimum.
"""
from collections import defaultdict
from pathlib import Path


def parse_domtblout(path):
    """Return dict: query(HMM name) -> (best_full_evalue, best_coverage).

    Coverage = summed HMM-coordinate span / HMM length.
    hmmsearch --domtblout columns (0-indexed): 0=target, 3=query, 5=qlen,
    6=full E-value, 15=hmm_from, 16=hmm_to.
    """
    agg: dict[str, list] = {}
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 17:
                continue
            query = parts[3]
            try:
                qlen = int(parts[5])
                full_e = float(parts[6])
                hmm_from = int(parts[15])
                hmm_to = int(parts[16])
            except ValueError:
                continue
            span = max(0, hmm_to - hmm_from + 1)
            if query not in agg:
                agg[query] = [qlen, full_e, span]
            else:
                if full_e < agg[query][1]:
                    agg[query][1] = full_e
                agg[query][2] += span

    result = {}
    for query, (qlen, min_e, covered) in agg.items():
        coverage = covered / qlen if qlen > 0 else 0.0
        result[query] = (min_e, coverage)
    return result


def load_cluster_membership(cluster_tsv):
    """Return (rep_to_members, member_to_rep) from an mmseqs cluster TSV (rep<TAB>member)."""
    rep_to_members = defaultdict(list)
    member_to_rep = {}
    with open(cluster_tsv) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            rep_to_members[rep].append(member)
            member_to_rep[member] = rep
    return rep_to_members, member_to_rep


def load_family_thresholds(path):
    """Return dict: rep_id -> threshold_evalue, from CALIBRATE_FAMILY_HMMS's output."""
    thresholds = {}
    with open(path) as fh:
        next(fh, None)  # header
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                try:
                    thresholds[parts[0]] = float(parts[1])
                except ValueError:
                    continue
    return thresholds


def family_presence_by_proteome(domtblout_paths, family_thresholds, default_evalue,
                                 min_coverage, domtblout_suffixes=('.family.domtblout',
                                                                    '.domtblout')):
    """Return dict: proteome_short -> set of family rep IDs "present" in it.

    A family is present in a proteome when its domtblout hit clears both the
    per-family E-value threshold (falling back to default_evalue when uncalibrated)
    and the minimum profile coverage.
    """
    presence = defaultdict(set)
    for dom_path in domtblout_paths:
        short = Path(dom_path).name
        for suffix in domtblout_suffixes:
            if short.endswith(suffix):
                short = short[:-len(suffix)]
                break
        for query, (evalue, coverage) in parse_domtblout(dom_path).items():
            threshold = family_thresholds.get(query, default_evalue)
            if evalue < threshold and coverage >= min_coverage:
                presence[short].add(query)
    return presence
