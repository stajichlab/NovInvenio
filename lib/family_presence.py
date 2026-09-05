"""Shared family-HMM presence-calling helpers for the novelty_discovery /
novelty_screen two-phase pathway (todo/novelty-discovery-screen.md).

Used by bin/novelty_presence_matrix.py (discovery phase, DISCOVERY_TARGET vs DISCOVERY_OUT) and
bin/novelty_screen.py (screen phase, families vs NEAR_INGROUP/BROAD_OUTGROUP) so both phases
call "present" the same way: present in a proteome if ANY target sequence there has a hit
with full-sequence E-value below the flat --default-family-evalue AND (profile coverage >= a
minimum fraction OR the aligned span is >= a minimum absolute residue count), both from that
SAME target (never mixing one target's E-value with another's coverage, and never pooling
coverage across targets -- see parse_domtblout).

The coverage-fraction floor alone (2026-09-03 review finding) turned out to punish long,
multi-domain proteins far more than it protects against its stated purpose: measured on a
real run, hits on HMMs >600 aa failed a 0.5 fraction floor 83% of the time, vs 7% for HMMs
<150 aa -- expected biology (a real ortholog of a multi-domain protein often only has one
domain conserved enough to align; the rest genuinely diverges past detection), not evidence
of promiscuous-domain false positives. The floor exists to stop a family HMM built from a
common, promiscuous domain (kinase, WD40, zinc finger, etc., shared across many unrelated
genes) from calling an unrelated outgroup protein present on the strength of that one shared
domain alone (ADR-0002 Q5). An absolute-residue alternative targets that concern directly --
a substantial aligned span (default 100 aa) is real evidence of homology regardless of what
fraction of a long multi-domain HMM it represents -- without penalizing long proteins purely
for being long. See CLAUDE.md/nextflow.config's `hmm_presence_min_residues`.

Presence is deliberately NOT gated by bin/calibrate_family_hmms.py's per-family
"negative-control" threshold. That threshold is derived from the best E-value the family's
own HMM scores against DISCOVERY_OUT, and presence-in-a-DISCOVERY_OUT-proteome was then
checked with `evalue < threshold` -- since the threshold IS the minimum observed E-value
across DISCOVERY_OUT, no DISCOVERY_OUT proteome could ever satisfy `E < min(E)`. That made
the outgroup-absence filter a structural no-op: every family with any DISCOVERY_OUT hit
(the vast majority of conserved families) was unconditionally called absent from every
DISCOVERY_OUT proteome, regardless of how strong the real orthology signal was (confirmed
against a real run: NCU00765/NCU00411/NCU01935 all had highly significant, annotated
DISCOVERY_OUT hits -- e.g. 3e-160 -- yet were called outgroup-absent because that very hit,
or an even better one from another DISCOVERY_OUT proteome, had set the threshold). See the
2026-09-03 review writeup. calibrate_family_hmms.py is kept only for the Nextflow process
interface / potential future refinement use; its output is not consumed here.
"""
from collections import defaultdict
from pathlib import Path


def _merged_span(intervals):
    """Total covered length of a set of (from, to) 1-based inclusive intervals, after
    merging overlaps -- so multiple domain hits on the same target aren't double-counted."""
    if not intervals:
        return 0
    ordered = sorted(intervals)
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return sum(end - start + 1 for start, end in merged)


def parse_domtblout(path, default_evalue, min_coverage, min_residues=0):
    """Return dict: query(HMM name) -> best qualifying full-sequence E-value.

    A query "qualifies" on a target when THAT SAME target's own hit clears both gates --
    full-sequence E-value < default_evalue AND (coverage (merged HMM-coordinate span /
    HMM length, on that target) >= min_coverage OR the merged aligned span itself is >=
    min_residues absolute aligned positions) -- never mixing one target's E-value with a
    different target's coverage, and never pooling coverage across different targets
    (either would let a family look present on evidence no single target actually
    supports). min_residues (default 0, i.e. disabled -- guarded explicitly, since a
    naive `covered >= 0` would trivially always be true) is the alternative to the
    fraction floor for long, multi-domain HMMs where a real but partial ortholog
    legitimately covers a small fraction of the full profile -- see module docstring. A
    query is present in this proteome if ANY target qualifies (mirrors
    bin/profile_to_matrix.py's semantics for the general mmseqs pathway) -- querying only
    the single best-E target's own coverage would wrongly call a family absent when its
    best-E hit happens to be a fragment but a different, slightly weaker-E target has
    full coverage (2026-09-03 review finding). A query with no qualifying target is
    omitted from the result.

    hmmsearch --domtblout columns (0-indexed): 0=target, 3=query, 5=qlen,
    6=full E-value, 15=hmm_from, 16=hmm_to.
    """
    per_target: dict[tuple, list] = {}
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 17:
                continue
            target = parts[0]
            query = parts[3]
            try:
                qlen = int(parts[5])
                full_e = float(parts[6])
                hmm_from = int(parts[15])
                hmm_to = int(parts[16])
            except ValueError:
                continue
            key = (query, target)
            if key not in per_target:
                per_target[key] = [qlen, full_e, []]
            per_target[key][2].append((hmm_from, hmm_to))

    result: dict[str, float] = {}
    for (query, _target), (qlen, full_e, intervals) in per_target.items():
        covered = _merged_span(intervals)
        coverage = covered / qlen if qlen > 0 else 0.0
        qualifies_coverage = coverage >= min_coverage or (min_residues > 0 and covered >= min_residues)
        if not (full_e < default_evalue and qualifies_coverage):
            continue
        best = result.get(query)
        if best is None or full_e < best:
            result[query] = full_e
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


def family_presence_by_proteome(domtblout_paths, default_evalue, min_coverage,
                                 min_residues=0,
                                 domtblout_suffixes=('.family.domtblout',
                                                      '.domtblout')):
    """Return dict: proteome_short -> set of family rep IDs "present" in it.

    A family is present in a proteome when any target sequence's own hit clears the
    flat E-value threshold (default_evalue -- see module docstring for why this is not
    per-family calibrated) and the coverage gate (fraction OR absolute residues) --
    see parse_domtblout.
    """
    presence = defaultdict(set)
    for dom_path in domtblout_paths:
        short = Path(dom_path).name
        for suffix in domtblout_suffixes:
            if short.endswith(suffix):
                short = short[:-len(suffix)]
                break
        for query in parse_domtblout(dom_path, default_evalue, min_coverage, min_residues):
            presence[short].add(query)
    return presence
