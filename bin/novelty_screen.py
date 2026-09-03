#!/usr/bin/env python3
"""
Phase 2 of the two-phase targeted novelty pipeline (todo/novelty-discovery-screen.md):
refine novelty_discovery's candidate list against a broader phylogenetic screen.

novelty_discovery (phase 1) calls a family "in-group specific" if its HMM has no hit in the
small DISCOVERY_OUT reference panel. This phase re-searches the SAME family HMMs against two
broader proteome sets and reclassifies each phase-1 candidate into one of three categories:

    target_specific  -- not found in NEAR_INGROUP or BROAD_OUTGROUP.  Unique to the target genome(s).
    clade_specific   -- found in NEAR_INGROUP, not found in BROAD_OUTGROUP.  Shared with close
                        relatives, but still absent from distant lineages.
    false_novelty    -- found in BROAD_OUTGROUP.  The family is widespread, not clade-specific;
                        removed from the screened candidate list (but kept, labelled, in the
                        screened matrix for visibility).

"Found in" a group means at least one hit clears the flat --default-family-evalue and
--min-coverage, both from the same target sequence, in at least one proteome of that group --
the same presence rule bin/novelty_presence_matrix.py applies for the DISCOVERY_OUT panel
(see lib/family_presence.py). Per-family calibration (--family-thresholds) is accepted for
Nextflow process interface compatibility but is NOT used to gate presence -- it is circular
(derived from the very hits it would judge); see lib/family_presence.py's module docstring.

Singleton candidates (targets with no multi-member family, issue #52) get the SAME
screening: the extended singleton query phase 1 searched (singletons + their own
paralogs, see bin/extend_singleton_query.py) is re-searched against NEAR_INGROUP/
BROAD_OUTGROUP and filtered with the identical significance + paralog-competition logic
(lib/singleton_presence.py), merged into the same near_in_presence/broad_out_presence sets
classify() reads. Without this, a singleton candidate always defaulted to target_specific
regardless of its true broader-set presence -- family HMM search structurally can't see it
(it has no family).

Every row from the discovery matrix is carried through -- extended with NEAR_INGROUP/BROAD_OUTGROUP
presence columns and a novelty_category -- matching the established convention that the
presence matrix holds every scored row while candidates.txt carries the filtered list (see
workflows/search.nf's BUILD_PRESENCE_MATRIX docstring). Rows that were not phase-1 candidates
get an empty novelty_category (they were never novelty candidates, so classifying them into
one of the three categories is meaningless), but still get real NEAR_INGROUP/BROAD_OUTGROUP presence
values -- the family HMM search covers every family, not just phase-1 survivors.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import parse_config  # noqa: E402
from family_presence import (  # noqa: E402
    family_presence_by_proteome,
    load_cluster_membership,
)
from singleton_presence import (  # noqa: E402
    load_paralog_info,
    parse_pairwise_tsv,
    proteome_short_from_hits_filename,
    score_singleton_hits,
)


def load_candidates(path):
    """Return {protein_id: source_proteome} from a candidates.txt (source::protein)."""
    candidates = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            source, _, protein = line.partition('::')
            candidates[protein if protein else source] = source if protein else ''
    return candidates


def read_matrix(path):
    """Return (header_fields, rows) from a TSV with a header line."""
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        rows = []
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            rows.append(line.split('\t'))
    return header, rows


def classify(rep, near_in_presence, broad_out_presence):
    """Return one of 'target_specific' | 'clade_specific' | 'false_novelty'."""
    in_broad_out = any(rep in present for present in broad_out_presence.values())
    if in_broad_out:
        return 'false_novelty'
    in_near_in = any(rep in present for present in near_in_presence.values())
    if in_near_in:
        return 'clade_specific'
    return 'target_specific'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--discovery-matrix', required=True, dest='discovery_matrix',
                    help='novelty_discovery presence_matrix.tsv')
    ap.add_argument('--discovery-candidates', required=True, dest='discovery_candidates',
                    help='novelty_discovery candidates.txt (source::protein)')
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv',
                    help='families_cluster.tsv (rep_id<TAB>member_id) from novelty_discovery')
    ap.add_argument('--near-in-domtblout', nargs='*', default=[], dest='near_in_domtblout',
                    help='hmmsearch domtblout files, family HMMs vs NEAR_INGROUP proteomes')
    ap.add_argument('--broad-out-domtblout', nargs='*', default=[], dest='broad_out_domtblout',
                    help='hmmsearch domtblout files, family HMMs vs BROAD_OUTGROUP proteomes')
    ap.add_argument('--family-thresholds', default=None, dest='family_thresholds',
                    help='Accepted for interface compatibility with the Nextflow process '
                         '(novelty_discovery\'s CALIBRATE_FAMILY_HMMS output) but NOT used '
                         'to gate presence -- see lib/family_presence.py\'s module '
                         'docstring. Presence always uses the flat --default-family-evalue.')
    ap.add_argument('--default-family-evalue', type=float, default=1e-3,
                    dest='default_family_evalue',
                    help='Default E-value threshold for families without calibration')
    ap.add_argument('--min-coverage', type=float, default=0.5, dest='min_coverage',
                    help='Minimum profile coverage for family presence')
    ap.add_argument('--near-in-singleton-hits', nargs='*', default=[],
                    dest='near_in_singleton_hits',
                    help='Parsed singleton-query pairwise hits vs NEAR_INGROUP proteomes')
    ap.add_argument('--broad-out-singleton-hits', nargs='*', default=[],
                    dest='broad_out_singleton_hits',
                    help='Parsed singleton-query pairwise hits vs BROAD_OUTGROUP proteomes')
    ap.add_argument('--paralog-cutoffs', nargs='+', default=[], dest='paralog_cutoffs',
                    help='paralog_cutoffs.tsv files from the DISCOVERY_TARGET self-vs-self '
                         'search (same ones novelty_presence_matrix.py used for phase 1) -- '
                         'supplies the paralog_protein_ID column for the paralog-competition '
                         'filter; the e-value column is not used (see --singleton-evalue).')
    ap.add_argument('--singleton-evalue', type=float, default=1e-5, dest='singleton_evalue',
                    help='Flat e-value significance cutoff for singleton hits, applied to '
                         'every hit (not per-query paralog-derived)')
    ap.add_argument('--paralog-competition-scope', choices=['proteome', 'target'],
                    default='target', dest='paralog_competition_scope',
                    help='Same semantics as bin/build_presence_matrix.py / '
                         "bin/novelty_presence_matrix.py (default: 'target', matches "
                         "nextflow.config's pipeline default)")
    ap.add_argument('--config', required=True, help='Analysis description CSV')
    ap.add_argument('--output-matrix', required=True, dest='output_matrix')
    ap.add_argument('--output-candidates', required=True, dest='output_candidates')
    args = ap.parse_args()

    header, rows = read_matrix(args.discovery_matrix)
    candidates = load_candidates(args.discovery_candidates)
    rep_to_members, member_to_rep = load_cluster_membership(args.cluster_tsv)
    singleton_ids = {rep for rep, members in rep_to_members.items() if len(members) == 1}

    samples = parse_config(args.config)
    near_in_shorts = sorted(s.short for s in samples if s.group == 'NEAR_INGROUP')
    broad_out_shorts = sorted(s.short for s in samples if s.group == 'BROAD_OUTGROUP')

    near_in_presence = family_presence_by_proteome(
        args.near_in_domtblout, args.default_family_evalue, args.min_coverage)
    broad_out_presence = family_presence_by_proteome(
        args.broad_out_domtblout, args.default_family_evalue, args.min_coverage)

    # Singleton screening (issue #52): same paralog-aware filtering as phase 1, merged
    # into the same presence dicts classify() reads -- a singleton and a family both
    # just contribute "this rep is present in this proteome" evidence.
    paralog_of = load_paralog_info(args.paralog_cutoffs)

    def _singleton_hits(paths):
        hits = []
        for path in paths:
            short = proteome_short_from_hits_filename(path)
            for query_id, target_id, evalue in parse_pairwise_tsv(path):
                hits.append((query_id, target_id, evalue, short))
        return hits

    near_in_singleton_presence, _ = score_singleton_hits(
        _singleton_hits(args.near_in_singleton_hits), singleton_ids,
        paralog_of, args.singleton_evalue, args.paralog_competition_scope)
    broad_out_singleton_presence, _ = score_singleton_hits(
        _singleton_hits(args.broad_out_singleton_hits), singleton_ids,
        paralog_of, args.singleton_evalue, args.paralog_competition_scope)

    for short, ids in near_in_singleton_presence.items():
        near_in_presence.setdefault(short, set()).update(ids)
    for short, ids in broad_out_singleton_presence.items():
        broad_out_presence.setdefault(short, set()).update(ids)

    pid_idx = header.index('protein_id')

    # The discovery matrix's own header already spans every proteome short ID in the config
    # (see bin/novelty_presence_matrix.py, which is not restricted to DISCOVERY_TARGET/DISCOVERY_OUT) --
    # NEAR_INGROUP/BROAD_OUTGROUP columns are usually already present there, always 0 (phase 1 never
    # searches them). Update those columns in place rather than appending duplicates; only
    # append a short ID that is genuinely new to this matrix.
    new_near_in = [s for s in near_in_shorts if s not in header]
    new_broad_out = [s for s in broad_out_shorts if s not in header]
    out_header = header + new_near_in + new_broad_out + ['novelty_category']

    out_rows = []
    n_by_category = {'target_specific': 0, 'clade_specific': 0, 'false_novelty': 0}
    screened_candidates = []
    for row in rows:
        pid = row[pid_idx]
        rep = member_to_rep.get(pid, pid)
        row_by_col = dict(zip(header, row))
        for s in near_in_shorts:
            row_by_col[s] = '1' if rep in near_in_presence.get(s, ()) else '0'
        for s in broad_out_shorts:
            row_by_col[s] = '1' if rep in broad_out_presence.get(s, ()) else '0'

        if pid in candidates:
            category = classify(rep, near_in_presence, broad_out_presence)
            n_by_category[category] += 1
            if category != 'false_novelty':
                source = candidates[pid]
                screened_candidates.append(f'{source}::{pid}' if source else pid)
        else:
            category = ''
        row_by_col['novelty_category'] = category
        out_rows.append([row_by_col.get(col, '0') for col in out_header])

    with open(args.output_matrix, 'w') as out:
        out.write('\t'.join(out_header) + '\n')
        for row in out_rows:
            out.write('\t'.join(row) + '\n')

    with open(args.output_candidates, 'w') as out:
        if screened_candidates:
            out.write('\n'.join(screened_candidates) + '\n')

    print(f"Screened {len(candidates)} phase-1 candidates: "
          f"{n_by_category['target_specific']} target-specific, "
          f"{n_by_category['clade_specific']} clade-specific, "
          f"{n_by_category['false_novelty']} false novelties removed", file=sys.stderr)


if __name__ == '__main__':
    main()
