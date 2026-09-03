#!/usr/bin/env python3
"""BUSCO presence-recovery metric for the family-parameter sweep (2026-09-03 follow-up to
ADR-0002 Q8 / issue #6).

bin/busco_family_recovery.py measures *clustering* quality (do a BUSCO's Complete copies
land in one gene family?) and is insensitive to --hmm_presence_evalue/--hmm_presence_cov/
--hmm_presence_min_residues -- those parameters only act downstream, at the hmmsearch
presence-calling step, after clustering has already happened. This script measures *that*
step directly: a single-copy BUSCO ortholog is, by definition, present in every species
BUSCO found a Complete copy in, so for each such species (other than the family's own
member proteomes, where presence is trivial by construction) the presence_matrix.tsv row
for that family SHOULD show 1. A 0 there is a real false-absence caused by the
presence-calling parameters (or by a domtblout/alignment quirk), not by clustering -- this
is the curation-free control that --hmm_presence_cov/--hmm_presence_min_residues sweeping
needs and bin/busco_family_recovery.py does not provide.

Only BUSCOs whose cluster-MEMBER copies (i.e. seed-group species) land in exactly one
profiled family are scored -- a copy from a non-member (outgroup) species is expected and
NOT a disqualifier (that's the population under test), but a BUSCO split across two
DIFFERENT families among its member copies would make "the family's presence row"
ambiguous, and that failure mode is clustering's to own, not presence-calling's.

IMPORTANT: pass --tables for OUTGROUP (or DISCOVERY_OUT) species, not (only) seed-group
species. A seed-group species' presence in its own family is trivial by clustering
construction (it IS a member), so it tests nothing about hmmsearch-based presence-calling
-- this script skips those pairs automatically and will report 0 scorable species-pairs
if every --tables species happens to be in the seed group.

Recovery is reported overall AND bucketed by the BUSCO marker's own alignment length
(full_table.tsv's Length column) -- because a metric aggregated only over BUSCO markers is
a claim about the length distribution of THIS BUSCO lineage set, not the whole proteome or
the whole novelty-candidate population, and BUSCO markers are hand-picked conserved
single-copy orthologs, not a random sample of gene lengths. Pass --reference-domtblout
(any FAMILY_HMMSEARCH domtblout from the same run) to also report how the *eligible*
BUSCO set's length distribution compares to the run's whole queried-family population, so a
sweep reader can see directly whether this run's BUSCO control is representative before
trusting its recovery number -- rather than assuming it, see the module docstring's
2026-09-03 finding that for a real N. crassa run BUSCO Complete gene lengths (median 460 aa,
IQR 304-702) closely tracked the whole proteome's (median 462 aa, IQR 302-686), with BUSCO
modestly underrepresenting the shortest (<150 aa) bucket (3.8% vs 6.7%) -- the one bucket
that was already unaffected by the coverage floor, so low-risk for this particular use.
"""
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from busco import (  # noqa: E402
    load_busco_map,
    load_member_to_rep,
    load_profiled_reps,
    read_matrix_rows,
)

LENGTH_BUCKETS = (
    ('<150aa', 0, 150),
    ('150-300aa', 150, 300),
    ('300-600aa', 300, 600),
    ('>600aa', 600, None),
)


def bucket_for(length):
    if length is None:
        return 'unknown'
    for name, lo, hi in LENGTH_BUCKETS:
        if length >= lo and (hi is None or length < hi):
            return name
    return 'unknown'  # pragma: no cover -- unreachable, LENGTH_BUCKETS covers [0, inf)


def recovered_buscos(busco_map, member_to_rep, min_species):
    """Return dict: busco_id -> rep_id, restricted to BUSCOs whose cluster-MEMBER copies
    (i.e. seed-group species -- non-members are expected and NOT a disqualifier here,
    unlike bin/busco_family_recovery.py's stricter "recovered" verdict) all land in
    exactly one profiled family. A non-member copy (e.g. an OUTGROUP/DISCOVERY_OUT
    species, never fed to clustering) is exactly the population this script exists to
    test via hmmsearch presence, not a clustering failure -- see module docstring."""
    result = {}
    for busco_id, by_species in busco_map.items():
        if len(by_species) < min_species:
            continue
        reps = {member_to_rep[pid] for pid in by_species.values() if pid in member_to_rep}
        if len(reps) == 1:
            result[busco_id] = next(iter(reps))
    return result


def any_row_for_rep(rep, member_to_rep, matrix_rows):
    """Return the matrix row (proteome_short -> '0'/'1') for any member of this family --
    every member's row carries the same family-level presence outside its own source
    proteome (see lib/family_presence.py), so any one representative row suffices."""
    if rep in matrix_rows:
        return matrix_rows[rep]
    for member, member_rep in member_to_rep.items():
        if member_rep == rep and member in matrix_rows:
            return matrix_rows[member]
    return None


def score_presence_recovery(busco_map, lengths, member_to_rep, matrix_rows, min_species):
    """Return (per_pair rows, summary dict, bucket_summary dict)."""
    rep_of_busco = recovered_buscos(busco_map, member_to_rep, min_species)

    per_pair = []
    hit = miss = no_row = 0
    bucket_totals = {name: [0, 0] for name, _, _ in LENGTH_BUCKETS}  # [hit, total]
    bucket_totals['unknown'] = [0, 0]

    for busco_id, rep in sorted(rep_of_busco.items()):
        row = any_row_for_rep(rep, member_to_rep, matrix_rows)
        by_species = busco_map[busco_id]
        length = lengths.get(busco_id)
        bucket = bucket_for(length)
        own_species = {sp for sp, pid in by_species.items()
                       if member_to_rep.get(pid) == rep and pid in matrix_rows}
        for species, protein_id in sorted(by_species.items()):
            if species in own_species:
                continue  # trivially present by construction, not a presence-calling test
            if row is None:
                no_row += 1
                verdict = 'no_row'
            else:
                present = row.get(species) == '1'
                if present:
                    hit += 1
                    verdict = 'present'
                else:
                    miss += 1
                    verdict = 'absent'
                bucket_totals[bucket][1] += 1
                if present:
                    bucket_totals[bucket][0] += 1
            per_pair.append({
                'busco_id': busco_id, 'species': species, 'family_rep': rep,
                'length': '' if length is None else length, 'length_bucket': bucket,
                'verdict': verdict,
            })

    scored = hit + miss
    summary = {
        'recovered_buscos_scored': len(rep_of_busco),
        'species_pairs_scored': scored,
        'species_pairs_present': hit,
        'species_pairs_absent': miss,
        'species_pairs_no_row': no_row,
        'presence_recovery_rate': (hit / scored) if scored else None,
    }
    bucket_summary = {}
    for name, (bhit, btotal) in bucket_totals.items():
        bucket_summary[name] = {
            'n': btotal,
            'recovery_rate': (bhit / btotal) if btotal else None,
        }
    return per_pair, summary, bucket_summary


def length_representativeness(eligible_lengths, reference_domtblout):
    """Compare the eligible-BUSCO length distribution to the whole run's queried-family
    (HMM qlen) length distribution from a reference domtblout -- see module docstring."""
    ref_lengths = []
    seen = set()
    with open(reference_domtblout) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            query = parts[3]
            if query in seen:
                continue
            try:
                qlen = int(parts[5])
            except ValueError:
                continue
            seen.add(query)
            ref_lengths.append(qlen)

    def summarize(lens):
        lens = sorted(x for x in lens if x is not None)
        n = len(lens)
        if n == 0:
            return {'n': 0, 'median': '', 'p25': '', 'p75': '', 'frac_gt600': '', 'frac_lt150': ''}
        return {
            'n': n,
            'median': lens[n // 2],
            'p25': lens[n // 4],
            'p75': lens[3 * n // 4],
            'frac_gt600': round(sum(1 for x in lens if x > 600) / n, 3),
            'frac_lt150': round(sum(1 for x in lens if x < 150) / n, 3),
        }

    return {'eligible_busco': summarize(eligible_lengths), 'reference': summarize(ref_lengths)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tables', nargs='*', default=[],
                    help='BUSCO full_table.tsv per proteome, as SHORT=path')
    ap.add_argument('--busco-map', default=None, dest='busco_map',
                    help='consolidated busco_id<TAB>species<TAB>protein_id[<TAB>status] TSV')
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv',
                    help='mmseqs *_cluster.tsv (rep<TAB>member)')
    ap.add_argument('--families', required=True, help='families.tsv (profiled families)')
    ap.add_argument('--matrix', required=True, help='presence_matrix.tsv from this run')
    ap.add_argument('--min-species', type=int, default=2, dest='min_species',
                    help='BUSCO must be Complete in >= this many species to be eligible')
    ap.add_argument('--reference-domtblout', default=None, dest='reference_domtblout',
                    help='optional FAMILY_HMMSEARCH domtblout from the same run, to report '
                         "whether the eligible BUSCO set's length distribution matches the "
                         "run's whole queried-family population (see module docstring)")
    ap.add_argument('--output', required=True, help='per (busco_id, species) results TSV')
    ap.add_argument('--summary', default=None,
                    help='summary TSV (default: *.summary.tsv next to --output)')
    args = ap.parse_args()

    if not args.tables and not args.busco_map:
        sys.exit('provide --tables SHORT=path ... and/or --busco-map')

    busco_map, lengths = load_busco_map(args.busco_map, args.tables)
    profiled_reps = load_profiled_reps(args.families)
    member_to_rep = load_member_to_rep(args.cluster_tsv, profiled_reps)
    matrix_rows = read_matrix_rows(args.matrix)

    per_pair, summary, bucket_summary = score_presence_recovery(
        busco_map, lengths, member_to_rep, matrix_rows, args.min_species)

    fields = ['busco_id', 'species', 'family_rep', 'length', 'length_bucket', 'verdict']
    with open(args.output, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(per_pair)

    summary_path = args.summary or (str(Path(args.output).with_suffix('')) + '.summary.tsv')
    with open(summary_path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['metric', 'value'])
        for k, v in summary.items():
            w.writerow([k, '' if v is None else v])
        for bucket, stats in bucket_summary.items():
            w.writerow([f'bucket_{bucket}_n', stats['n']])
            w.writerow([f'bucket_{bucket}_recovery_rate',
                       '' if stats['recovery_rate'] is None else stats['recovery_rate']])
        if args.reference_domtblout:
            rep = length_representativeness(
                (lengths.get(bid) for bid in recovered_buscos(
                    busco_map, member_to_rep, args.min_species)),
                args.reference_domtblout)
            for group, stats in rep.items():
                for k, v in stats.items():
                    w.writerow([f'length_repr_{group}_{k}', v])

    rate = summary['presence_recovery_rate']
    print(f"BUSCO presence recovery: "
          f"{'n/a' if rate is None else f'{rate:.3f}'} "
          f"({summary['species_pairs_present']}/{summary['species_pairs_scored']} "
          f"species-pairs, {summary['recovered_buscos_scored']} recovered BUSCOs)",
          file=sys.stderr)
    if summary['species_pairs_scored'] == 0 and summary['recovered_buscos_scored'] > 0:
        print("NOTE: 0 scorable species-pairs despite recovered BUSCOs -- this metric is "
              "only informative for species whose proteins are NOT themselves family/cluster "
              "members (presence there is hmmsearch-based, the thing this parameter "
              "controls). If every --tables species is in the seed group (DISCOVERY_TARGET "
              "for novelty_discovery, or the ingroup for --cluster_tool mmseqs), presence in "
              "all of them is trivial by clustering construction and there is nothing left "
              "to test. Pass BUSCO tables for OUTGROUP/DISCOVERY_OUT species instead (or in "
              "addition).", file=sys.stderr)


if __name__ == '__main__':
    main()
