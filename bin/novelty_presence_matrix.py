#!/usr/bin/env python3
"""Build the novelty discovery presence matrix from two evidence sources.

1. Family HMM search: multi-member target families are profiled (famsa +
   hmmbuild) and searched via hmmsearch against all proteomes (target +
   DISCOVERY_OUT).  A family is "present" in a proteome if the domtblout has a
   hit with full-sequence E-value below the (optionally calibrated)
   per-family threshold AND profile coverage >= --min-coverage.

2. Singleton pairwise search: singletons (target proteins not in any
   multi-member family) are searched via phmmer/diamond/blast against the
   DISCOVERY_OUT proteomes.  A singleton is "present" in a proteome if the
   pairwise search reports a hit with E-value below the singleton
   threshold.  The singleton's source proteome is always present.

The script merges both sources into a single presence_matrix.tsv with
columns:

    protein_id, source_proteome, <sorted proteome shorts>

And a candidates.txt with lines:

    source_proteome::protein_id

A protein/family is a novelty candidate if it is present in >=
--target-min-frac of the target proteomes AND absent from all DISCOVERY_OUT
proteomes (--disc-out-max-frac, default 0.0 = strictly absent).
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import parse_config  # noqa: E402
from family_presence import (  # noqa: E402
    load_cluster_membership,
    load_family_thresholds,
    parse_domtblout,
)


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def parse_pairwise_tsv(path):
    """Return list of (query_id, target_id, evalue) tuples from parsed hits TSV."""
    hits = []
    with open(path) as fh:
        header = fh.readline()
        del header
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            try:
                evalue = float(parts[2])
            except ValueError:
                continue
            hits.append((parts[0], parts[1], evalue))
    return hits


def load_protein_map(path):
    """Return dict: protein_id -> proteome_short."""
    mapping = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--family-domtblout', nargs='+', default=[],
                    dest='family_domtblout',
                    help='Per-proteome hmmsearch domtblout files for family HMMs')
    ap.add_argument('--singleton-hits', nargs='+', default=[],
                    dest='singleton_hits',
                    help='Parsed pairwise hits TSVs for singletons')
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv',
                    help='mmseqs *_cluster.tsv (rep_id<TAB>member_id)')
    ap.add_argument('--protein-map', required=True, dest='protein_map',
                    help='protein_id<TAB>proteome_short TSV')
    ap.add_argument('--config', required=True, help='Analysis description CSV')
    ap.add_argument('--family-thresholds', default=None, dest='family_thresholds',
                    help='Calibrated family thresholds TSV (rep_id<TAB>threshold)')
    ap.add_argument('--default-family-evalue', type=float, default=1e-3,
                    dest='default_family_evalue',
                    help='Default E-value threshold for families without calibration')
    ap.add_argument('--singleton-evalue', type=float, default=1e-5,
                    dest='singleton_evalue',
                    help='E-value threshold for singleton pairwise hits')
    ap.add_argument('--min-coverage', type=float, default=0.5,
                    dest='min_coverage',
                    help='Minimum profile coverage for family presence')
    ap.add_argument('--target-min-frac', type=float, default=0.75,
                    dest='target_min_frac',
                    help='Minimum fraction of target proteomes a family/protein must be present in')
    ap.add_argument('--disc-out-max-frac', type=float, default=0.0,
                    dest='disc_out_max_frac',
                    help='Max fraction of DISCOVERY_OUT proteomes a candidate may be present in')
    ap.add_argument('--output-matrix', required=True, dest='output_matrix')
    ap.add_argument('--output-candidates', required=True, dest='output_candidates')
    ap.add_argument('--output-evalues', default=None, dest='output_evalues',
                    help='Optional sidecar TSV, same shape as --output-matrix, holding the '
                         'family-HMM or singleton hit e-value per (protein, proteome) cell '
                         '(empty for absent/self/unsearched-in-phase-1 columns) -- report-only '
                         'evidence, does not affect candidate calling.')
    args = ap.parse_args()

    samples = parse_config(args.config)
    short_to_group = {s.short: s.group for s in samples}
    all_shorts = sorted(short_to_group.keys())

    rep_to_members, member_to_rep = load_cluster_membership(args.cluster_tsv)
    protein_to_proteome = load_protein_map(args.protein_map)

    # Identify multi-member families (>= 2 members) and singletons
    multi_reps = {rep for rep, members in rep_to_members.items() if len(members) >= 2}
    singleton_reps = {rep for rep, members in rep_to_members.items() if len(members) == 1}

    # --- Family HMM presence ---
    family_thresholds = {}
    if args.family_thresholds:
        family_thresholds = load_family_thresholds(args.family_thresholds)

    # family_presence[proteome_short] = set of family rep IDs present
    # family_evalue[(proteome_short, rep_id)] = the qualifying hit's full-seq e-value
    family_presence = defaultdict(set)
    family_evalue = {}

    for dom_path in args.family_domtblout:
        short = Path(dom_path).name
        for suffix in ('.family.domtblout', '.domtblout'):
            if short.endswith(suffix):
                short = short[:-len(suffix)]
                break
        hits = parse_domtblout(dom_path)
        for query, (evalue, coverage) in hits.items():
            threshold = family_thresholds.get(query, args.default_family_evalue)
            if evalue < threshold and coverage >= args.min_coverage:
                family_presence[short].add(query)
                family_evalue[(short, query)] = evalue

    # --- Singleton pairwise presence ---
    # singleton_presence[proteome_short] = set of singleton protein IDs present
    # singleton_evalue[(proteome_short, protein_id)] = best (lowest) qualifying hit e-value
    singleton_presence = defaultdict(set)
    singleton_evalue = {}

    for hits_path in args.singleton_hits:
        name = Path(hits_path).name
        for suffix in ('.parsed.tsv', '.tsv', '.parsed'):
            if name.endswith(suffix):
                name = name[:-len(suffix)]
                break
        # PARSE_HITS names its output '<query_id>_vs_<target_id>.parsed.tsv' (the same
        # convention workflows/search.nf uses) -- the singleton search's query_id is always
        # the literal 'singletons', so the proteome short is whatever follows '_vs_'.
        short = name.rsplit('_vs_', 1)[1] if '_vs_' in name else name
        for query_id, target_id, evalue in parse_pairwise_tsv(hits_path):
            if evalue < args.singleton_evalue:
                singleton_presence[short].add(query_id)
                prev = singleton_evalue.get((short, query_id))
                if prev is None or evalue < prev:
                    singleton_evalue[(short, query_id)] = evalue

    # --- Build combined presence matrix ---
    # Collect all proteins: family members + singletons
    all_proteins = set()

    # Family members
    for rep in multi_reps:
        for member in rep_to_members[rep]:
            all_proteins.add(member)

    # Singletons (proteins in clusters of size 1)
    for rep in singleton_reps:
        all_proteins.add(rep)

    # Build presence: protein_id -> {proteome_short: 0/1}
    # Build evalues in parallel: protein_id -> {proteome_short: e-value or ''} (report-only
    # evidence; the source proteome's own cell is always '' since presence there isn't a hit).
    presence = {p: {s: 0 for s in all_shorts} for p in all_proteins}
    evalues = {p: {s: '' for s in all_shorts} for p in all_proteins}

    # Family members: present in proteomes where family HMM is present, plus source proteome
    for rep in multi_reps:
        family_present = set()
        for short, present_reps in family_presence.items():
            if rep in present_reps:
                family_present.add(short)
        for member in rep_to_members[rep]:
            # A fresh copy per member: each member is present in every proteome the family
            # HMM hit, plus (always) its own source proteome -- never another member's.
            present_proteomes = set(family_present)
            source = protein_to_proteome.get(member, '')
            if source:
                present_proteomes.add(source)
            if member in presence:
                for sp in present_proteomes:
                    if sp in presence[member]:
                        presence[member][sp] = 1
                        if sp != source:
                            evalues[member][sp] = family_evalue.get((sp, rep), '')

    # Singletons: present in proteomes where pairwise hit was found, plus source proteome
    for rep in singleton_reps:
        present_proteomes = set()
        for short, present_ids in singleton_presence.items():
            if rep in present_ids:
                present_proteomes.add(short)
        source = protein_to_proteome.get(rep, '')
        if source:
            present_proteomes.add(source)
        if rep in presence:
            for sp in present_proteomes:
                if sp in presence[rep]:
                    presence[rep][sp] = 1
                    if sp != source:
                        evalues[rep][sp] = singleton_evalue.get((sp, rep), '')

    # --- Write presence matrix ---
    with open(args.output_matrix, 'w') as out:
        out.write('protein_id\tsource_proteome\t' + '\t'.join(all_shorts) + '\n')
        for protein in sorted(all_proteins):
            source = protein_to_proteome.get(protein, '')
            vals = '\t'.join(str(presence[protein][s]) for s in all_shorts)
            out.write(f'{protein}\t{source}\t{vals}\n')

    if args.output_evalues:
        with open(args.output_evalues, 'w') as out:
            out.write('protein_id\tsource_proteome\t' + '\t'.join(all_shorts) + '\n')
            for protein in sorted(all_proteins):
                source = protein_to_proteome.get(protein, '')
                vals = '\t'.join(str(evalues[protein][s]) for s in all_shorts)
                out.write(f'{protein}\t{source}\t{vals}\n')

    # --- Filter to novelty candidates ---
    target_shorts = [s for s in all_shorts if short_to_group.get(s) == 'DISCOVERY_TARGET']
    disc_out_shorts = [s for s in all_shorts if short_to_group.get(s) == 'DISCOVERY_OUT']

    candidates = []
    for protein in sorted(all_proteins):
        target_present = sum(presence[protein][s] for s in target_shorts)
        disc_out_present = sum(presence[protein][s] for s in disc_out_shorts)

        target_frac = target_present / len(target_shorts) if target_shorts else 0
        disc_out_frac = disc_out_present / len(disc_out_shorts) if disc_out_shorts else 0

        if target_frac >= args.target_min_frac and disc_out_frac <= args.disc_out_max_frac:
            source = protein_to_proteome.get(protein, '')
            candidates.append(f'{source}::{protein}')

    with open(args.output_candidates, 'w') as out:
        if candidates:
            out.write('\n'.join(candidates) + '\n')

    print(f"Built presence matrix: {len(all_proteins)} proteins × {len(all_shorts)} proteomes; "
          f"{len(candidates)} novelty candidates", file=sys.stderr)


if __name__ == '__main__':
    main()
