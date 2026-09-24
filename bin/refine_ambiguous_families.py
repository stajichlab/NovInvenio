#!/usr/bin/env python3
"""Tier R: split ambiguous mmseqs families using a targeted within-family diamond
search, so a small, cheap refinement pass can be measured against Tier C (raw
cluster membership) and Tier P (full pairwise) -- see notes/superpowers/specs/
2026-09-13-cluster-vs-pairwise-sensitivity-design.md for the full design and the
verified definition of "ambiguous family" this module implements.

A family is ambiguous iff:
  1. it is not in --oversized-families (excluded -- refinement cost is O(k^2) per
     family and a pathologically large cluster would dominate the whole run), AND
  2. at least one seed-group species contributes >=2 members (species-duplication --
     the only signal available before any HMM exists; confirmed exact on HEX1:
     Neurospora crassa contributes both HEX1_NEUCR and IF5A_NEUCR to one family), AND
  3. it is already a near-miss novelty candidate in the existing (Tier C+H)
     presence_matrix.tsv: seed-group presence fraction >= --ingroup-min-frac AND
     other-group presence fraction > --other-max-frac (i.e. currently rejected only
     because of other-group presence).

--query-group picks the seed group, as in build_presence_matrix.py: IN (default)
for the gain direction (families/, seeded from the ingroup), OUT for the loss
direction (loss_families/, seeded from the outgroup). With OUT the roles swap:
the outgroup is the seed group and the ingroup is the other group (issue #164).

Ambiguous families are refined by a single within-family diamond all-vs-all (over
the union of every ambiguous family's members, not one diamond call per family --
see main()), then split via find_paralog_edges_to_cut() + connected_components(),
using each member's own genome-local registered paralog (from the already-published
self_hits/<Short>.paralog_cutoffs.tsv, one per genome -- no new self-search) as the
forced split boundary. Non-ambiguous families pass through unmodified so the output
has the exact same rep<TAB>member / families.tsv-index contract as the pipeline's
own families_cluster.tsv + families.tsv, and score_controls.py's existing
--presence-mode cluster_membership path scores it unchanged.
"""
import argparse
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import INGROUP_ROLES, OUTGROUP_ROLES, parse_config  # noqa: E402
from hits import PARSERS, open_input  # noqa: E402

DEFAULT_DIAMOND_EVALUE = 1e-5


# --------------------------------------------------------------------- ambiguity
def seed_and_other_groups(samples, query_group):
    """(seed_ids, other_ids) for --query-group IN (gain) or OUT (loss)."""
    ingroup_ids = {s.short for s in samples if s.group in INGROUP_ROLES}
    outgroup_ids = {s.short for s in samples if s.group in OUTGROUP_ROLES}
    if query_group == 'IN':
        return ingroup_ids, outgroup_ids
    if query_group == 'OUT':
        return outgroup_ids, ingroup_ids
    raise ValueError(f'query_group must be IN or OUT, not {query_group!r}')


def species_of(protein_id, protein_to_proteome):
    return protein_to_proteome.get(protein_id)


def has_species_duplication(members, protein_to_proteome):
    species = [species_of(m, protein_to_proteome) for m in members]
    species = [s for s in species if s is not None]
    return len(species) != len(set(species))


def family_fractions(rep, members, presence_by_protein, ingroup_ids, outgroup_ids):
    """(ingroup_frac, outgroup_frac) present, from the run's ACTUAL Tier C+H presence
    (presence_matrix.tsv), looked up via the family representative's own row -- per the
    spec's ambiguous-family definition ("near-miss novelty candidate in the EXISTING
    presence_matrix.tsv"), NOT from raw cluster membership. Raw membership is always
    ingroup-only under this pipeline's clustering (ADR-0002), so an outgroup-fraction
    check against it can never be > 0 -- this was a real bug found running against real
    data (0/28623 families ever flagged ambiguous); ledgered 2026-09-13.
    """
    presence = presence_by_protein.get(rep, {})
    ing = sum(presence.get(s, 0) for s in ingroup_ids) / len(ingroup_ids) if ingroup_ids else 0.0
    out = sum(presence.get(s, 0) for s in outgroup_ids) / len(outgroup_ids) if outgroup_ids else 0.0
    return ing, out


def detect_ambiguous_families(fam_members, protein_to_proteome, presence_by_protein,
                              ingroup_ids, outgroup_ids, oversized_reps,
                              ingroup_min_frac, other_max_frac):
    ambiguous = set()
    for rep, members in fam_members.items():
        if rep in oversized_reps:
            continue
        if not has_species_duplication(members, protein_to_proteome):
            continue
        ing_frac, out_frac = family_fractions(rep, members, presence_by_protein,
                                              ingroup_ids, outgroup_ids)
        if ing_frac >= ingroup_min_frac and out_frac > other_max_frac:
            ambiguous.add(rep)
    return ambiguous


def load_profiled_reps(families_tsv):
    """representative_id set for families that were actually profiled (>= min members) --
    i.e. the pipeline's own families.tsv, which is narrower than raw --cluster-tsv
    membership. Mirrors score_controls.py's load_profiled_reps()."""
    reps = set()
    with open(families_tsv) as fh:
        fh.readline()  # header: family_index, representative_id, n_members
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[1]:
                reps.add(parts[1])
    return reps


def filter_to_profiled(fam_members, profiled_reps):
    """Restrict fam_members (rep -> [member, ...]) to reps that were actually profiled
    by the real pipeline (in profiled_reps) -- without this, Tier R's family universe
    silently includes unprofiled singletons/small clusters the real pipeline excluded,
    a real bug found running against real data (POS_LAH's singleton cluster was
    incorrectly resurrected this way)."""
    return {rep: members for rep, members in fam_members.items() if rep in profiled_reps}


# ---------------------------------------------------------------- graph splitting
def load_paralog_map(cutoff_files):
    """protein_ID -> paralog_protein_ID, from one or more per-genome
    self_hits/<Short>.paralog_cutoffs.tsv files (already published by Tier P --
    see parse_self_hits.py; no new self-search)."""
    paralog_of = {}
    for path in cutoff_files:
        with open(path) as fh:
            header = fh.readline()
            del header
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    paralog_of[parts[0]] = parts[1]
    return paralog_of


def build_within_family_edges(hit_lines, evalue_cutoff):
    """{frozenset({a, b}), ...} for every diamond hit pair passing evalue_cutoff."""
    edges = set()
    for hit in hit_lines:
        if hit.query_id == hit.target_id:
            continue
        if hit.evalue < evalue_cutoff:
            edges.add(frozenset({hit.query_id, hit.target_id}))
    return edges


def find_paralog_edges_to_cut(members, paralog_map):
    """Edges to remove before connected components: X-P where P is X's own
    registered within-genome paralog and P is a co-member of the same family --
    mirrors build_presence_matrix.py's target-scope paralog-competition test
    (the one that already rescues HEX-1 in Tier P), applied as a forced graph
    split instead of a per-hit disqualification."""
    member_set = set(members)
    cuts = set()
    for m in members:
        p = paralog_map.get(m)
        if p and p in member_set and p != m:
            cuts.add(frozenset({m, p}))
    return cuts


def connected_components(members, edges):
    """List of member-id lists, one per connected component (a lone member with no
    surviving edge becomes its own singleton component)."""
    parent = {m: m for m in members}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for edge in edges:
        a, b = tuple(edge)
        if a in parent and b in parent:
            union(a, b)

    groups = defaultdict(list)
    for m in members:
        groups[find(m)].append(m)
    return list(groups.values())


def split_family(members, edges, paralog_map):
    """members -> list of subfamily member-lists, after removing forced-split
    (registered-paralog) edges and taking connected components on what remains."""
    cuts = find_paralog_edges_to_cut(members, paralog_map)
    remaining = edges - cuts
    return connected_components(members, remaining)


# ------------------------------------------------------------------- sequences
def extract_needed_sequences(pep_paths, needed_ids):
    """protein_id -> sequence, scanning each --Short.pep.fa once, keeping only ids
    in needed_ids (the union of every ambiguous family's members) -- avoids loading
    whole proteomes when only a small fraction of their proteins are needed."""
    seqs = {}
    remaining = set(needed_ids)
    for path in pep_paths.values():
        if not remaining:
            break
        current_id, chunks = None, []
        with open(path) as fh:
            for line in fh:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if current_id in remaining:
                        seqs[current_id] = ''.join(chunks)
                        remaining.discard(current_id)
                    current_id = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line)
            if current_id in remaining:
                seqs[current_id] = ''.join(chunks)
                remaining.discard(current_id)
    return seqs


def write_fasta(seqs, path):
    with open(path, 'w') as fh:
        for pid, seq in seqs.items():
            fh.write(f'>{pid}\n{seq}\n')


def run_diamond_within_family(fasta_path, cpus=1):
    """One diamond makedb + blastp all-vs-all over the combined ambiguous-family
    FASTA (not per-family -- process-startup overhead would dominate at ~1000+
    families). Returns the path to the raw outfmt-6 hits file."""
    db_path = fasta_path.with_suffix('.dmnd')
    subprocess.run(['diamond', 'makedb', '--in', str(fasta_path), '--db', str(db_path)],
                    check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    hits_path = fasta_path.with_suffix('.hits.tsv')
    subprocess.run([
        'diamond', 'blastp', '--very-sensitive', '--threads', str(cpus),
        '--query', str(fasta_path), '--db', str(db_path),
        '--outfmt', '6', 'qseqid', 'sseqid', 'evalue', 'bitscore',
        '--out', str(hits_path),
    ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return hits_path


# --------------------------------------------------------------------- output
def write_refined_families(fam_members, ambiguous_reps, subfamilies_by_rep,
                           out_cluster_tsv, out_families_tsv):
    """Same rep<TAB>member / families.tsv-index shape as the pipeline's own
    families_cluster.tsv + families.tsv, so score_controls.py's existing
    --cluster-tsv/--families loaders consume this file unmodified. Non-ambiguous
    families pass through verbatim; ambiguous ones are replaced by their split
    subfamilies, each keyed by its own first member as the new representative."""
    with open(out_cluster_tsv, 'w') as cfh, open(out_families_tsv, 'w') as ffh:
        ffh.write('family_index\trepresentative_id\tn_members\n')
        idx = 0
        for rep, members in fam_members.items():
            if rep not in ambiguous_reps:
                idx += 1
                ffh.write(f'fam_{idx:06d}\t{rep}\t{len(members)}\n')
                for m in members:
                    cfh.write(f'{rep}\t{m}\n')
                continue
            for subfamily in subfamilies_by_rep.get(rep, [members]):
                new_rep = subfamily[0]
                idx += 1
                ffh.write(f'fam_{idx:06d}\t{new_rep}\t{len(subfamily)}\n')
                for m in subfamily:
                    cfh.write(f'{new_rep}\t{m}\n')


# ------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv',
                    help="run's families_cluster.tsv (rep<TAB>member)")
    ap.add_argument('--families', required=True,
                    help="run's families.tsv (profiled families index)")
    ap.add_argument('--matrix', required=True,
                    help="run's presence_matrix.tsv (for protein_id -> source_proteome "
                         "and the near-miss-novelty check)")
    ap.add_argument('--oversized-families', default=None, dest='oversized_families',
                    help="run's oversized_families.tsv (excluded from refinement)")
    ap.add_argument('--config', required=True, help='Analysis description CSV')
    ap.add_argument('--query-group', choices=['IN', 'OUT'], default='IN', dest='query_group',
                    help='Seed group of --cluster-tsv: IN (default, gain direction, '
                         'families/) or OUT (loss direction, loss_families/)')
    ap.add_argument('--pep', nargs='+', required=True,
                    help='Short=path.pep.fa pairs for every seed-group proteome')
    ap.add_argument('--self-hits', nargs='+', required=True, dest='self_hits',
                    help='self_hits/<Short>.paralog_cutoffs.tsv files (already '
                         'published by the pairwise run for these same genomes)')
    ap.add_argument('--ingroup-min-frac', type=float, default=0.75, dest='ingroup_min_frac',
                    help='Min seed-group presence fraction (the ingroup for '
                         '--query-group IN, the outgroup for OUT)')
    ap.add_argument('--other-max-frac', type=float, default=0.0, dest='other_max_frac')
    ap.add_argument('--diamond-evalue', type=float, default=DEFAULT_DIAMOND_EVALUE,
                    dest='diamond_evalue')
    ap.add_argument('--cpus', type=int, default=1)
    ap.add_argument('--tmp-dir', default=None, dest='tmp_dir')
    ap.add_argument('--output-cluster-tsv', required=True, dest='output_cluster_tsv')
    ap.add_argument('--output-families', required=True, dest='output_families')
    args = ap.parse_args()

    import pandas as pd  # local import: only main() needs it, unlike the pure fns above

    samples = parse_config(args.config)
    seed_ids, other_ids = seed_and_other_groups(samples, args.query_group)
    pep_paths = {}
    for pair in args.pep:
        short, path = pair.split('=', 1)
        pep_paths[short] = Path(path)

    matrix = pd.read_csv(args.matrix, sep='\t')
    protein_to_proteome = dict(zip(matrix['protein_id'], matrix['source_proteome']))
    proteome_cols = [c for c in matrix.columns if c not in ('protein_id', 'source_proteome')]
    presence_by_protein = matrix.set_index('protein_id')[proteome_cols].astype(int).to_dict('index')

    fam_members = defaultdict(list)
    with open(args.cluster_tsv) as fh:
        for line in fh:
            rep, member = line.rstrip('\n').split('\t')[:2]
            fam_members[rep].append(member)
    fam_members = dict(fam_members)

    profiled_reps = load_profiled_reps(args.families)
    fam_members = filter_to_profiled(fam_members, profiled_reps)

    oversized_reps = set()
    if args.oversized_families:
        with open(args.oversized_families) as fh:
            header = fh.readline()
            del header
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if parts and parts[0]:
                    oversized_reps.add(parts[0])

    ambiguous = detect_ambiguous_families(
        fam_members, protein_to_proteome, presence_by_protein, seed_ids, other_ids,
        oversized_reps, args.ingroup_min_frac, args.other_max_frac)
    print(f'{len(ambiguous)} / {len(fam_members)} families flagged ambiguous', file=sys.stderr)

    paralog_map = load_paralog_map(args.self_hits)

    tmp_dir = Path(args.tmp_dir) if args.tmp_dir else Path(tempfile.mkdtemp())
    tmp_dir.mkdir(parents=True, exist_ok=True)
    needed_ids = {m for rep in ambiguous for m in fam_members[rep]}
    seqs = extract_needed_sequences(pep_paths, needed_ids)
    combined_fasta = tmp_dir / 'ambiguous_family_members.fa'
    write_fasta(seqs, combined_fasta)

    subfamilies_by_rep = {}
    if ambiguous:
        hits_path = run_diamond_within_family(combined_fasta, cpus=args.cpus)
        with open_input(hits_path) as fh:
            all_edges = build_within_family_edges(PARSERS['diamond'](fh), args.diamond_evalue)
        for rep in ambiguous:
            members = fam_members[rep]
            member_set = set(members)
            family_edges = {e for e in all_edges if e <= member_set}
            subfamilies_by_rep[rep] = split_family(members, family_edges, paralog_map)

    write_refined_families(fam_members, ambiguous, subfamilies_by_rep,
                           args.output_cluster_tsv, args.output_families)


if __name__ == '__main__':
    main()
