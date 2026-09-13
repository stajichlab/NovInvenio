#!/usr/bin/env python3
"""Tier R: split ambiguous mmseqs families using a targeted within-family diamond
search, so a small, cheap refinement pass can be measured against Tier C (raw
cluster membership) and Tier P (full pairwise) -- see notes/superpowers/specs/
2026-09-13-cluster-vs-pairwise-sensitivity-design.md for the full design and the
verified definition of "ambiguous family" this module implements.

A family is ambiguous iff:
  1. it is not in --oversized-families (excluded -- refinement cost is O(k^2) per
     family and a pathologically large cluster would dominate the whole run), AND
  2. at least one ingroup species contributes >=2 members (species-duplication --
     the only signal available before any HMM exists; confirmed exact on HEX1:
     Neurospora crassa contributes both HEX1_NEUCR and IF5A_NEUCR to one family), AND
  3. it is already a near-miss novelty candidate in the existing (Tier C+H)
     presence_matrix.tsv: ingroup presence fraction >= --ingroup-min-frac AND
     outgroup presence fraction > --other-max-frac (i.e. currently rejected only
     because of outgroup presence).

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
def species_of(protein_id, protein_to_proteome):
    return protein_to_proteome.get(protein_id)


def has_species_duplication(members, protein_to_proteome):
    species = [species_of(m, protein_to_proteome) for m in members]
    species = [s for s in species if s is not None]
    return len(species) != len(set(species))


def family_fractions(rep, members, protein_to_proteome, ingroup_ids, outgroup_ids):
    """(ingroup_frac, outgroup_frac) present, from raw cluster membership -- mirrors
    build_cluster_membership_presence()'s notion of presence (member's own genome
    counts as present), i.e. this is Tier C's own presence call for this family,
    reused here purely to find near-miss novelty candidates worth refining."""
    proteomes = {species_of(m, protein_to_proteome) for m in members}
    proteomes.discard(None)
    ing = len(proteomes & ingroup_ids) / len(ingroup_ids) if ingroup_ids else 0.0
    out = len(proteomes & outgroup_ids) / len(outgroup_ids) if outgroup_ids else 0.0
    return ing, out


def detect_ambiguous_families(fam_members, protein_to_proteome, ingroup_ids, outgroup_ids,
                              oversized_reps, ingroup_min_frac, other_max_frac):
    ambiguous = set()
    for rep, members in fam_members.items():
        if rep in oversized_reps:
            continue
        if not has_species_duplication(members, protein_to_proteome):
            continue
        ing_frac, out_frac = family_fractions(rep, members, protein_to_proteome,
                                              ingroup_ids, outgroup_ids)
        if ing_frac >= ingroup_min_frac and out_frac > other_max_frac:
            ambiguous.add(rep)
    return ambiguous


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


if __name__ == '__main__':
    pass  # main() wired in Task 4
