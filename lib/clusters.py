"""Gene-family membership derived from mmseqs easy-cluster output.

mmseqs easy-cluster's *_cluster.tsv lists one (rep_id, member_id) pair per
line, with the rep also listed as its own member.  build_families() collapses
that into {rep_id: [member_ids]}, keeping only clusters with more than one
member — a singleton candidate has no "subgroup" to report.
"""
from collections import defaultdict


def read_cluster_tsv(path):
    """Return {member_protein_id: rep_protein_id} from an mmseqs *_cluster.tsv."""
    member_to_rep = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            member_to_rep[member] = rep
    return member_to_rep


def build_families(member_to_rep):
    """Return {rep_id: sorted [member_ids]} for clusters with >1 member."""
    families = defaultdict(list)
    for member, rep in member_to_rep.items():
        families[rep].append(member)
    return {rep: sorted(members) for rep, members in families.items() if len(members) > 1}
