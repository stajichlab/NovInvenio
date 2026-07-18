"""Gene-family membership derived from mmseqs easy-cluster output.

mmseqs easy-cluster's *_cluster.tsv lists one (rep_id, member_id) pair per
line, with the rep also listed as its own member.  build_families() collapses
that into {rep_id: [member_ids]}, keeping only clusters with more than one
member — a singleton candidate has no "subgroup" to report.
"""
from collections import defaultdict
from pathlib import Path


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


class FamilyIndex:
    """Maps protein_id -> family index from an mmseqs *_cluster.tsv, and
    accumulates per-family source-proteome sets across a single pass over
    matrix rows.

    Used by lib/report_data.py's payload builders to group candidates
    recovered independently in several species under one family rather than
    as unrelated rows.  Singleton clusters carry no family (index_of returns
    -1) — there is nothing to collapse.
    """

    def __init__(self, cluster_tsv):
        families = {}
        if cluster_tsv and Path(cluster_tsv).exists():
            families = build_families(read_cluster_tsv(cluster_tsv))
        self._families = families
        self._member_to_rep = {m: rep for rep, members in families.items() for m in members}
        self._reps = sorted(families.keys())
        self._rep_index = {rep: i for i, rep in enumerate(self._reps)}
        self._species: list[set] = [set() for _ in self._reps]

    def index_of(self, protein_id: str, source_proteome: str = '') -> int:
        """Return this protein's family index (-1 if not part of a family).

        Also records *source_proteome* against the family, so payload() can
        report how many distinct species independently recovered it.
        """
        rep = self._member_to_rep.get(protein_id)
        i = self._rep_index.get(rep, -1) if rep else -1
        if i >= 0 and source_proteome:
            self._species[i].add(source_proteome)
        return i

    def payload(self) -> list[dict]:
        return [
            {'rep': rep, 'size': len(self._families[rep]), 'species': sorted(self._species[i])}
            for i, rep in enumerate(self._reps)
        ]
