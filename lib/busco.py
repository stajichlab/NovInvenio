"""Shared BUSCO full_table.tsv / cluster-membership loading for the family-profile sweep
(ADR-0002 Q8, issue #6). Used by bin/busco_family_recovery.py (clustering-quality metric)
and bin/busco_presence_recovery.py (presence-calling-quality metric, 2026-09-03) so both
read BUSCO input the same way.
"""
import csv
from collections import defaultdict


def parse_busco_full_table(path, species):
    """Yield (busco_id, species, protein_id, length) for Complete rows of a BUSCO
    full_table.tsv. `length` is the BUSCO marker's alignment length (column 5, `Length`)
    -- a readily-available proxy for the gene's size, used to check whether a metric
    computed from Complete BUSCOs is representative of the whole proteome's length
    distribution (a real risk: BUSCO markers are hand-picked conserved single-copy
    orthologs, not a random sample of genes) rather than assumed representative.
    """
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 5:
                continue
            busco_id, status, sequence = f[0], f[1], f[2]
            if status != 'Complete':
                continue
            # BUSCO may append :start-end to the sequence id; strip to the bare protein id.
            protein_id = sequence.split(':')[0]
            if not protein_id:
                continue
            try:
                length = int(f[4])
            except ValueError:
                length = None
            yield busco_id, species, protein_id, length


def load_busco_map(map_path=None, tables=None):
    """Return (mapping, lengths).

    mapping: busco_id -> {species: protein_id}, for Complete single-copy BUSCOs.
    lengths: busco_id -> alignment length (int), only populated from --tables (raw
    full_table.tsv rows carry it; a consolidated --busco-map may not).

    From a consolidated --busco-map (busco_id<TAB>species<TAB>protein_id[...]) or from raw
    per-species full_table.tsv files given as SHORT=path.
    """
    mapping: dict[str, dict[str, str]] = defaultdict(dict)
    lengths: dict[str, int] = {}
    if map_path:
        with open(map_path) as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) < 3 or f[0].startswith('#'):
                    continue
                # a status column, if present, must be Complete
                if len(f) >= 4 and f[3] and f[3] != 'Complete':
                    continue
                mapping[f[0]][f[1]] = f[2]
    for spec in tables or []:
        if '=' not in spec:
            raise SystemExit(f'--tables entries must be SHORT=path (got: {spec!r})')
        short, path = spec.split('=', 1)
        for busco_id, species, protein_id, length in parse_busco_full_table(path, short):
            mapping[busco_id][species] = protein_id
            if length is not None:
                lengths.setdefault(busco_id, length)
    return mapping, lengths


def load_profiled_reps(families_tsv):
    """representative_id set for families that were actually profiled (>= min members)."""
    reps = set()
    with open(families_tsv) as fh:
        fh.readline()  # header: family_index, representative_id, n_members
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[1]:
                reps.add(parts[1])
    return reps


def load_member_to_rep(cluster_tsv, keep_reps):
    """member_id -> rep_id, restricted to families in keep_reps (i.e. profiled)."""
    member_to_rep = {}
    with open(cluster_tsv) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            if rep in keep_reps:
                member_to_rep[member] = rep
    return member_to_rep


def load_rep_to_members(cluster_tsv, keep_reps):
    """rep_id -> list of member_ids, restricted to families in keep_reps."""
    rep_to_members: dict[str, list] = defaultdict(list)
    with open(cluster_tsv) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            if rep in keep_reps:
                rep_to_members[rep].append(member)
    return rep_to_members


def read_matrix_rows(matrix_tsv):
    """Return dict: protein_id -> {proteome_short: '0'/'1', ...} from a presence_matrix.tsv
    (columns: protein_id, source_proteome, <sorted proteome shorts>)."""
    rows = {}
    with open(matrix_tsv, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        proteome_cols = [c for c in reader.fieldnames or []
                         if c not in ('protein_id', 'source_proteome')]
        for row in reader:
            rows[row['protein_id']] = {c: row.get(c, '0') for c in proteome_cols}
    return rows
