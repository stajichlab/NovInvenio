"""Shared paralog-aware singleton hit scoring (issue #52).

Applies the same filter logic bin/build_presence_matrix.py uses for the classic pairwise
pathway to a singleton pairwise search's hits:

  1. Significance filter: hit e-value < default_evalue (flat, applied to every hit).
  2. Paralog-competition filter: disqualify a hit the query's paralog beats -- scope is
     'proteome' (the paralog's best hit anywhere in the target proteome) or 'target'
     (only on the same target protein), matching build_presence_matrix.py's
     --paralog-competition-scope semantics.

Used by bin/novelty_presence_matrix.py (phase 1, vs DISCOVERY_OUT) and
bin/novelty_screen.py (phase 2, vs NEAR_INGROUP/BROAD_OUTGROUP) so a singleton's
definition of "present" stays identical across both phases -- filter 2 is what catches
cross-reactivity with a conserved paralog (e.g. NCU08332/HEX-1 vs eIF5A) that a flat
e-value threshold alone can't tell apart from a real ortholog.

(2026-09-03: filter 1 used to be a per-query "paralog-cutoff" -- hit e-value must beat the
query's own within-proteome paralog e-value (from a self-vs-self search reported down to
E=100, not the usual significance floor), falling back to default_evalue when no paralog
was detected. Dropped: this swung both too tight -- any protein with a close in-genome
paralog got an unreachable cutoff (measured on real N. crassa self-search data: 28-48% of
proteins had a "cutoff" tighter than 1e-50, rejecting real orthologs regardless of whether
the paralog was actually relevant to that hit -- e.g. the reported NCU00765/NCU00411/
NCU01935 false-absence findings) -- and too loose -- 12-28% of proteins with no real
paralog still picked up self-search noise (E up to 100) as their "cutoff", looser than the
intended default. Filter 2 already does the actual paralogy test as a direct head-to-head
comparison instead of an absolute-magnitude proxy, so filter 1 was redundant on top of it,
not protective.)
"""
import csv
from collections import defaultdict
from pathlib import Path


def parse_pairwise_tsv(path):
    """Return list of (query_id, target_id, evalue) tuples from a parsed hits TSV.

    Only the first three columns (query_id, target_id, evalue) are used; bitscore/
    query_proteome/target_proteome, if present, are ignored.
    """
    hits = []
    with open(path) as fh:
        next(fh, None)  # header
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


def proteome_short_from_hits_filename(path):
    """Derive the target proteome short from a singleton-search hits filename.

    PARSE_HITS names its output '<query_id>_vs_<target_id>.parsed.tsv' (the same
    convention workflows/search.nf uses) -- a singleton/context search's query_id is
    always a fixed literal (e.g. 'singletons'), so the proteome short is whatever
    follows the last '_vs_'.
    """
    name = Path(path).name
    for suffix in ('.parsed.tsv', '.tsv', '.parsed'):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    return name.rsplit('_vs_', 1)[1] if '_vs_' in name else name


def load_paralog_info(cutoff_files):
    """Return paralog_of: protein_id -> paralog_protein_id (or None if not detected), from
    paralog_cutoffs.tsv files -- feeds filter 2 (paralog-competition) only; see module
    docstring for why the per-query e-value column is no longer used as a significance
    cutoff (filter 1)."""
    paralog_of = {}
    for path in cutoff_files:
        with open(path, newline='') as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                pid = row.get('protein_ID')
                if not pid:
                    continue
                paralog_of[pid] = row.get('paralog_protein_ID') or None
    return paralog_of


def score_singleton_hits(hits, singleton_ids, paralog_of,
                         default_evalue, competition_scope='proteome'):
    """Filter a singleton search's hits down to qualifying presence calls.

    hits: iterable of (query_id, target_id, evalue, proteome_short) tuples -- covers
    both true singletons and, where detected, their own within-genome paralog (searched
    together so filter 2 has something to compare against; see
    bin/extend_singleton_query.py).

    singleton_ids: set of query_ids that are true singletons. Only these are ever
    scored/returned -- a paralog present in `hits` purely to support the comparison is
    never itself emitted unless it independently is also in singleton_ids.

    default_evalue: flat significance cutoff (filter 1), applied to every hit regardless
    of query -- see module docstring for why this is no longer per-query paralog-derived.

    Returns (presence, evalue):
      presence[proteome_short] = set of singleton_ids present
      evalue[(proteome_short, singleton_id)] = best (lowest) qualifying hit e-value
    """
    hits = list(hits)

    best_ev = {}
    for query_id, target_id, ev, short in hits:
        key = (query_id, short if competition_scope == 'proteome' else target_id)
        if key not in best_ev or ev < best_ev[key]:
            best_ev[key] = ev

    presence = defaultdict(set)
    evalue = {}
    for query_id, target_id, ev, short in hits:
        if query_id not in singleton_ids:
            continue

        if not (ev < default_evalue):
            continue

        paralog_id = paralog_of.get(query_id)
        if paralog_id:
            key = (paralog_id, short if competition_scope == 'proteome' else target_id)
            paralog_ev = best_ev.get(key)
            if paralog_ev is not None and paralog_ev < ev:
                continue  # disqualified: the paralog explains this hit better

        presence[short].add(query_id)
        prev = evalue.get((short, query_id))
        if prev is None or ev < prev:
            evalue[(short, query_id)] = ev

    return presence, evalue
