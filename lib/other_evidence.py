"""Other-group protein-hit evidence record (issue #208).

Builds the long-format sidecar `*.other_evidence.tsv.gz`: one row per (query protein,
other-group proteome) cell that has at least one hit below the significance floor,
including hits that filter 2 (paralog competition) or filter 3 (coverage floor) removed.
Evidence only. Nothing here changes a presence call.

Shared by bin/build_presence_matrix.py and lib/singleton_presence.py so both pathways
write the same columns.
"""
import math

import pandas as pd

EVIDENCE_COLUMNS = [
    'protein_id', 'source_proteome', 'other_proteome', 'evidence',
    'best_target_id', 'best_evalue', 'best_bitscore', 'best_qcov',
    'max_qcov', 'n_hits', 'status',
    'paralog_id', 'paralog_evalue', 'paralog_delta',
]

# Order of precedence when no hit in a cell is kept: the status of the cell's best hit.
STATUS_KEPT = 'kept'


def _delta(best_ev, paralog_ev):
    """log10(best) - log10(paralog); +inf when the paralog E-value is 0; '' if missing."""
    if best_ev is None or paralog_ev is None or pd.isna(best_ev) or pd.isna(paralog_ev):
        return ''
    if paralog_ev <= 0:
        return math.inf
    if best_ev <= 0:
        return -math.inf
    return math.log10(best_ev) - math.log10(paralog_ev)


def build_other_evidence(hits, status, other_ids, paralog_of, best_ev, scope):
    """Return the evidence DataFrame (EVIDENCE_COLUMNS).

    hits       DataFrame of significant hits (after filter 1) with columns query_id,
               query_proteome, target_id, target_proteome, evalue, bitscore and,
               optionally, qcov.
    status     Series aligned to hits.index: 'kept', 'paralog_filtered' or 'coverage_floor'.
    other_ids  proteomes judged as the other group.
    paralog_of protein_id -> paralog protein id (or None).
    best_ev    {(query_proteome, query_id, key): best evalue}, key = target_id for scope
               'target' or target_proteome for scope 'proteome' (as build_presence_matrix).
    """
    if hits.empty:
        return pd.DataFrame(columns=EVIDENCE_COLUMNS)
    h = hits[hits['target_proteome'].isin(other_ids)].copy()
    if h.empty:
        return pd.DataFrame(columns=EVIDENCE_COLUMNS)
    h['status'] = status.reindex(h.index)
    h['qcov_n'] = (pd.to_numeric(h['qcov'], errors='coerce') if 'qcov' in h.columns
                   else pd.Series(float('nan'), index=h.index))
    rows = []
    for (qp, qid, tp), g in h.groupby(['query_proteome', 'query_id', 'target_proteome'],
                                      sort=True):
        best = g.loc[g['evalue'].idxmin()]
        cell_status = STATUS_KEPT if (g['status'] == STATUS_KEPT).any() else best['status']
        pid = paralog_of.get(qid)
        pid = None if pid is None or pd.isna(pid) else pid
        pev = None
        if pid is not None:
            key = best['target_id'] if scope == 'target' else tp
            pev = best_ev.get((qp, pid, key))
        max_q = g['qcov_n'].max()
        rows.append({
            'protein_id': qid, 'source_proteome': qp, 'other_proteome': tp,
            'evidence': 'protein',
            'best_target_id': best['target_id'],
            'best_evalue': best['evalue'], 'best_bitscore': best['bitscore'],
            'best_qcov': '' if pd.isna(best['qcov_n']) else best['qcov_n'],
            'max_qcov': '' if pd.isna(max_q) else max_q,
            'n_hits': len(g), 'status': cell_status,
            'paralog_id': pid or '',
            'paralog_evalue': '' if pev is None else pev,
            'paralog_delta': _delta(best['evalue'], pev),
        })
    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)
