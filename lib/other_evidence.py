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


# --- per-candidate summary (bin/make_novelties.py) --------------------------

NOVELTY_COLUMNS = ['other_protein_max_qcov', 'other_protein_n_filtered',
                   'other_protein_min_paralog_evalue', 'other_tblastn_max_cov',
                   'other_tblastn_n_genomes', 'other_signal_class']


def _num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def load_protein_evidence(path):
    """Return {protein_id: [row dict, ...]} from a *.other_evidence.tsv(.gz).

    None when the file is missing or 0 bytes (a stub): the evidence was not measured. A
    header-only file returns {} (measured, nothing found).
    """
    import csv
    import gzip
    import os
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    opener = gzip.open if str(path).endswith('.gz') else open
    out = {}
    with opener(path, 'rt', newline='') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            out.setdefault(row['protein_id'], []).append(row)
    return out


def load_tblastn_coverage(path):
    """Return {protein_id: [row dict, ...]} from tblastn_summary.coverage.tsv(.gz); None if a stub."""
    return load_protein_evidence(path)


def signal_class(c, threshold, has_hit):
    """'none' (no hit), 'domain_only' (C < T), 'broad' (C >= T); '' when a hit has no coverage."""
    if not has_hit:
        return 'none'
    if c is None:
        return ''
    return 'broad' if c >= threshold else 'domain_only'


def candidate_signal(protein_rows, tblastn_rows, threshold):
    """Per-candidate NOVELTY_COLUMNS values (strings) from its evidence rows.

    protein_rows / tblastn_rows may be None when that sidecar was not supplied ("not
    measured": its columns stay blank). C = max(protein max_qcov, tblastn coverage x 100)
    over the evidence that exists.
    """
    out = {k: '' for k in NOVELTY_COLUMNS}
    qcovs, has_hit = [], False
    if protein_rows is not None:
        qs = [_num(r.get('max_qcov')) for r in protein_rows]
        qs = [q for q in qs if q is not None]
        out['other_protein_max_qcov'] = f'{max(qs):g}' if qs else ''
        filtered = [r for r in protein_rows if r.get('status') != STATUS_KEPT]
        out['other_protein_n_filtered'] = str(len(filtered))
        pe = [_num(r.get('paralog_evalue')) for r in filtered]
        pe = [p for p in pe if p is not None]
        out['other_protein_min_paralog_evalue'] = f'{min(pe):g}' if pe else ''
        qcovs += qs
        has_hit |= bool(protein_rows)
    if tblastn_rows is not None:
        cs = [_num(r.get('query_span_cov')) for r in tblastn_rows]
        cs = [c for c in cs if c is not None]
        out['other_tblastn_max_cov'] = f'{max(cs):.3f}' if cs else ''
        out['other_tblastn_n_genomes'] = str(len({r['genome'] for r in tblastn_rows}))
        qcovs += [c * 100 for c in cs]
        has_hit |= bool(tblastn_rows)
    if protein_rows is None and tblastn_rows is None:
        return out
    out['other_signal_class'] = signal_class(max(qcovs) if qcovs else None, threshold, has_hit)
    return out
