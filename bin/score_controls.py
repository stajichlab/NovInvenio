#!/usr/bin/env python3
"""Score a family-profile run against curated biological controls (ADR-0002 Q8).

Reads a controls CSV (`configs/controls/<clade>.controls.csv`; see
`configs/controls/README.md`), resolves each control's anchor to a **gene family** in a
given run, reads that family's novelty call from the run's `presence_matrix.tsv`, and
compares it to the control's `expected_call`:

  * positive controls (expected_call = novel) measure **recall / sensitivity** — a true
    lineage-specific gene that *should* be flagged.
  * negative controls (expected_call = core) measure the **false-novelty / FP rate** — a
    conserved gene that must **not** be flagged (BUSCO single-copy orthologs are free
    negatives).

One recall + FP number per run, so the parameter sweep (issue #6) can score each grid
point on biology, not just candidate counts.

Anchor resolution (`anchor_type` column):
  * protein_id — direct: the protein's gene family (via the mmseqs cluster membership).
  * fasta      — a sequence under configs/controls/seqs/: hmmsearch the family HMM db
                 against it; the best-scoring family is the anchor's family (survives
                 re-annotation / ID drift).
  * busco      — a BUSCO id resolved to a protein via --busco-map (busco_id<TAB>protein_id),
                 then treated like protein_id. Without --busco-map, busco rows are reported
                 unresolved (the BUSCO recovery machinery lands in issue #6).

The "call" is the same novelty predicate the pathway itself applies (profile_to_matrix.py
keep-rule): a family is called **novel** iff its presence fraction within the query group
(ingroup for novelty) is >= --ingroup-min-frac AND its presence fraction in the other
group is <= --other-max-frac. Anything else is scored as **not novel** (a negative
control's expected "core").

Placeholder rows (control_id starting `EXAMPLE_`, or an anchor still wrapped in <...>) are
skipped so the shipped template is inert until curated.
"""
import argparse
import csv
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import parse_config  # noqa: E402

META_COLS = ('protein_id', 'source_proteome')


# --------------------------------------------------------------------------- controls
def is_placeholder(row: dict) -> bool:
    """A template placeholder row that should be skipped (not yet curated)."""
    cid = (row.get('control_id') or '').strip()
    anchor = (row.get('anchor') or '').strip()
    if cid.startswith('EXAMPLE_'):
        return True
    if '<' in anchor or '>' in anchor:
        return True
    return False


def load_controls(path):
    """Yield curated control rows (dicts), skipping placeholders."""
    kept, skipped = [], 0
    with open(path) as fh:
        for row in csv.DictReader(fh):
            if is_placeholder(row):
                skipped += 1
                continue
            kept.append(row)
    return kept, skipped


# ---------------------------------------------------------------------- family lookup
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


def load_family_membership(cluster_tsv, keep_reps):
    """Return (member_to_rep, rep_to_members) restricted to profiled families."""
    member_to_rep = {}
    rep_to_members = defaultdict(list)
    with open(cluster_tsv) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            if rep in keep_reps:
                member_to_rep[member] = rep
                rep_to_members[rep].append(member)
    return member_to_rep, rep_to_members


def family_presence_vector(matrix, members, proteome_cols):
    """OR the presence rows of a family's members → {proteome_short: 0/1}.

    All members share the family HMM's per-outgroup presence; each member row also marks
    its own source proteome present. The union across member rows is the family's full
    presence vector.
    """
    sub = matrix[matrix['protein_id'].isin(members)]
    if sub.empty:
        return None
    present = (sub[proteome_cols].max(axis=0) > 0).astype(int)
    return present.to_dict()


def family_call(presence, ingroup_ids, outgroup_ids, ingroup_min_frac, other_max_frac):
    """'novel' or 'not_novel' for a family's presence vector (novelty predicate)."""
    if not ingroup_ids:
        return 'not_novel'
    q = sum(presence.get(s, 0) for s in ingroup_ids) / len(ingroup_ids)
    o = (sum(presence.get(s, 0) for s in outgroup_ids) / len(outgroup_ids)
         if outgroup_ids else 0.0)
    novel = (q >= ingroup_min_frac) and (o <= other_max_frac)
    return 'novel' if novel else 'not_novel'


# --------------------------------------------------------------------- anchor resolve
def resolve_protein_anchor(protein_id, member_to_rep):
    """protein_id → family rep, or None if the protein is in no profiled family."""
    return member_to_rep.get(protein_id)


def best_family_from_domtblout(path):
    """Best (lowest full-seq E-value) family rep in an hmmsearch --domtblout, or None.

    Column 4 (0-indexed 3) is the query = family HMM name = family representative id
    (hmmbuild -n <rep>). Column 7 (idx 6) is the full-sequence E-value.
    """
    best_rep, best_e = None, None
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            f = line.split()
            if len(f) < 7:
                continue
            query = f[3]
            try:
                full_e = float(f[6])
            except ValueError:
                continue
            if best_e is None or full_e < best_e:
                best_rep, best_e = query, full_e
    return best_rep


def resolve_fasta_anchor(fasta_path, profiles_hmm, cpus=1):
    """hmmsearch the family HMM db against a sequence anchor → its best family rep."""
    if not profiles_hmm:
        return None
    with tempfile.NamedTemporaryFile(suffix='.domtblout', delete=False) as tf:
        dom = tf.name
    try:
        subprocess.run(
            ['hmmsearch', '--cpu', str(cpus), '--domtblout', dom, profiles_hmm, fasta_path],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        return best_family_from_domtblout(dom)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    finally:
        Path(dom).unlink(missing_ok=True)


def load_busco_map(path):
    """busco_id → protein_id, from a two-column TSV (no header)."""
    mapping = {}
    if not path:
        return mapping
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[0]:
                mapping[parts[0]] = parts[1]
    return mapping


# ---------------------------------------------------------------------------- scoring
def resolve_anchor(row, member_to_rep, busco_map, controls_dir, profiles_hmm, cpus):
    """Resolve a control row's anchor to a family rep (or None). Returns (rep, note)."""
    atype = (row.get('anchor_type') or '').strip()
    anchor = (row.get('anchor') or '').strip()
    if atype == 'protein_id':
        rep = resolve_protein_anchor(anchor, member_to_rep)
        return rep, ('' if rep else 'protein not in any profiled family')
    if atype == 'busco':
        pid = busco_map.get(anchor)
        if not pid:
            return None, 'no --busco-map entry (see issue #6)'
        rep = resolve_protein_anchor(pid, member_to_rep)
        return rep, ('' if rep else 'busco protein not in any profiled family')
    if atype == 'fasta':
        fpath = anchor if Path(anchor).is_absolute() else str(controls_dir / anchor)
        if not Path(fpath).exists():
            return None, f'fasta anchor not found: {fpath}'
        rep = resolve_fasta_anchor(fpath, profiles_hmm, cpus)
        return rep, ('' if rep else 'no family hmmsearch hit (or hmmsearch unavailable)')
    return None, f'unknown anchor_type: {atype!r}'


def score_controls(controls, matrix, member_to_rep, rep_to_members, samples,
                   ingroup_min_frac, other_max_frac, busco_map, controls_dir,
                   profiles_hmm, cpus):
    proteome_cols = [c for c in matrix.columns if c not in META_COLS]
    ingroup_ids = [s.short for s in samples if s.group == 'IN' and s.short in proteome_cols]
    outgroup_ids = [s.short for s in samples if s.group == 'OUT' and s.short in proteome_cols]

    results = []
    for row in controls:
        cid = (row.get('control_id') or '').strip()
        cls = (row.get('class') or '').strip().lower()
        expected = (row.get('expected_call') or '').strip().lower()
        rep, note = resolve_anchor(row, member_to_rep, busco_map, controls_dir,
                                   profiles_hmm, cpus)

        actual = 'unresolved'
        if rep is not None:
            presence = family_presence_vector(matrix, rep_to_members.get(rep, []),
                                              proteome_cols)
            if presence is None:
                rep, note = None, 'family has no matrix rows'
            else:
                actual = family_call(presence, ingroup_ids, outgroup_ids,
                                     ingroup_min_frac, other_max_frac)

        # outcome: hit/miss for positives, tn/fp for negatives, unresolved otherwise.
        if actual == 'unresolved':
            outcome = 'unresolved'
        elif cls == 'positive':
            outcome = 'hit' if actual == 'novel' else 'miss'
        elif cls == 'negative':
            outcome = 'fp' if actual == 'novel' else 'tn'
        else:
            outcome = 'n/a'

        results.append({
            'control_id': cid,
            'class': cls,
            'expected_call': expected,
            'anchor_type': (row.get('anchor_type') or '').strip(),
            'resolved_family': rep or '',
            'actual_call': actual,
            'outcome': outcome,
            'note': note,
        })
    return results


def summarize(results):
    pos = [r for r in results if r['class'] == 'positive']
    neg = [r for r in results if r['class'] == 'negative']
    pos_res = [r for r in pos if r['actual_call'] != 'unresolved']
    neg_res = [r for r in neg if r['actual_call'] != 'unresolved']
    hits = sum(1 for r in pos_res if r['outcome'] == 'hit')
    fps = sum(1 for r in neg_res if r['outcome'] == 'fp')
    recall = (hits / len(pos_res)) if pos_res else None
    fp_rate = (fps / len(neg_res)) if neg_res else None
    return {
        'n_controls': len(results),
        'n_positive': len(pos),
        'n_positive_resolved': len(pos_res),
        'positive_hits': hits,
        'recall': recall,
        'n_negative': len(neg),
        'n_negative_resolved': len(neg_res),
        'negative_fp': fps,
        'fp_rate': fp_rate,
        'n_unresolved': sum(1 for r in results if r['actual_call'] == 'unresolved'),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--controls', required=True, help='configs/controls/<clade>.controls.csv')
    ap.add_argument('--matrix', required=True, help="run's presence_matrix.tsv")
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv',
                    help='mmseqs *_cluster.tsv (rep_id<TAB>member_id)')
    ap.add_argument('--families', required=True,
                    help='families.tsv (profiled families index)')
    ap.add_argument('--config', required=True, help='analysis description CSV (IN/OUT groups)')
    ap.add_argument('--profiles', default=None,
                    help='family_profiles.hmm — required only for fasta anchors')
    ap.add_argument('--busco-map', default=None, dest='busco_map',
                    help='busco_id<TAB>protein_id TSV — required only for busco anchors (#6)')
    ap.add_argument('--ingroup-min-frac', type=float, default=0.75, dest='ingroup_min_frac')
    ap.add_argument('--other-max-frac', type=float, default=0.0, dest='other_max_frac')
    ap.add_argument('--cpus', type=int, default=1, help='hmmsearch --cpu for fasta anchors')
    ap.add_argument('--output', required=True, help='per-control results TSV')
    ap.add_argument('--summary', default=None,
                    help='optional summary TSV (defaults next to --output as *.summary.tsv)')
    args = ap.parse_args()

    controls, n_skipped = load_controls(args.controls)
    matrix = pd.read_csv(args.matrix, sep='\t')
    samples = parse_config(args.config)
    profiled_reps = load_profiled_reps(args.families)
    member_to_rep, rep_to_members = load_family_membership(args.cluster_tsv, profiled_reps)
    busco_map = load_busco_map(args.busco_map)
    controls_dir = Path(args.controls).resolve().parent

    results = score_controls(
        controls, matrix, member_to_rep, rep_to_members, samples,
        args.ingroup_min_frac, args.other_max_frac, busco_map, controls_dir,
        args.profiles, args.cpus,
    )

    fields = ['control_id', 'class', 'expected_call', 'anchor_type',
              'resolved_family', 'actual_call', 'outcome', 'note']
    with open(args.output, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(results)

    summary = summarize(results)
    summary['n_placeholder_skipped'] = n_skipped
    summary_path = args.summary or (str(Path(args.output).with_suffix('')) + '.summary.tsv')
    with open(summary_path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['metric', 'value'])
        for k, v in summary.items():
            w.writerow([k, '' if v is None else v])

    recall = summary['recall']
    fp_rate = summary['fp_rate']
    print(f"controls: {summary['n_controls']} scored, {n_skipped} placeholders skipped, "
          f"{summary['n_unresolved']} unresolved", file=sys.stderr)
    print(f"recall (positives): "
          f"{'n/a' if recall is None else f'{recall:.3f}'} "
          f"({summary['positive_hits']}/{summary['n_positive_resolved']})", file=sys.stderr)
    print(f"fp-rate (negatives): "
          f"{'n/a' if fp_rate is None else f'{fp_rate:.3f}'} "
          f"({summary['negative_fp']}/{summary['n_negative_resolved']})", file=sys.stderr)


if __name__ == '__main__':
    main()
