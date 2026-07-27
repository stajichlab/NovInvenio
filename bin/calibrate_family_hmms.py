#!/usr/bin/env python3
"""Calibrate family HMM thresholds using the DISC_OUT panel as a negative control.

For each profiled family HMM, we search it against the DISC_OUT proteomes
(the discovery outgroup panel) and record the highest-scoring (lowest E-value)
hit.  The per-family threshold is set to that E-value, so that the family is
only called "present" in a downstream proteome if the hit is at least as
strong as the best false-positive observed in the negative control.

If a family HMM has no hit in any DISC_OUT proteome, the threshold falls back
to the global E-value parameter (--default-evalue).

Input:
  --domtblout  Per-proteome hmmsearch domtblout files (DISC_OUT proteomes)
  --families   families.tsv from extract_family_seqs.py (family_index, rep_id, n_members)

Output:
  TSV file: representative_id, threshold_evalue, calibration_source
  where calibration_source is 'negative_control' or 'global_default'.
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path


def parse_domtblout_best_evalue(path):
    """Return dict: query(HMM name) -> minimum full-sequence E-value.

    hmmsearch --domtblout columns (0-indexed): 0=target, 3=query, 6=full E-value.
    """
    best_e: dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            query = parts[3]
            try:
                evalue = float(parts[6])
            except ValueError:
                continue
            if query not in best_e or evalue < best_e[query]:
                best_e[query] = evalue
    return best_e


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--domtblout', nargs='+', required=True,
                    help='Per-proteome hmmsearch domtblout files (DISC_OUT scan)')
    ap.add_argument('--families', required=True,
                    help='families.tsv from extract_family_seqs.py')
    ap.add_argument('--default-evalue', type=float, default=1e-3,
                    dest='default_evalue',
                    help='Fallback E-value threshold for families with no outgroup hit')
    ap.add_argument('--output', required=True,
                    help='Output TSV: representative_id, threshold_evalue, calibration_source')
    args = ap.parse_args()

    # Collect the best (lowest) E-value per family HMM across all DISC_OUT proteomes
    family_best_e: dict[str, float] = {}
    for dom_path in args.domtblout:
        for query, evalue in parse_domtblout_best_evalue(dom_path).items():
            if query not in family_best_e or evalue < family_best_e[query]:
                family_best_e[query] = evalue

    # Read families.tsv to get all profiled family representative IDs
    reps = []
    with open(args.families) as fh:
        fh.readline()  # skip header
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                reps.append(parts[1])

    with open(args.output, 'w') as out:
        out.write('representative_id\tthreshold_evalue\tcalibration_source\n')
        for rep in reps:
            best_e = family_best_e.get(rep)
            if best_e is not None:
                out.write(f'{rep}\t{best_e:g}\tnegative_control\n')
            else:
                out.write(f'{rep}\t{args.default_evalue:g}\tglobal_default\n')

    n_calibrated = sum(1 for r in reps if family_best_e.get(r) is not None)
    print(f"Calibrated {n_calibrated}/{len(reps)} families from negative control, "
          f"{len(reps) - n_calibrated} using global default", file=sys.stderr)


if __name__ == '__main__':
    main()
