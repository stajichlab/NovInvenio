#!/usr/bin/env python3
"""
Render a self-contained interactive HTML report from the pipeline outputs.

Combines the annotated presence matrix, the TBLASTN summary and the per-species
novelties tables into a single HTML file with an embedded JSON payload.  The page
has no external dependencies, so it can be copied off the cluster and opened from
a file:// URL.

Example:
  make_report.py \
      --matrix results/pezio4_asco/presence_matrix.function.tsv \
      --config configs/pezio4_asco.csv \
      --tblastn_summary results/pezio4_asco/tblastn_summary.tsv \
      --novelties results/pezio4_asco/novelties.*.tsv \
      --candidates_fa results/pezio4_asco/candidates.fa \
      --cluster_tsv results/pezio4_asco/clusters/clusters_cluster.tsv \
      --evalues results/pezio4_asco/presence_matrix.evalues.tsv \
      --output results/pezio4_asco/novelties.html
"""
import argparse
import html
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import parse_config  # noqa: E402
from report_data import build_payload  # noqa: E402
from report_template import HTML_TEMPLATE  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--matrix', required=True,
                    help='presence_matrix.function.tsv (annotated) or presence_matrix.tsv')
    ap.add_argument('--config', required=True,
                    help='Analysis config CSV (GROUP, Species, Strain, Short, TaxonGroup)')
    ap.add_argument('--tblastn_summary',
                    help='tblastn_summary.tsv (optional; TBLASTN columns are omitted without it)')
    ap.add_argument('--novelties', nargs='*', default=[],
                    help='novelties.<SHORT>.tsv files (optional; novelty status is '
                         'recomputed from the matrix when omitted)')
    ap.add_argument('--candidates_fa',
                    help='candidates.fa, used to fill in protein sequences (optional)')
    ap.add_argument('--cluster_tsv',
                    help='mmseqs easy-cluster *_cluster.tsv (rep -> member), used to group '
                         'candidates into gene families across ingroup species (optional)')
    ap.add_argument('--evalues',
                    help='presence_matrix.evalues.tsv sidecar (optional; same shape as '
                         '--matrix, e-value per proteome cell instead of 0/1) -- shown in the '
                         'detail panel to help validate a presence call. Missing/empty file '
                         'means no e-value evidence available for this run.')
    ap.add_argument('--context_matrix',
                    help='context_presence.tsv (optional; issue #48) -- NEAR_INGROUP/'
                         'BROAD_OUTGROUP presence for the candidate list only, report-only, '
                         'never affects novelty calling. Missing/empty means no context '
                         'evidence available (e.g. --cluster_tool other than pairwise, or no '
                         'NEAR_INGROUP/BROAD_OUTGROUP rows in the config).')
    ap.add_argument('--context_evalues',
                    help='context_presence.evalues.tsv sidecar for --context_matrix')
    ap.add_argument('--project', default=None,
                    help='Project name shown in the report title (default: matrix parent dir)')
    ap.add_argument('--ingroup_min_frac', type=float, default=0.75,
                    help='Fraction of ingroup proteomes required for a novelty (default: 0.75)')
    ap.add_argument('--sequences', choices=['novelties', 'all', 'none'], default='novelties',
                    help='Which proteins carry an embedded sequence. Sequences dominate the '
                         'file size, so "all" is opt-in (default: novelties)')
    ap.add_argument('--method', default='pairwise',
                    help='Search pathway that produced --matrix, labelling this run in the '
                         'cross-method support column (default: pairwise)')
    ap.add_argument('--support_matrix',
                    help="The *other* pathway's presence_matrix.tsv. When given, each row's "
                         'support column records which method(s) call it novel — the '
                         'cross-method concordance check (ADR-0002 Phase 2)')
    ap.add_argument('--support_method', default=None,
                    help='Label for --support_matrix\'s pathway (default: mmseqs)')
    ap.add_argument('--output', required=True, help='Output HTML file')
    args = ap.parse_args()

    if not (0.0 < args.ingroup_min_frac <= 1.0):
        sys.exit('--ingroup_min_frac must be in (0, 1]')

    samples = parse_config(args.config)
    project = args.project or Path(args.matrix).resolve().parent.name

    # Nextflow expands a glob into separate tokens, but a shell that finds no
    # match passes the pattern through literally — drop anything unresolved.
    novelty_paths = [p for p in args.novelties if Path(p).exists()]
    if args.novelties and not novelty_paths:
        print('WARNING: no --novelties files matched; deriving novelty status from the matrix',
              file=sys.stderr)

    payload = build_payload(
        matrix_path=args.matrix,
        config_samples=samples,
        tblastn_path=args.tblastn_summary,
        novelty_paths=novelty_paths,
        candidates_fa=args.candidates_fa,
        cluster_tsv=args.cluster_tsv,
        evalues_path=args.evalues,
        context_matrix_path=args.context_matrix,
        context_evalues_path=args.context_evalues,
        ingroup_min_frac=args.ingroup_min_frac,
        project=project,
        sequences=args.sequences,
        method=args.method,
        support_matrix=args.support_matrix,
        support_method=args.support_method,
    )

    # separators: drop the whitespace json.dumps adds after every delimiter —
    # worth several hundred KB across ~20k rows.
    payload_json = json.dumps(payload, separators=(',', ':'))
    # The payload lives in a <script> block, so a literal "</script" inside any
    # annotation string would close it early.
    payload_json = payload_json.replace('</', '<\\/')

    doc = (HTML_TEMPLATE
           .replace('__PROJECT_TITLE__', html.escape(project))
           .replace('/*__PAYLOAD__*/', payload_json))

    out = Path(args.output)
    if out.parent != Path(''):
        out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(doc, encoding='utf-8')

    n_rows = len(payload['rows'])
    n_nov = sum(1 for r in payload['rows'] if r[payload['fields'].index('nov')])
    n_seq = sum(1 for r in payload['rows'] if r[payload['fields'].index('seq')])
    size_mb = out.stat().st_size / 1e6
    print(f'Wrote {out} ({size_mb:.1f} MB): {n_rows} proteins, {n_nov} novelty candidates, '
          f'{n_seq} sequences embedded', file=sys.stderr)


if __name__ == '__main__':
    main()
