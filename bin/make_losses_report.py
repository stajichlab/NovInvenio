#!/usr/bin/env python3
"""
Render a self-contained interactive HTML report of candidate gene losses —
genes present in the outgroup but absent from the ingroup, prioritized as
possible lineage-specific loss pathways.

Reads the loss-search outputs (see workflows/loss_search.nf): a presence
matrix built with --query-group OUT, and (optionally) a TBLASTN summary of
the loss candidates against ingroup genomic DNA. The page has no external
dependencies, so it can be copied off the cluster and opened from a
file:// URL.

Example:
  make_losses_report.py \
      --matrix results/pezio4_asco/loss_presence_matrix.function.tsv \
      --config configs/pezio4_asco.csv \
      --tblastn_summary results/pezio4_asco/loss_tblastn_summary.tsv \
      --cluster_tsv results/pezio4_asco/clusters/loss_clusters_cluster.tsv \
      --output results/pezio4_asco/losses.html
"""
import argparse
import html
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import parse_config  # noqa: E402
from losses_report_template import LOSSES_HTML_TEMPLATE  # noqa: E402
from report_data import build_losses_payload  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--matrix', required=True,
                    help='loss_presence_matrix.function.tsv (annotated) or the plain matrix')
    ap.add_argument('--config', required=True,
                    help='Analysis config CSV (GROUP, Species, Strain, Short, TaxonGroup)')
    ap.add_argument('--tblastn_summary',
                    help='loss_tblastn_summary.tsv — TBLASTN of loss candidates against '
                         'ingroup genomic DNA (optional; the ingroup-TBLASTN column is '
                         'omitted without it)')
    ap.add_argument('--cluster_tsv',
                    help='loss mmseqs easy-cluster *_cluster.tsv (rep -> member), used to '
                         'group candidates into gene families across outgroup species (optional)')
    ap.add_argument('--outgroup_min_frac', type=float, default=0.75,
                    help='Minimum fraction of outgroup proteomes a loss candidate must be '
                         'present in (conservation threshold). Applied here to trim the '
                         'annotated matrix, which is the full presence table rather than the '
                         'filtered candidate list; pass the same value the loss search used '
                         '(default: 0.75)')
    ap.add_argument('--loss_ingroup_max_frac', type=float, default=0.0,
                    help='Maximum fraction of ingroup proteomes a loss candidate may still be '
                         'present in (0.0 = strictly absent). Must match build_presence_matrix.py '
                         "--other-max-frac from the loss search so reported rows match the "
                         'clustered candidate set (default: 0.0)')
    ap.add_argument('--project', default=None,
                    help='Project name shown in the report title (default: matrix parent dir)')
    ap.add_argument('--output', required=True, help='Output HTML file')
    args = ap.parse_args()

    if not (0.0 < args.outgroup_min_frac <= 1.0):
        sys.exit('--outgroup_min_frac must be in (0, 1]')
    if not (0.0 <= args.loss_ingroup_max_frac < 1.0):
        sys.exit('--loss_ingroup_max_frac must be in [0, 1)')

    samples = parse_config(args.config)
    project = args.project or Path(args.matrix).resolve().parent.name

    payload = build_losses_payload(
        matrix_path=args.matrix,
        config_samples=samples,
        tblastn_path=args.tblastn_summary,
        cluster_tsv=args.cluster_tsv,
        outgroup_min_frac=args.outgroup_min_frac,
        loss_ingroup_max_frac=args.loss_ingroup_max_frac,
        project=project,
    )

    payload_json = json.dumps(payload, separators=(',', ':'))
    payload_json = payload_json.replace('</', '<\\/')

    doc = (LOSSES_HTML_TEMPLATE
           .replace('__PROJECT_TITLE__', html.escape(project))
           .replace('/*__PAYLOAD__*/', payload_json))

    out = Path(args.output)
    if out.parent != Path(''):
        out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(doc, encoding='utf-8')

    n_rows = len(payload['rows'])
    n_flagged = sum(1 for r in payload['rows'] if r[payload['fields'].index('tb_hit')])
    size_mb = out.stat().st_size / 1e6
    print(f'Wrote {out} ({size_mb:.2f} MB): {n_rows} loss candidates '
          f'({n_flagged} flagged by ingroup TBLASTN)', file=sys.stderr)


if __name__ == '__main__':
    main()
