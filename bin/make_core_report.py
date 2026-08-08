#!/usr/bin/env python3
"""
Render a self-contained interactive HTML report of CORE genes — proteins
present in nearly every sampled proteome (ingroup and outgroup alike).

Needs no new search or annotation step: it re-reads the same annotated
presence matrix used by make_report.py and asks the opposite question of the
novelty report. The page has no external dependencies, so it can be copied
off the cluster and opened from a file:// URL.

Example:
  make_core_report.py \
      --matrix results/pezio4_asco/presence_matrix.function.tsv \
      --config configs/pezio4_asco.csv \
      --cluster_tsv results/pezio4_asco/clusters/clusters_cluster.tsv \
      --output results/pezio4_asco/core.html
"""
import argparse
import html
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import parse_config  # noqa: E402
from core_report_template import CORE_HTML_TEMPLATE  # noqa: E402
from gff3_genes import resolve_gff3_paths  # noqa: E402
from report_data import build_core_payload  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--matrix', required=True,
                    help='presence_matrix.function.tsv (annotated) or presence_matrix.tsv')
    ap.add_argument('--config', required=True,
                    help='Analysis config CSV (GROUP, Species, Strain, Short, TaxonGroup)')
    ap.add_argument('--cluster_tsv',
                    help='mmseqs easy-cluster *_cluster.tsv (rep -> member), used to group '
                         'candidates into gene families across ingroup species (optional)')
    ap.add_argument('--core_min_frac', type=float, default=0.95,
                    help='Minimum presence fraction across all proteomes (ingroup + outgroup) '
                         'for a gene to count as core (default: 0.95)')
    ap.add_argument('--project', default=None,
                    help='Project name shown in the report title (default: matrix parent dir)')
    ap.add_argument('--data_dir',
                    help='Directory containing the FASTA/GFF3 files referenced by --config '
                         '(optional; used only to resolve each species\' GFF3 column value '
                         'for the report\'s chrom/start columns).')
    ap.add_argument('--output', required=True, help='Output HTML file')
    args = ap.parse_args()

    if not (0.0 < args.core_min_frac <= 1.0):
        sys.exit('--core_min_frac must be in (0, 1]')

    samples = parse_config(args.config)
    project = args.project or Path(args.matrix).resolve().parent.name

    payload = build_core_payload(
        matrix_path=args.matrix,
        config_samples=samples,
        cluster_tsv=args.cluster_tsv,
        core_min_frac=args.core_min_frac,
        project=project,
        gff3_paths=resolve_gff3_paths(samples, args.data_dir),
    )

    # separators: drop the whitespace json.dumps adds after every delimiter.
    payload_json = json.dumps(payload, separators=(',', ':'))
    # The payload lives in a <script> block, so a literal "</script" inside any
    # annotation string would close it early.
    payload_json = payload_json.replace('</', '<\\/')

    doc = (CORE_HTML_TEMPLATE
           .replace('__PROJECT_TITLE__', html.escape(project))
           .replace('/*__PAYLOAD__*/', payload_json))

    out = Path(args.output)
    if out.parent != Path(''):
        out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(doc, encoding='utf-8')

    n_rows = len(payload['rows'])
    size_mb = out.stat().st_size / 1e6
    print(f'Wrote {out} ({size_mb:.2f} MB): {n_rows} core genes '
          f'(>= {args.core_min_frac:.0%} presence)', file=sys.stderr)


if __name__ == '__main__':
    main()
