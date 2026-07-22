#!/usr/bin/env python3
"""Render a publication-quality PDF summary of a NovInvenio run (issue #19).

A printable companion to the interactive HTML reports: it reuses the SAME payloads
lib/report_data.py builds (no new search or annotation), renders a fixed figure set with
matplotlib, and writes a multi-page PDF via PdfPages. Regenerable standalone from an
existing results directory, exactly like bin/make_report.py.

Example:
  make_pdf_report.py \
      --matrix results/pezizo5/presence_matrix.function.tsv \
      --config configs/pezizo5.csv \
      --tblastn_summary results/pezizo5/tblastn_summary.tsv \
      --cluster_tsv results/pezizo5/candidate_families_cluster.tsv \
      --loss_matrix results/pezizo5/loss_presence_matrix.function.tsv \
      --loss_cluster_tsv results/pezizo5/loss_candidate_families_cluster.tsv \
      --output view/pezizo5/summary.pdf
"""
import argparse
import sys
from pathlib import Path

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import parse_config  # noqa: E402
from report_data import build_payload, build_losses_payload  # noqa: E402
import figures  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--matrix', required=True,
                    help='presence_matrix.function.tsv (annotated) or presence_matrix.tsv')
    ap.add_argument('--config', required=True, help='analysis config CSV')
    ap.add_argument('--tblastn_summary', help='tblastn_summary.tsv (optional)')
    ap.add_argument('--cluster_tsv', help='candidate family cluster TSV (optional; gene families)')
    ap.add_argument('--ingroup_min_frac', type=float, default=0.75)
    ap.add_argument('--project', default=None,
                    help='project name shown on the title page (default: matrix parent dir)')
    # Loss-direction inputs (optional — adds the losses figure).
    ap.add_argument('--loss_matrix', help='loss_presence_matrix.function.tsv (optional)')
    ap.add_argument('--loss_tblastn_summary', help='loss_tblastn_summary.tsv (optional)')
    ap.add_argument('--loss_cluster_tsv', help='loss candidate family cluster TSV (optional)')
    ap.add_argument('--outgroup_min_frac', type=float, default=0.75)
    ap.add_argument('--loss_ingroup_max_frac', type=float, default=0.0)
    ap.add_argument('--output', required=True, help='output PDF file')
    args = ap.parse_args()

    samples = parse_config(args.config)
    project = args.project or Path(args.matrix).resolve().parent.name

    payload = build_payload(
        matrix_path=args.matrix,
        config_samples=samples,
        tblastn_path=args.tblastn_summary,
        cluster_tsv=args.cluster_tsv,
        ingroup_min_frac=args.ingroup_min_frac,
        project=project,
        sequences='none',  # the PDF never needs embedded sequences
    )

    losses_payload = None
    if args.loss_matrix and Path(args.loss_matrix).exists():
        losses_payload = build_losses_payload(
            matrix_path=args.loss_matrix,
            config_samples=samples,
            tblastn_path=args.loss_tblastn_summary,
            cluster_tsv=args.loss_cluster_tsv,
            outgroup_min_frac=args.outgroup_min_frac,
            loss_ingroup_max_frac=args.loss_ingroup_max_frac,
            project=project,
        )

    figs = figures.build_all(payload, losses_payload)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out) as pdf:
        for fig in figs:
            pdf.savefig(fig)
            plt.close(fig)
        meta = pdf.infodict()
        meta['Title'] = f'{project} — NovInvenio summary'
        meta['Creator'] = 'NovInvenio make_pdf_report.py'

    size_kb = out.stat().st_size / 1024
    print(f'Wrote {out} ({size_kb:.0f} KB): {len(figs)} figures', file=sys.stderr)


if __name__ == '__main__':
    main()
