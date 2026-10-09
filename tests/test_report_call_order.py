"""main.nf passes the REPORT workflow's inputs positionally. A swapped pair stages one input
under another's name (2026-10-08: query_lowcov received the TBLASTN coverage .gz and
MAKE_REPORT died with a UnicodeDecodeError). Check main.nf's call order against `take:`."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# main.nf argument -> the REPORT take name it must feed
EXPECT = {
    "novelty_query_lowcov": "query_lowcov",
    "novelty_other_evidence": "other_evidence",
    "novelty_tblastn_cov": "tblastn_coverage",
    "loss_other_evidence": "loss_other_evidence",
    "loss_tblastn_cov": "loss_tblastn_coverage",
}


def take_names():
    text = (ROOT / "workflows" / "report.nf").read_text()
    block = text.split("take:")[1].split("main:")[0]
    names = []
    for line in block.splitlines():
        m = re.match(r"\s{4}(\w+)\s*(//.*)?$", line)
        if m:
            names.append(m.group(1))
    return names


def call_args():
    text = (ROOT / "main.nf").read_text()
    body = text.split("    REPORT(\n")[1].split("\n    )")[0]
    return [a.strip().rstrip(",") for a in body.splitlines() if a.strip() and not a.strip().startswith("//")]


def test_report_call_order_matches_take_order():
    takes, args = take_names(), call_args()
    assert len(takes) == len(args)
    for arg, want in EXPECT.items():
        assert takes.index(want) == args.index(arg), (arg, want)
