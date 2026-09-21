"""bin/parse_hits.py: normalises raw phmmer/diamond/blast output to the common TSV
build_presence_matrix.py/context_presence.py/singleton_presence.py consume.

Issue #129: the wider diamond/blast --outfmt (length/pident/qcov/scov/qlen/slen) must
survive this normalisation step, or widening the raw tool output accomplishes nothing --
this is the step that was silently discarding them before this change, since it wrote a
hardcoded 6-column row. The new columns are carried through, not yet USED by any filter
(that is Option 4 in the paralog-rescue spec, deliberately still future work); a
narrow/phmmer hit leaves them BLANK, matching the existing blank-cell convention already
used by build_presence_matrix.py's --output-evalues/--output-targets sidecars.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / 'bin' / 'parse_hits.py'


def run(tmp_path, input_text, fmt):
    inp = tmp_path / 'raw.tsv'
    inp.write_text(input_text)
    out = tmp_path / 'parsed.tsv'
    subprocess.run([sys.executable, str(SCRIPT), '--input', str(inp), '--format', fmt,
                    '--output', str(out), '--evalue', '1', '--query-proteome', 'Q',
                    '--target-proteome', 'T'], check=True, capture_output=True, text=True)
    return pd.read_csv(out, sep='\t', dtype=str)


DIAMOND_WIDE = "q1\tt1\t1e-20\t150\t42\t35.5\t80\t90\t100\t45\n"
DIAMOND_NARROW = "q1\tt1\t1e-20\t150\n"


def test_narrow_diamond_input_gets_blank_metric_columns(tmp_path):
    df = run(tmp_path, DIAMOND_NARROW, 'diamond')
    assert list(df.columns) == ['query_id', 'target_id', 'evalue', 'bitscore',
                                 'query_proteome', 'target_proteome',
                                 'length', 'pident', 'qcov', 'scov', 'qlen', 'slen']
    row = df.iloc[0]
    for col in ('length', 'pident', 'qcov', 'scov', 'qlen', 'slen'):
        assert row[col] == '' or pd.isna(row[col]), f'{col} should be blank, got {row[col]!r}'


def test_wide_diamond_input_carries_the_metrics_through(tmp_path):
    df = run(tmp_path, DIAMOND_WIDE, 'diamond')
    row = df.iloc[0]
    assert row['length'] == '42'
    assert row['pident'] == '35.5'
    assert row['qcov'] == '80.0'
    assert row['scov'] == '90.0'
    assert row['qlen'] == '100'
    assert row['slen'] == '45'


def test_phmmer_input_gets_blank_metric_columns(tmp_path):
    tblout = (
        "#\n"
        "geneA - query1 - 1.2e-50 200.0 0.0 2.4e-50 199.0 0.0 1.0 1 0 0 1 1 1 1 desc\n"
    )
    df = run(tmp_path, tblout, 'phmmer')
    row = df.iloc[0]
    for col in ('length', 'pident', 'qcov', 'scov', 'qlen', 'slen'):
        assert row[col] == '' or pd.isna(row[col])


def test_original_six_columns_are_unchanged_by_position_or_value(tmp_path):
    # A caller reading this file by column NAME (every real consumer, via pandas) or by
    # the first four raw values must see identical results whether the metrics columns
    # are populated or blank.
    df = run(tmp_path, DIAMOND_WIDE, 'diamond')
    row = df.iloc[0]
    assert (row['query_id'], row['target_id'], row['evalue'], row['bitscore'],
            row['query_proteome'], row['target_proteome']) == \
           ('q1', 't1', '1e-20', '150.0', 'Q', 'T')
