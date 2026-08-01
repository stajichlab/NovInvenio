import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))

from split_fasta_chunks import split_fasta_chunks  # noqa: E402


def _fasta(n):
    return ''.join(f'>seq{i}\nAAAA\n' for i in range(n))


def test_splits_into_fixed_size_chunks(tmp_path):
    inp = tmp_path / 'in.fa'
    inp.write_text(_fasta(7))
    outdir = tmp_path / 'out'

    split_fasta_chunks(str(inp), chunk_size=3, outdir=str(outdir), prefix='c_')

    chunks = sorted(outdir.glob('c_*.fa'))
    assert [c.name for c in chunks] == ['c_000000.fa', 'c_000001.fa', 'c_000002.fa']
    counts = [chunk.read_text().count('>') for chunk in chunks]
    assert counts == [3, 3, 1]
    # Every input sequence appears exactly once across the chunks.
    all_ids = set()
    for chunk in chunks:
        all_ids.update(line[1:].strip() for line in chunk.read_text().splitlines()
                       if line.startswith('>'))
    assert all_ids == {f'seq{i}' for i in range(7)}


def test_exact_multiple_of_chunk_size(tmp_path):
    inp = tmp_path / 'in.fa'
    inp.write_text(_fasta(6))
    outdir = tmp_path / 'out'

    split_fasta_chunks(str(inp), chunk_size=3, outdir=str(outdir), prefix='c_')

    chunks = sorted(outdir.glob('c_*.fa'))
    assert len(chunks) == 2
    assert all(chunk.read_text().count('>') == 3 for chunk in chunks)


def test_empty_input_emits_one_empty_chunk(tmp_path):
    inp = tmp_path / 'in.fa'
    inp.write_text('')
    outdir = tmp_path / 'out'

    split_fasta_chunks(str(inp), chunk_size=3, outdir=str(outdir), prefix='c_')

    chunks = sorted(outdir.glob('c_*.fa'))
    assert [c.name for c in chunks] == ['c_000000.fa']
    assert chunks[0].read_text() == ''
