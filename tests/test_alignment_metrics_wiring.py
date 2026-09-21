"""Structural checks that the wider --outfmt (issue #129) reached every DIAMOND_SEARCH/
DIAMOND_SELF/BLAST_SEARCH/BLAST_SELF invocation -- and specifically NOT phmmer (#150,
a separate follow-up: --tblout carries no alignment geometry, and --domtblout is a
structurally different, per-domain output shape, not just a wider row of the same
format) and NOT diamond makedb (does not accept search-time options).
"""
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DIAMOND = (REPO / 'modules' / 'diamond.nf').read_text()
SELF_SEARCH = (REPO / 'modules' / 'self_search.nf').read_text()
BLAST = (REPO / 'modules' / 'blast.nf').read_text()
PHMMER = (REPO / 'modules' / 'phmmer.nf').read_text()

DIAMOND_WIDE_TAIL = 'length pident qcovhsp scovhsp qlen slen'
BLAST_WIDE_TAIL = 'length pident qcovhsp qlen slen'   # no scovhsp equivalent in blastp


def test_diamond_search_emits_the_wide_diamond_columns():
    block = DIAMOND.split('process DIAMOND_SEARCH', 1)[1]
    assert f'qseqid sseqid evalue bitscore {DIAMOND_WIDE_TAIL}' in block


def test_diamond_makedb_is_untouched():
    makedb = DIAMOND.split('process DIAMOND_SEARCH', 1)[0]
    assert 'qcovhsp' not in makedb


def test_diamond_self_emits_the_wide_diamond_columns():
    block = SELF_SEARCH.split('process DIAMOND_SELF', 1)[1].split('process BLAST_SELF')[0]
    assert f'qseqid sseqid evalue bitscore {DIAMOND_WIDE_TAIL}' in block


def test_diamond_self_stays_at_very_sensitive():
    # The wider outfmt must not disturb the deliberate --very-sensitive choice
    # (HEX-1/eIF-5A paralog detection).
    block = SELF_SEARCH.split('process DIAMOND_SELF', 1)[1].split('process BLAST_SELF')[0]
    assert '--very-sensitive' in block


def test_blast_search_emits_the_wide_blast_columns():
    block = BLAST.split('process BLAST_SEARCH', 1)[1] if 'process BLAST_SEARCH' in BLAST else BLAST
    assert f'qseqid sseqid evalue bitscore {BLAST_WIDE_TAIL}' in block
    # blastp has no subject-coverage keyword -- must not appear anywhere in this block
    assert 'scovhsp' not in block


def test_blast_self_emits_the_wide_blast_columns():
    block = SELF_SEARCH.split('process BLAST_SELF', 1)[1]
    assert f'qseqid sseqid evalue bitscore {BLAST_WIDE_TAIL}' in block
    assert 'scovhsp' not in block


def test_phmmer_search_is_not_touched_by_this_issue():
    assert 'qcovhsp' not in PHMMER
    assert '--tblout' in PHMMER, 'still per-sequence tblout -- #150 tracks moving to domtblout'


def test_phmmer_self_is_not_touched_by_this_issue():
    block = SELF_SEARCH.split('process PHMMER_SELF', 1)[1].split('process DIAMOND_SELF')[0]
    assert 'qcovhsp' not in block
    assert '--tblout' in block
