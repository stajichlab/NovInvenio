"""
Assemble the interactive report payload from pipeline outputs.

The payload is a JSON-serialisable dict embedded into report.html by
bin/make_report.py.  Rows are stored column-wise as plain arrays (see ROW_FIELDS)
rather than dicts, because a per-row dict repeats every key ~20k times and roughly
triples the embedded JSON size.

No I/O happens at import time.
"""
import csv
import re
from pathlib import Path

# Order of the per-protein arrays in payload['rows'].  The browser reads these
# positionally via the same names exported in payload['fields'].
ROW_FIELDS = [
    'id',        # protein_id
    'src',       # index into payload['proteomes']
    'pres',      # bitstring over payload['proteomes'] order
    'tb',        # bitstring over payload['tblastn_genomes'] order
    'gene',      # gene_name
    'prod',      # product_description
    'fsrc',      # index into payload['fsources'], or -1
    'sprot',     # Best_Swissprot
    'pfam_n',    # Pfam_Names (comma-separated)
    'pfam_a',    # Pfam_Accessions (comma-separated)
    'pfam_e',    # Pfam_Evalues (comma-separated)
    'nov',       # 1 if this protein is a novelty candidate for its source proteome
    'seq',       # protein sequence ('' when not loaded)
]

# Trailing transcript/protein suffixes that separate a FungiDB gene ID from the
# per-transcript protein ID used in the proteome FASTAs.
#   Afu1g01620-T-p1      -> Afu1g01620
#   NCU00499-t26_1-p1    -> NCU00499
_GENE_ID_SUFFIX = re.compile(r'-[Tt][^-]*(-p\d+)?$|-p\d+$')


def gene_id_from_protein_id(protein_id: str) -> str:
    """Strip transcript/protein suffixes to recover a likely source gene ID."""
    return _GENE_ID_SUFFIX.sub('', protein_id)


def read_matrix(path: str | Path) -> tuple[list[str], list[dict]]:
    """Return (fieldnames, rows) from an annotated presence matrix TSV."""
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        return list(reader.fieldnames or []), list(reader)


def read_tblastn_summary(path: str | Path | None) -> tuple[list[str], dict[str, set]]:
    """Return (genome_ids, {protein_id: set_of_genome_ids_with_hits})."""
    if not path or not Path(path).exists():
        return [], {}
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        genomes = [c for c in (reader.fieldnames or []) if c != 'protein_id']
        hits = {}
        for row in reader:
            hit = {g for g in genomes if row.get(g, '0') == '1'}
            if hit:
                hits[row['protein_id']] = hit
    return genomes, hits


def read_novelties(paths: list[str | Path]) -> tuple[set[str], dict[str, str]]:
    """Return (novelty_protein_ids, {protein_id: sequence}) from novelties.<SHORT>.tsv.

    The sequence map is only populated when the files carry a protein_sequence
    column; older pipeline runs omit it.
    """
    ids: set[str] = set()
    seqs: dict[str, str] = {}
    for path in paths:
        with open(path, newline='') as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                pid = row.get('protein_id', '')
                if not pid:
                    continue
                ids.add(pid)
                seq = (row.get('protein_sequence') or '').strip()
                if seq:
                    seqs[pid] = seq
    return ids, seqs


def read_fasta_seqs(path: str | Path | None) -> dict[str, str]:
    """Minimal FASTA reader — {id: sequence}, ID taken up to the first whitespace.

    Deliberately not Bio.SeqIO: this only needs raw strings, and candidates.fa
    can hold tens of thousands of records.
    """
    if not path or not Path(path).exists():
        return {}
    seqs: dict[str, str] = {}
    current, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if current:
                    seqs[current] = ''.join(chunks)
                current, chunks = line[1:].split(None, 1)[0], []
            elif current:
                chunks.append(line.strip())
    if current:
        seqs[current] = ''.join(chunks)
    return seqs


def derive_novelties(rows, ingroup_ids, outgroup_ids, ingroup_min_frac) -> set[str]:
    """Recompute novelty status when no novelties.<SHORT>.tsv files are supplied.

    Mirrors bin/make_novelties.py as it is wired in workflows/summarize.nf, i.e.
    with the TBLASTN filter skipped: originates from an ingroup species, present
    in >= ingroup_min_frac of the ingroup, absent from every outgroup proteome.
    """
    novel = set()
    if not ingroup_ids:
        return novel
    for row in rows:
        if row.get('source_proteome', '') not in ingroup_ids:
            continue
        present = sum(1 for c in ingroup_ids if row.get(c, '0') == '1')
        if present / len(ingroup_ids) < ingroup_min_frac:
            continue
        if any(row.get(c, '0') == '1' for c in outgroup_ids):
            continue
        novel.add(row.get('protein_id', ''))
    return novel


def build_payload(
    matrix_path,
    config_samples,
    tblastn_path=None,
    novelty_paths=None,
    candidates_fa=None,
    ingroup_min_frac=0.75,
    project='NovInvenio',
    sequences='novelties',
) -> dict:
    """Build the embedded report payload.

    sequences: 'novelties' (default), 'all', or 'none' — which rows carry a
    protein_sequence.  Sequences dominate payload size, so 'all' is opt-in.
    """
    header, rows = read_matrix(matrix_path)

    # Only proteomes that actually have a column in the matrix are shown; ingroup
    # first so the heatmap's column blocks read left-to-right IN then OUT.
    ingroup = [s for s in config_samples if s.group == 'IN' and s.short in header]
    outgroup = [s for s in config_samples if s.group == 'OUT' and s.short in header]
    proteomes = ingroup + outgroup
    shorts = [s.short for s in proteomes]
    ingroup_ids = [s.short for s in ingroup]
    outgroup_ids = [s.short for s in outgroup]

    if not shorts:
        raise ValueError(
            f'No proteome columns from the config are present in {matrix_path}. '
            f'Matrix columns: {header}'
        )

    tb_genomes, tb_hits = read_tblastn_summary(tblastn_path)
    # Keep the TBLASTN block ordered like the outgroup block, and drop any genome
    # the config does not describe.
    tb_genomes = [g for g in outgroup_ids if g in tb_genomes]

    novelty_ids, novelty_seqs = read_novelties(novelty_paths or [])
    if not novelty_paths:
        novelty_ids = derive_novelties(rows, ingroup_ids, outgroup_ids, ingroup_min_frac)

    seq_lookup = dict(novelty_seqs)
    if sequences != 'none' and candidates_fa:
        # candidates.fa covers every candidate, not just novelties; it does not
        # override sequences already read from the novelties tables.
        for pid, seq in read_fasta_seqs(candidates_fa).items():
            seq_lookup.setdefault(pid, seq)

    fsources: list[str] = []
    fsource_idx: dict[str, int] = {}
    out_rows = []

    for row in rows:
        pid = row.get('protein_id', '')
        src = row.get('source_proteome', '')
        pres = ''.join('1' if row.get(s, '0') == '1' else '0' for s in shorts)
        hit_genomes = tb_hits.get(pid, set())
        tb = ''.join('1' if g in hit_genomes else '0' for g in tb_genomes)

        fsrc = row.get('function_source', '') or ''
        if fsrc:
            if fsrc not in fsource_idx:
                fsource_idx[fsrc] = len(fsources)
                fsources.append(fsrc)
            fsrc_i = fsource_idx[fsrc]
        else:
            fsrc_i = -1

        is_nov = 1 if pid in novelty_ids else 0
        if sequences == 'all' or (sequences == 'novelties' and is_nov):
            seq = seq_lookup.get(pid, '')
        else:
            seq = ''

        out_rows.append([
            pid,
            shorts.index(src) if src in shorts else -1,
            pres,
            tb,
            row.get('gene_name', '') or '',
            row.get('product_description', '') or '',
            fsrc_i,
            row.get('Best_Swissprot', '') or '',
            row.get('Pfam_Names', '') or '',
            row.get('Pfam_Accessions', '') or '',
            row.get('Pfam_Evalues', '') or '',
            is_nov,
            seq,
        ])

    return {
        'project': project,
        'ingroup_min_frac': ingroup_min_frac,
        'fields': ROW_FIELDS,
        'proteomes': [
            {
                'short': s.short,
                'species': s.species,
                'strain': s.strain,
                'group': s.group,
                'taxon': s.taxon_group,
            }
            for s in proteomes
        ],
        'tblastn_genomes': tb_genomes,
        'fsources': fsources,
        'rows': out_rows,
    }
