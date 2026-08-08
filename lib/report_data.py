"""
Assemble the interactive report payload from pipeline outputs.

The payload is a JSON-serialisable dict embedded into novelties.html by
bin/make_report.py.  Rows are stored column-wise as plain arrays (see ROW_FIELDS)
rather than dicts, because a per-row dict repeats every key ~20k times and roughly
triples the embedded JSON size.

No I/O happens at import time.
"""
import csv
from pathlib import Path

from clusters import FamilyIndex
from config_parser import INGROUP_ROLES, OUTGROUP_ROLES  # noqa: F401 (re-exported)
from gff3_genes import (  # noqa: F401 (gene_id_from_protein_id re-exported for
    gene_id_from_protein_id,                              # backward compatibility --
    load_gff3_index,                                       # see gff3_genes.py's
    lookup_gene_position,                                  # docstring for why it moved)
)

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
    'fam',       # index into payload['families'], or -1 if not part of a multi-member cluster
    'seq',       # protein sequence ('' when not loaded)
    'support',   # which search method(s) called this protein novel — e.g. 'pairwise',
                 # 'mmseqs', or 'pairwise+mmseqs' (cross-method concordance). '' if no
                 # method calls it. Single-method runs carry just the one method name.
    'category',  # novelty_category from novelty_screen.py (issue #27): 'target_specific',
                 # 'clade_specific', 'false_novelty', or '' — either not a phase-1 candidate,
                 # or this run's matrix predates/doesn't come from the novelty_discovery /
                 # novelty_screen pathway (--cluster_tool pairwise/mmseqs never populate it).
    'ev',        # search-hit e-values (issue #44), aligned position-for-position with 'pres'
                 # — comma-separated, one entry per payload['proteomes'] column, '' where
                 # there's no e-value evidence (absent, self-sourced, or the run's pathway
                 # doesn't track e-values). Report-only: never affects presence/novelty calls.
                 # 'pres'/'ev' both extend to cover context columns too (issue #48) —
                 # payload['proteomes'] entries tagged {'context': true} are appended after
                 # the scored ingroup/outgroup columns, so 'pres'/'ev' stay one bit/entry per
                 # payload['proteomes'] column throughout.
    'chrom',     # GFF3-derived chromosome/scaffold/contig name for this row's source
                 # protein, '' when the source species has no --config GFF3 column value,
                 # the file can't be resolved/parsed, or the protein ID has no match in it.
    'start',     # GFF3-derived 1-based gene/mRNA start coordinate (int), or null/None when
                 # unknown (same conditions as 'chrom'). Per-protein-record as currently
                 # modeled -- see CLAUDE.md's GFF3 chrom/start note on splice isoforms.
]

# gene_id_from_protein_id() / _GENE_ID_SUFFIX moved to lib/gff3_genes.py --
# lookup_gene_position() (gff3_genes.py) needs it, and this module needs
# lookup_gene_position(), so it lives on the side that doesn't import back.
# Imported above and re-exported for any external caller still importing it
# from here (e.g. tests/test_report_data.py).


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


def read_evalues(path: str | Path | None) -> dict[str, dict[str, str]]:
    """Return {protein_id: {proteome_short: evalue_str}} from an evalues sidecar TSV.

    Same shape as the presence matrix (protein_id, source_proteome, <proteome columns>)
    but cells hold the qualifying hit's e-value (or '') instead of 0/1 — see
    bin/build_presence_matrix.py / bin/novelty_presence_matrix.py's --output-evalues.
    Missing, unreadable, or empty/header-only files (e.g. EMPTY_EVALUES_STUB for pathways
    that don't track e-values yet) return an empty dict — callers treat that as "no
    evidence available" rather than an error.
    """
    if not path or not Path(path).exists() or not Path(path).stat().st_size:
        return {}
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        cols = [c for c in (reader.fieldnames or []) if c not in ('protein_id', 'source_proteome')]
        out: dict[str, dict[str, str]] = {}
        for row in reader:
            pid = row.get('protein_id', '')
            if not pid:
                continue
            out[pid] = {c: row.get(c, '') or '' for c in cols}
    return out


def read_context(
    matrix_path: str | Path | None, evalues_path: str | Path | None
) -> tuple[list[str], dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    """Return (context_shorts, presence_by_pid, evalue_by_pid) from a context_presence.tsv
    pair (issue #48) -- NEAR_INGROUP/BROAD_OUTGROUP presence for the candidate list only,
    produced by bin/context_presence.py. Report-only: these columns never determine
    novelty. Missing/empty/header-only input (e.g. --cluster_tool other than pairwise, or
    a config with no NEAR_INGROUP/BROAD_OUTGROUP rows) returns ([], {}, {}) -- callers
    treat that as "no context evidence for this run".
    """
    if not matrix_path or not Path(matrix_path).exists() or not Path(matrix_path).stat().st_size:
        return [], {}, {}
    with open(matrix_path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        shorts = [c for c in (reader.fieldnames or []) if c not in ('protein_id', 'source_proteome')]
        presence: dict[str, dict[str, str]] = {}
        for row in reader:
            pid = row.get('protein_id', '')
            if not pid:
                continue
            presence[pid] = {s: row.get(s, '0') or '0' for s in shorts}
    if not shorts:
        return [], {}, {}
    return shorts, presence, read_evalues(evalues_path)


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


def read_support_novelties(
    support_matrix_path, config_samples, ingroup_min_frac
) -> set[str]:
    """Novelty protein_ids called by a *second* pathway's presence matrix.

    Cross-method concordance (ADR-0002 Phase 2): the other pathway's run
    (pairwise vs mmseqs) emits the same matrix contract, so its novelty set is
    recomputed with derive_novelties() — the same TBLASTN-skipped predicate the
    primary side uses when no novelties.<SHORT>.tsv is supplied — and compared
    per protein_id.  Returns an empty set when the path is missing.
    """
    if not support_matrix_path or not Path(support_matrix_path).exists():
        return set()
    header, rows = read_matrix(support_matrix_path)
    ingroup_ids = [s.short for s in config_samples if s.group in INGROUP_ROLES and s.short in header]
    outgroup_ids = [s.short for s in config_samples if s.group in OUTGROUP_ROLES and s.short in header]
    return derive_novelties(rows, ingroup_ids, outgroup_ids, ingroup_min_frac)


def _chrom_start(
    protein_id: str,
    source_short: str,
    gff3_paths: dict[str, str],
    gff3_cache: dict[str, dict],
) -> tuple[str, int | None]:
    """('' , None) unless source_short has a resolved GFF3 that yields a hit.

    gff3_cache is shared across every row in one payload build so each
    species' GFF3 file is parsed at most once (see gff3_genes.load_gff3_index).
    """
    path = gff3_paths.get(source_short)
    if not path:
        return '', None
    index = load_gff3_index(path, gff3_cache)
    hit = lookup_gene_position(protein_id, index)
    if hit is None:
        return '', None
    chrom, start = hit
    return chrom, start


def build_payload(
    matrix_path,
    config_samples,
    tblastn_path=None,
    novelty_paths=None,
    candidates_fa=None,
    cluster_tsv=None,
    evalues_path=None,
    context_matrix_path=None,
    context_evalues_path=None,
    ingroup_min_frac=0.75,
    project='NovInvenio',
    sequences='novelties',
    method='pairwise',
    support_matrix=None,
    support_method=None,
    gff3_paths=None,
) -> dict:
    """Build the embedded report payload.

    sequences: 'novelties' (default), 'all', or 'none' — which rows carry a
    protein_sequence.  Sequences dominate payload size, so 'all' is opt-in.

    cluster_tsv: mmseqs easy-cluster *_cluster.tsv (rep -> member), used to
    group candidates from different ingroup species into a gene family — a
    lineage-specific gene recovered independently in several ingroup taxa
    otherwise shows up as unrelated rows with no way to tell they're the same
    locus.  Clusters with a single member carry no family (fam == -1); there
    is nothing to collapse.

    method / support_matrix / support_method: cross-method concordance
    (ADR-0002 Phase 2).  `method` labels the pathway that produced matrix_path
    ('pairwise' or 'mmseqs').  When `support_matrix` (the *other* pathway's
    presence matrix) is given, each row's `support` field records which
    method(s) call that protein novel — 'pairwise', 'mmseqs' or
    'pairwise+mmseqs'.  Multi-method rows are high-confidence; disagreements
    flag threshold/boundary sensitivity.  Single-method runs leave `support` as
    just `method` for novelty rows.

    evalues_path (issue #44): optional presence_matrix.evalues.tsv sidecar —
    report-only search-hit e-value evidence, surfaced in the detail panel to
    help judge whether a presence call is a strong or marginal hit. Missing or
    empty means no evidence available; never affects presence/novelty calls.

    context_matrix_path / context_evalues_path (issue #48): optional
    context_presence.tsv pair — NEAR_INGROUP/BROAD_OUTGROUP presence for the
    candidate list only (bin/context_presence.py, --cluster_tool pairwise
    only). Appended to payload['proteomes'] as extra columns tagged
    {'context': True}, extending 'pres'/'ev' to cover them — but they are
    never counted toward ingroup/outgroup novelty stats.

    gff3_paths: optional {short: resolved GFF3 file path} (see
    lib/gff3_genes.resolve_gff3_paths) — supplies each row's 'chrom'/'start'
    from the source proteome's GFF3. A short with no entry (no --config GFF3
    value, or it didn't resolve under --data_dir) just leaves those rows'
    chrom/start blank; never an error.
    """
    header, rows = read_matrix(matrix_path)
    evalue_lookup = read_evalues(evalues_path)
    gff3_paths = gff3_paths or {}
    gff3_cache: dict[str, dict] = {}

    fam_index = FamilyIndex(cluster_tsv)

    # Only proteomes that actually have a column in the matrix are shown; ingroup
    # first so the heatmap's column blocks read left-to-right IN then OUT.
    ingroup = [s for s in config_samples if s.group in INGROUP_ROLES and s.short in header]
    outgroup = [s for s in config_samples if s.group in OUTGROUP_ROLES and s.short in header]
    proteomes = ingroup + outgroup
    shorts = [s.short for s in proteomes]
    ingroup_ids = [s.short for s in ingroup]
    outgroup_ids = [s.short for s in outgroup]

    # Context columns (issue #48): a NEAR_INGROUP/BROAD_OUTGROUP short that already has a
    # *scored* column in the main matrix (e.g. another producer's fix broadened matrix
    # columns via OUTGROUP_ROLES, even though it was never actually searched) must not be
    # appended a second time -- that would double-count it in payload['proteomes'] and
    # silently fold zero-only "context" evidence into the scored outgroup stats above.
    context_shorts, context_presence, context_evalue_lookup = read_context(
        context_matrix_path, context_evalues_path
    )
    context_samples = {s.short: s for s in config_samples if s.short in context_shorts}
    context_shorts = [s for s in context_shorts if s in context_samples and s not in shorts]

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

    # Cross-method concordance: which method(s) call each protein novel.
    support_method = support_method or 'mmseqs'
    support_novelty_ids = read_support_novelties(support_matrix, config_samples,
                                                 ingroup_min_frac)
    methods = [method] + ([support_method] if support_matrix else [])

    seq_lookup = dict(novelty_seqs)
    if sequences != 'none' and candidates_fa:
        # candidates.fa covers every candidate, not just novelties; it does not
        # override sequences already read from the novelties tables.
        for pid, seq in read_fasta_seqs(candidates_fa).items():
            seq_lookup.setdefault(pid, seq)

    fsources: list[str] = []
    fsource_idx: dict[str, int] = {}
    out_rows = []
    categories: set[str] = set()

    for row in rows:
        pid = row.get('protein_id', '')
        src = row.get('source_proteome', '')
        row_context_pres = context_presence.get(pid, {})
        row_context_ev = context_evalue_lookup.get(pid, {})
        pres = ''.join(
            ['1' if row.get(s, '0') == '1' else '0' for s in shorts] +
            ['1' if row_context_pres.get(s, '0') == '1' else '0' for s in context_shorts]
        )
        row_evalues = evalue_lookup.get(pid, {})
        ev = ','.join(
            [row_evalues.get(s, '') for s in shorts] +
            [row_context_ev.get(s, '') for s in context_shorts]
        )
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

        fam_i = fam_index.index_of(pid, src)

        # support = the methods that call this protein novel, primary method first.
        row_methods = ([method] if is_nov else [])
        if support_matrix and pid in support_novelty_ids:
            row_methods.append(support_method)
        support = '+'.join(row_methods)

        # novelty_category (issue #27/#28): only present when the matrix came from the
        # novelty_discovery/novelty_screen pathway; '' for pairwise/mmseqs runs.
        category = row.get('novelty_category', '') or ''
        if category:
            categories.add(category)

        chrom, start = _chrom_start(pid, src, gff3_paths, gff3_cache)

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
            fam_i,
            seq,
            support,
            category,
            ev,
            chrom,
            start,
        ])

    return {
        'project': project,
        'ingroup_min_frac': ingroup_min_frac,
        'fields': ROW_FIELDS,
        'methods': methods,
        'support_method': support_method if support_matrix else None,
        'proteomes': [
            {
                'short': s.short,
                'species': s.species,
                'strain': s.strain,
                'group': s.group,
                'taxon': s.taxon_group,
            }
            for s in proteomes
        ] + [
            {
                'short': context_samples[short].short,
                'species': context_samples[short].species,
                'strain': context_samples[short].strain,
                'group': context_samples[short].group,
                'taxon': context_samples[short].taxon_group,
                'context': True,
            }
            for short in context_shorts
        ],
        'tblastn_genomes': tb_genomes,
        'fsources': fsources,
        'novelty_categories': sorted(categories),
        'has_evalues': bool(evalue_lookup),
        'has_context': bool(context_shorts),
        'families': fam_index.payload(),
        'rows': out_rows,
    }


# Order of the per-protein arrays in a CORE payload's payload['rows'].
CORE_ROW_FIELDS = [
    'id',        # protein_id
    'src',       # index into payload['proteomes']
    'frac',      # presence fraction across every proteome column (ingroup + outgroup)
    'gene',      # gene_name
    'prod',      # product_description
    'fsrc',      # index into payload['fsources'], or -1
    'sprot',     # Best_Swissprot
    'pfam_n',    # Pfam_Names (comma-separated)
    'pfam_a',    # Pfam_Accessions (comma-separated)
    'fam',       # index into payload['families'], or -1 if not part of a multi-member cluster
    'chrom',     # GFF3-derived chromosome/scaffold/contig name (see ROW_FIELDS' 'chrom')
    'start',     # GFF3-derived 1-based start coordinate, int or null (see ROW_FIELDS' 'start')
]


def build_core_payload(
    matrix_path,
    config_samples,
    cluster_tsv=None,
    core_min_frac=0.95,
    project='NovInvenio',
    gff3_paths=None,
) -> dict:
    """Build the embedded payload for the CORE (near-universal genes) report.

    Unlike build_payload(), this needs no new search or annotation step — it
    re-reads the same annotated presence matrix and asks the opposite
    question: instead of "present in the ingroup, absent from the outgroup",
    "present in essentially every proteome, ingroup and outgroup alike".

    A row only exists in the matrix when it originated as an ingroup query
    (see build_presence_matrix.py), so this necessarily reports the *ingroup
    view* of conservation: a gene both ingroup and outgroup share, discovered
    from the ingroup side.  A truly core gene surfaces once per ingroup
    species that independently recovered it as a qualifying hit, so rows are
    grouped into families via the same mmseqs clustering used for novelty
    candidates (cluster_tsv), for a fuller picture of the redundancy.

    gff3_paths: optional {short: resolved GFF3 file path} — see build_payload()'s
    docstring; supplies each row's 'chrom'/'start' from its (always ingroup)
    source proteome's GFF3.
    """
    header, rows = read_matrix(matrix_path)

    fam_index = FamilyIndex(cluster_tsv)
    gff3_paths = gff3_paths or {}
    gff3_cache: dict[str, dict] = {}

    ingroup = [s for s in config_samples if s.group in INGROUP_ROLES and s.short in header]
    outgroup = [s for s in config_samples if s.group in OUTGROUP_ROLES and s.short in header]
    proteomes = ingroup + outgroup
    shorts = [s.short for s in proteomes]
    ingroup_ids = {s.short for s in ingroup}

    if not shorts:
        raise ValueError(
            f'No proteome columns from the config are present in {matrix_path}. '
            f'Matrix columns: {header}'
        )

    fsources: list[str] = []
    fsource_idx: dict[str, int] = {}
    out_rows = []

    for row in rows:
        src = row.get('source_proteome', '')
        if src not in ingroup_ids:
            continue

        present = sum(1 for s in shorts if row.get(s, '0') == '1')
        frac = present / len(shorts)
        if frac < core_min_frac:
            continue

        pid = row.get('protein_id', '')

        fsrc = row.get('function_source', '') or ''
        if fsrc:
            if fsrc not in fsource_idx:
                fsource_idx[fsrc] = len(fsources)
                fsources.append(fsrc)
            fsrc_i = fsource_idx[fsrc]
        else:
            fsrc_i = -1

        fam_i = fam_index.index_of(pid, src)
        chrom, start = _chrom_start(pid, src, gff3_paths, gff3_cache)

        out_rows.append([
            pid,
            shorts.index(src),
            round(frac, 4),
            row.get('gene_name', '') or '',
            row.get('product_description', '') or '',
            fsrc_i,
            row.get('Best_Swissprot', '') or '',
            row.get('Pfam_Names', '') or '',
            row.get('Pfam_Accessions', '') or '',
            fam_i,
            chrom,
            start,
        ])

    return {
        'project': project,
        'core_min_frac': core_min_frac,
        'fields': CORE_ROW_FIELDS,
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
        'fsources': fsources,
        'families': fam_index.payload(),
        'rows': out_rows,
    }


# Order of the per-protein arrays in a LOSSES payload's payload['rows'].
LOSSES_ROW_FIELDS = [
    'id',          # protein_id (an outgroup protein)
    'src',         # index into payload['proteomes']
    'frac',        # presence fraction across outgroup proteome columns only
    'out_breadth', # # distinct outgroup species carrying any member of this gene family
                   #   (or, for a singleton, this protein's own outgroup presence count)
    'in_retained', # # distinct ingroup species that still retain any member of this family
                   #   (or, for a singleton, this protein's own ingroup presence count) — 0 is a
                   #   clean loss; >0 means "nearly missing" (allowed by loss_ingroup_max_frac)
    'gene',        # gene_name (from the outgroup side's own annotation)
    'prod',        # product_description
    'fsrc',        # index into payload['fsources'], or -1
    'sprot',       # Best_Swissprot
    'pfam_n',      # Pfam_Names (comma-separated)
    'pfam_a',      # Pfam_Accessions (comma-separated)
    'tb_hit',      # 1 if TBLASTN found this outgroup protein in an ingroup genome
    'tb_genomes',  # comma-separated ingroup genome IDs with a TBLASTN hit
    'fam',         # index into payload['families'], or -1 if not part of a multi-member cluster
    'chrom',       # GFF3-derived chromosome/scaffold/contig name, resolved against the
                   # *outgroup* protein's GFF3 (that's where the gene is) — see ROW_FIELDS
    'start',       # GFF3-derived 1-based start coordinate, int or null — see ROW_FIELDS
]


def build_losses_payload(
    matrix_path,
    config_samples,
    tblastn_path=None,
    cluster_tsv=None,
    outgroup_min_frac=0.75,
    loss_ingroup_max_frac=0.0,
    project='NovInvenio',
    gff3_paths=None,
) -> dict:
    """Build the embedded payload for the LOSSES (candidate gene-loss) report.

    matrix_path is loss_presence_matrix.function.tsv — built by the same
    build_presence_matrix.py machinery as the novelty matrix, but with
    --query-group OUT (see workflows/loss_search.nf): rows are sourced from
    outgroup proteins.

    The *matrix* is the full presence/absence table (every outgroup query
    with any qualifying hit) — LOSS_ANNOTATE runs on it, not on the filtered
    candidate list. So, exactly like build_core_payload() and
    derive_novelties(), this recomputes the candidate predicate here rather
    than trusting the matrix's row set:

      - present in >= outgroup_min_frac of the outgroup  (conserved), AND
      - present in <= loss_ingroup_max_frac of the ingroup  (lost / nearly lost).

    loss_ingroup_max_frac mirrors build_presence_matrix.py's --other-max-frac
    (default 0.0 = strictly absent from the ingroup); pass the same value the
    loss search used so the reported rows match loss_candidates.txt and hence
    the mmseqs families built from it.

    Beyond the per-protein outgroup fraction (frac), each row also carries two
    *family-level* aggregates — the biological unit of a loss is the gene
    family, not one outgroup protein:
      - out_breadth: how many distinct outgroup species carry any member of the
        family (conservation breadth of the whole family).
      - in_retained: how many ingroup species still retain any member (0 = clean
        loss across the ingroup; >0 only when loss_ingroup_max_frac was raised).
    Singletons (no multi-member cluster) fall back to their own row's counts.

    tblastn_path, if given, is loss_tblastn_summary.tsv — TBLASTN of the
    loss candidates against ingroup *genomic* DNA (not just the annotated
    proteins), the same false-positive check the novelty report runs in the
    opposite direction: an absence-of-hit in the protein search does not
    prove absence in the genome, it may just be a missed gene model. A hit
    here does not remove the row (this is reporting-only, mirroring
    make_novelties.py's --skip_tblastn_filter default) — it downgrades the
    row's priority instead, surfaced via tb_hit/tb_genomes.

    gff3_paths: optional {short: resolved GFF3 file path} — see build_payload()'s
    docstring. Rows here are sourced from outgroup proteins, so 'chrom'/'start'
    are resolved against the row's own (outgroup) source proteome's GFF3 — the
    same field used everywhere else in this module, since that's already where
    the loss-candidate gene actually lives.
    """
    header, rows = read_matrix(matrix_path)

    fam_index = FamilyIndex(cluster_tsv)
    gff3_paths = gff3_paths or {}
    gff3_cache: dict[str, dict] = {}

    ingroup = [s for s in config_samples if s.group in INGROUP_ROLES and s.short in header]
    outgroup = [s for s in config_samples if s.group in OUTGROUP_ROLES and s.short in header]
    proteomes = ingroup + outgroup
    shorts = [s.short for s in proteomes]
    ingroup_ids = [s.short for s in ingroup]
    outgroup_ids = [s.short for s in outgroup]
    outgroup_index = {s.short for s in outgroup}

    if not shorts:
        raise ValueError(
            f'No proteome columns from the config are present in {matrix_path}. '
            f'Matrix columns: {header}'
        )

    tb_genomes, tb_hits = read_tblastn_summary(tblastn_path)

    # Pass 1: keep only rows that clear the loss predicate, and accumulate each
    # gene family's outgroup breadth / ingroup retention across its members.
    kept = []
    fam_out_species: dict[int, set] = {}
    fam_in_species: dict[int, set] = {}
    for row in rows:
        src = row.get('source_proteome', '')
        if src not in outgroup_index:
            continue

        out_present = [s for s in outgroup_ids if row.get(s, '0') == '1']
        in_present = [s for s in ingroup_ids if row.get(s, '0') == '1']
        out_frac = len(out_present) / len(outgroup_ids) if outgroup_ids else 0.0
        in_frac = len(in_present) / len(ingroup_ids) if ingroup_ids else 0.0
        if out_frac < outgroup_min_frac or in_frac > loss_ingroup_max_frac:
            continue

        pid = row.get('protein_id', '')
        fam_i = fam_index.index_of(pid, src)
        if fam_i >= 0:
            fam_out_species.setdefault(fam_i, set()).update(out_present)
            fam_in_species.setdefault(fam_i, set()).update(in_present)
        kept.append((row, pid, src, out_frac, out_present, in_present, fam_i))

    # Pass 2: emit rows, resolving the family-level aggregates.
    fsources: list[str] = []
    fsource_idx: dict[str, int] = {}
    out_rows = []
    for row, pid, src, out_frac, out_present, in_present, fam_i in kept:
        fsrc = row.get('function_source', '') or ''
        if fsrc:
            if fsrc not in fsource_idx:
                fsource_idx[fsrc] = len(fsources)
                fsources.append(fsrc)
            fsrc_i = fsource_idx[fsrc]
        else:
            fsrc_i = -1

        if fam_i >= 0:
            out_breadth = len(fam_out_species[fam_i])
            in_retained = len(fam_in_species[fam_i])
        else:
            out_breadth = len(out_present)
            in_retained = len(in_present)

        hit_genomes = sorted(tb_hits.get(pid, set()))
        chrom, start = _chrom_start(pid, src, gff3_paths, gff3_cache)

        out_rows.append([
            pid,
            shorts.index(src),
            round(out_frac, 4),
            out_breadth,
            in_retained,
            row.get('gene_name', '') or '',
            row.get('product_description', '') or '',
            fsrc_i,
            row.get('Best_Swissprot', '') or '',
            row.get('Pfam_Names', '') or '',
            row.get('Pfam_Accessions', '') or '',
            1 if hit_genomes else 0,
            ','.join(hit_genomes),
            fam_i,
            chrom,
            start,
        ])

    return {
        'project': project,
        'outgroup_min_frac': outgroup_min_frac,
        'loss_ingroup_max_frac': loss_ingroup_max_frac,
        'n_ingroup': len(ingroup_ids),
        'n_outgroup': len(outgroup_ids),
        'fields': LOSSES_ROW_FIELDS,
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
        'families': fam_index.payload(),
        'rows': out_rows,
    }
