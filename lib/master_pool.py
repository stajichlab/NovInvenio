"""Master species pool: one row per species, already collapsed to its
representative assembly (see repr_assignments.tsv join in
bin/build_master_pool.py), with absolute Protein/DNA paths and a
fixed-width Lineage. See docs/superpowers/specs/2026-09-05-config-builder-design.md
("Data model" -> "Master pool").
"""
import csv
import re
from dataclasses import dataclass

from lineage import RANK_NAMES

MASTER_POOL_FIELDS = ['Species', 'Strain', 'ProteinPath', 'DNAPath', 'Lineage', 'NCBI_TaxID']


@dataclass
class MasterSample:
    species: str
    strain: str
    protein_path: str
    dna_path: str
    lineage: list[str]  # exactly len(RANK_NAMES); '' for an unrecorded rank
    ncbi_taxid: str


def load_master_pool(path) -> list[MasterSample]:
    samples = []
    with open(path, newline='') as fh:
        for row in csv.DictReader(fh):
            tokens = row['Lineage'].split(';')
            if len(tokens) != len(RANK_NAMES):
                raise SystemExit(
                    f"{path}: Lineage for {row['Species']!r} has {len(tokens)} fields, "
                    f"expected {len(RANK_NAMES)} ({','.join(RANK_NAMES)})"
                )
            samples.append(MasterSample(
                species=row['Species'].strip(),
                strain=row.get('Strain', '').strip(),
                protein_path=row['ProteinPath'].strip(),
                dna_path=row.get('DNAPath', '').strip(),
                lineage=[t.strip() for t in tokens],
                ncbi_taxid=row.get('NCBI_TaxID', '').strip(),
            ))

    seen = set()
    for s in samples:
        if s.species in seen:
            raise SystemExit(f"{path}: Duplicate Species {s.species!r}")
        seen.add(s.species)
    return samples


def make_short(species: str, used: set[str]) -> str:
    """Same shape as bin/convert_bfd_samples.py::make_short (genus[:4] +
    epithet[:4], numeric-suffix on collision) so Shorts look like existing
    configs. Determinism is the CALLER's responsibility -- see
    assign_shorts, which always processes the full pool in a fixed order.
    """
    parts = species.split()
    genus = re.sub(r'[^A-Za-z]', '', parts[0])[:4] if parts else 'Sp'
    epithet = re.sub(r'[^A-Za-z]', '', parts[1])[:4] if len(parts) > 1 else ''
    base = (genus[:1].upper() + genus[1:].lower() + epithet.lower())[:8] or 'Sp'
    short = base
    n = 2
    while short in used:
        short = f"{base[:8 - len(str(n))]}{n}"
        n += 1
    used.add(short)
    return short


def assign_shorts(samples: list[MasterSample]) -> dict[str, str]:
    """Species -> Short, deterministic: always sorts the given samples by
    Species name before assigning, so the same species set yields the same
    Shorts regardless of input order or which batch/subset is being
    rendered. Callers MUST pass the full master pool (not a per-batch
    subset) so a Short's meaning is stable across separate renders -- see
    the spec's "Short is not stored in the master pool" note.
    """
    used: set[str] = set()
    mapping = {}
    for s in sorted(samples, key=lambda s: s.species):
        mapping[s.species] = make_short(s.species, used)
    return mapping
