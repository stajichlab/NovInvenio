import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Sample:
    group: str          # 'IN' or 'OUT'
    species: str
    strain: str
    protein: str        # filename, resolved relative to protein_dir
    dna: str            # filename, resolved relative to dna_dir (may be empty)
    short: str          # unique ≤8-char ID used throughout as the proteome key
    taxon_group: str


def parse_config(config_path: str | Path) -> list[Sample]:
    samples = []
    with open(config_path) as fh:
        for row in csv.DictReader(fh):
            samples.append(Sample(
                group=row['GROUP'].strip(),
                species=row['Species'].strip(),
                strain=row.get('Strain', '').strip(),
                protein=row['Protein'].strip(),
                dna=row.get('DNA', '').strip(),
                short=row['Short'].strip(),
                taxon_group=row['TaxonGroup'].strip(),
            ))
    return samples


def get_ingroup(samples: list[Sample]) -> list[Sample]:
    return [s for s in samples if s.group == 'IN']


def get_outgroup(samples: list[Sample]) -> list[Sample]:
    return [s for s in samples if s.group == 'OUT']


def short_to_group(samples: list[Sample]) -> dict[str, str]:
    """Map Short ID → GROUP (IN/OUT)."""
    return {s.short: s.group for s in samples}
