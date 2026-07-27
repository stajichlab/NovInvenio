import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Union


@dataclass
class Sample:
    group: str          # 'IN', 'OUT', 'TARGET', 'DISC_OUT', 'NEAR_IN', 'BROAD_OUT'
    species: str
    strain: str
    protein: str        # filename, resolved relative to protein_dir
    dna: str            # filename, resolved relative to dna_dir (may be empty)
    short: str          # unique ≤8-char ID used throughout as the proteome key
    taxon_group: str


# All recognised group labels. IN/OUT are the classic pairwise/mmseqs pathway
# roles; TARGET/DISC_OUT/NEAR_IN/BROAD_OUT are used by the novelty_discovery /
# novelty_screen workflow.
GROUPS = {'IN', 'OUT', 'TARGET', 'DISC_OUT', 'NEAR_IN', 'BROAD_OUT'}


def parse_config(config_path: Union[str, Path]) -> list[Sample]:
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


def get_group(samples: list[Sample], group: str) -> list[Sample]:
    """Return all samples belonging to the given GROUP label."""
    return [s for s in samples if s.group == group]


def get_ingroup(samples: list[Sample]) -> list[Sample]:
    return get_group(samples, 'IN')


def get_outgroup(samples: list[Sample]) -> list[Sample]:
    return get_group(samples, 'OUT')


def get_target(samples: list[Sample]) -> list[Sample]:
    """Return TARGET samples for novelty_discovery."""
    return get_group(samples, 'TARGET')


def get_disc_out(samples: list[Sample]) -> list[Sample]:
    """Return DISC_OUT samples for novelty_discovery."""
    return get_group(samples, 'DISC_OUT')


def get_near_in(samples: list[Sample]) -> list[Sample]:
    """Return NEAR_IN samples for novelty_screen."""
    return get_group(samples, 'NEAR_IN')


def get_broad_out(samples: list[Sample]) -> list[Sample]:
    """Return BROAD_OUT samples for novelty_screen."""
    return get_group(samples, 'BROAD_OUT')


def short_to_group(samples: list[Sample]) -> dict[str, str]:
    """Map Short ID → GROUP (IN/OUT/TARGET/DISC_OUT/NEAR_IN/BROAD_OUT)."""
    return {s.short: s.group for s in samples}
