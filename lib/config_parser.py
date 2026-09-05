import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Union


@dataclass
class Sample:
    group: str          # 'IN', 'OUT', 'DISCOVERY_TARGET', 'DISCOVERY_OUT', 'NEAR_INGROUP',
                         # 'BROAD_OUTGROUP' -- always normalized to these canonical spellings
                         # by parse_config(), even if the source CSV used an old alias.
    species: str
    strain: str
    protein: str        # filename, resolved relative to protein_dir
    dna: str            # filename, resolved relative to dna_dir (may be empty)
    short: str          # unique ≤8-char ID used throughout as the proteome key
    taxon_group: str
    gff3: str = ''      # optional filename, resolved relative to data_dir (also checked
                         # under pep/, dna/, genome/, scaffolds/, and gff3/ subdirs) --
                         # empty means no chrom/start data for this species' report rows
    source_db: str = ''  # optional genome-database key used to build a per-gene linkout in
                         # the HTML reports. Accepted forms: 'fungidb',
                         # 'mycocosm:<portal>', 'ensemblfungi:<species>',
                         # 'veupathdb:<project>', or any URL template containing '{gene}'.
                         # Empty means the report falls back to a model-organism FungiDB
                         # link (if the annotation came from one) or an NCBI Protein search.
                         # Report-only: never affects any presence/novelty call.
    taxid: str = ''      # optional NCBI taxonomy ID. Report-only -- turns the species name
                         # in a report's detail panel into an NCBI Taxonomy link; without it
                         # the link falls back to a by-name search.


# All recognised group labels. IN/OUT are the classic pairwise/mmseqs pathway roles;
# DISCOVERY_TARGET/DISCOVERY_OUT/NEAR_INGROUP/BROAD_OUTGROUP are used by the
# novelty_discovery/novelty_screen workflow (todo/novelty-discovery-screen.md).
GROUPS = {'IN', 'OUT', 'DISCOVERY_TARGET', 'DISCOVERY_OUT', 'NEAR_INGROUP', 'BROAD_OUTGROUP'}

# Original novelty_discovery/novelty_screen labels (issues #24-#29), renamed for clarity
# (todo/rename-novelty-discovery-group-labels.md) -- still accepted in a config CSV's GROUP
# column and normalized to the canonical spelling above by parse_config(), so existing
# configs keep working unchanged.
GROUP_ALIASES = {
    'TARGET': 'DISCOVERY_TARGET',
    'DISC_OUT': 'DISCOVERY_OUT',
    'NEAR_IN': 'NEAR_INGROUP',
    'BROAD_OUT': 'BROAD_OUTGROUP',
}

# Coarse ingroup/outgroup banding shared by every downstream consumer that only
# cares "is this the query side or the reference side" (report payload builders,
# make_novelties.py) rather than the fine-grained role (used by the workflows
# themselves via get_discovery_target()/get_discovery_out()/etc. above).
# 'NEAR_INGROUP'/'BROAD_OUTGROUP' fold into the outgroup band here; the finer
# near/broad distinction is carried by the novelty_category column (issue #28),
# not a third band.
INGROUP_ROLES = {'IN', 'DISCOVERY_TARGET'}
OUTGROUP_ROLES = {'OUT', 'DISCOVERY_OUT', 'NEAR_INGROUP', 'BROAD_OUTGROUP'}


def parse_config(config_path: Union[str, Path]) -> list[Sample]:
    samples = []
    with open(config_path) as fh:
        for row in csv.DictReader(fh):
            group = row['GROUP'].strip()
            group = GROUP_ALIASES.get(group, group)
            samples.append(Sample(
                group=group,
                species=row['Species'].strip(),
                strain=row.get('Strain', '').strip(),
                protein=row['Protein'].strip(),
                dna=row.get('DNA', '').strip(),
                short=row['Short'].strip(),
                taxon_group=row['TaxonGroup'].strip(),
                gff3=(row.get('GFF3') or '').strip(),
                source_db=(row.get('SourceDB') or '').strip(),
                taxid=(row.get('NCBI_TaxID') or '').strip(),
            ))
    return samples


def get_group(samples: list[Sample], group: str) -> list[Sample]:
    """Return all samples belonging to the given GROUP label."""
    return [s for s in samples if s.group == group]


def get_ingroup(samples: list[Sample]) -> list[Sample]:
    return get_group(samples, 'IN')


def get_outgroup(samples: list[Sample]) -> list[Sample]:
    return get_group(samples, 'OUT')


def get_discovery_target(samples: list[Sample]) -> list[Sample]:
    """Return DISCOVERY_TARGET samples for novelty_discovery."""
    return get_group(samples, 'DISCOVERY_TARGET')


def get_discovery_out(samples: list[Sample]) -> list[Sample]:
    """Return DISCOVERY_OUT samples for novelty_discovery."""
    return get_group(samples, 'DISCOVERY_OUT')


def get_near_ingroup(samples: list[Sample]) -> list[Sample]:
    """Return NEAR_INGROUP samples for novelty_screen."""
    return get_group(samples, 'NEAR_INGROUP')


def get_broad_outgroup(samples: list[Sample]) -> list[Sample]:
    """Return BROAD_OUTGROUP samples for novelty_screen."""
    return get_group(samples, 'BROAD_OUTGROUP')


def short_to_group(samples: list[Sample]) -> dict[str, str]:
    """Map Short ID → GROUP (IN/OUT/DISCOVERY_TARGET/DISCOVERY_OUT/NEAR_INGROUP/BROAD_OUTGROUP)."""
    return {s.short: s.group for s in samples}
