"""Species -> SourceDB (per-gene database linkout) loader.

SourceDB is the config CSV's existing optional, report-only column (see
lib/report_common.py's genomeDbLink()): "fungidb", "mycocosm:<portal>",
"ensemblfungi:<species>", "veupathdb:<project>", "ncbipep", or a URL
template containing "{gene}". This module loads a small hand-curated seed file
(config_support/source_db.csv: Species,SourceDB,notes) mapping a species
name to one of those values, for the targeted config-builder renderer to
merge in at render time -- mirroring how config_support/traits/ is loaded
and merged, since neither is derivable from Fungi_BFD/repr_assignments.tsv
data alone.

A SourceDB value is only semi-trusted (see CLAUDE.md's "The interactive
report" section): the report layer's own genomeDbLink() already requires
an http(s) scheme before building a link from a raw URL template, but this
loader hard-errors on an obviously malformed value at load time too, so a
data-entry mistake in the seed file is caught here rather than silently
reaching a rendered config.
"""
import csv

_KNOWN_PREFIXES = ('mycocosm:', 'ensemblfungi:', 'veupathdb:')


def _is_valid(value: str) -> bool:
    if value in ('fungidb', 'ncbipep'):
        return True
    if value.startswith(_KNOWN_PREFIXES) and len(value.split(':', 1)[1]) > 0:
        return True
    if '{gene}' in value and value.startswith(('http://', 'https://')):
        return True
    return False


def load_source_db(path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with open(path, newline='') as fh:
        for row in csv.DictReader(fh):
            species = row['Species'].strip()
            value = row['SourceDB'].strip()
            if not _is_valid(value):
                raise SystemExit(
                    f"{path}: {species!r} has an invalid SourceDB value {value!r} -- "
                    "must be 'fungidb', 'ncbipep', 'mycocosm:<portal>', "
                    "'ensemblfungi:<species>', 'veupathdb:<project>', or an http(s) "
                    "URL template containing '{gene}'"
                )
            mapping[species] = value
    return mapping
