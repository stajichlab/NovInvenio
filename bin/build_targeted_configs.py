#!/usr/bin/env python3
"""
Render one or more NovInvenio analysis config CSVs from a small YAML batch
spec: a focal species per study, an ingroup-companion picker (mode:
nearest/trait/explicit), and a named, reusable outgroup pool. See
docs/superpowers/specs/2026-09-05-config-builder-design.md ("Selection spec
format", "Renderer") for the full design.

Usage:
    bin/build_targeted_configs.py \\
        --master-pool config_support/master_pool.csv \\
        --trait-definitions config_support/traits/trait_definitions.yaml \\
        --traits config_support/traits/traits.csv \\
        --batch-spec configs/batches/mucoromycota_focal_v1.yaml \\
        --outdir configs/ \\
        --link-dir data/mucoromycota_focal_v1 \\
        --source-db config_support/source_db.csv

--link-dir is what makes the rendered config directly usable with main.nf's
--data_dir: without it, Protein/DNA columns carry the master pool's absolute
paths, which main.nf's resolve_fa() cannot resolve (see _row()'s docstring).

--source-db populates the config's optional SourceDB column (per-gene database
linkout) from a small Species,SourceDB seed CSV -- see lib/source_db.py. A
species not in the seed file gets an empty SourceDB, not an error.
"""
import argparse
import csv
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from lineage import RANK_NAMES  # noqa: E402
from master_pool import assign_shorts, load_master_pool  # noqa: E402
from targeted_selection import (  # noqa: E402
    DEFAULT_SCOPE_RANK,
    select_nearest,
    select_trait,
)
from source_db import load_source_db  # noqa: E402
from trait_data import load_trait_definitions, load_traits  # noqa: E402

CONFIG_FIELDS = ['GROUP', 'Species', 'Strain', 'Protein', 'DNA', 'Short', 'TaxonGroup', 'NCBI_TaxID', 'SourceDB']


def _by_species(pool):
    return {s.species: s for s in pool}


def _resolve_members(names, by_species, context):
    missing = [n for n in names if n not in by_species]
    if missing:
        raise SystemExit(f"{context}: species not found in master pool: {missing}")
    return names


def _link_into(link_dir, subdir, src_path):
    """Symlink src_path into <link_dir>/<subdir>/<basename>, creating the
    subdir if needed and skipping if a symlink with that name already
    exists -- mirrors bin/convert_bfd_samples.py's --link-dir pattern.
    Returns the basename to write into the config's Protein/DNA column.
    """
    src = Path(src_path)
    target_dir = Path(link_dir) / subdir
    target_dir.mkdir(parents=True, exist_ok=True)
    link = target_dir / src.name
    if not link.exists():
        link.symlink_to(src.resolve())
    return src.name


def _row(group, sample, short, taxon_group, link_dir=None, source_db_map=None):
    # NOTE: without --link-dir, ProteinPath/DNAPath (absolute paths from the
    # master pool) are written straight into the Protein/DNA columns -- this
    # keeps existing fixture-based tests passing without needing a link dir,
    # but the resulting config is NOT main.nf-ready: main.nf's resolve_fa()
    # builds candidates as file("${params.data_dir}/${basename}"), so an
    # absolute path here resolves to a nonexistent path. Pass --link-dir to
    # get a config with basenames plus a small per-batch symlink directory
    # main.nf's --data_dir can actually point at.
    if link_dir is not None:
        protein = _link_into(link_dir, 'pep', sample.protein_path)
        dna = _link_into(link_dir, 'dna', sample.dna_path)
    else:
        protein = sample.protein_path
        dna = sample.dna_path
    return {
        'GROUP': group, 'Species': sample.species, 'Strain': sample.strain,
        'Protein': protein, 'DNA': dna, 'Short': short,
        'TaxonGroup': taxon_group, 'NCBI_TaxID': sample.ncbi_taxid,
        'SourceDB': (source_db_map or {}).get(sample.species, ''),
    }


def render_batch(master_pool_path, trait_definitions_path, traits_path, batch_spec_path, outdir, link_dir=None, source_db_path=None) -> list[dict]:
    pool = load_master_pool(master_pool_path)
    by_species = _by_species(pool)
    short_map = assign_shorts(pool)

    definitions = load_trait_definitions(trait_definitions_path)
    traits_by_species = load_traits(traits_path, definitions)
    source_db_map = load_source_db(source_db_path) if source_db_path else {}

    with open(batch_spec_path) as fh:
        spec = yaml.safe_load(fh)

    batch_name = spec['batch']
    outgroup_pools = {
        name: _resolve_members(cfg['members'], by_species, f"outgroup_pools.{name}")
        for name, cfg in spec.get('outgroup_pools', {}).items()
    }

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summaries = []

    for study in spec['studies']:
        focal_name = study['focal']
        if focal_name not in by_species:
            raise SystemExit(f"study focal {focal_name!r} not found in master pool")
        focal = by_species[focal_name]

        pool_name = study['outgroup_pool']
        if pool_name not in outgroup_pools:
            raise SystemExit(f"study for {focal_name!r}: outgroup_pool {pool_name!r} not defined")
        outgroup_members = outgroup_pools[pool_name]
        if focal_name in outgroup_members:
            raise SystemExit(
                f"study for {focal_name!r}: focal species is inside its own outgroup pool {pool_name!r}"
            )
        outgroup_set = set(outgroup_members)

        extra = study['ingroup_extra']
        mode = extra['mode']
        n = extra.get('n')
        scope_rank = extra.get('scope_rank', DEFAULT_SCOPE_RANK)

        if mode == 'nearest':
            # excluded= is passed straight through so exclusion happens BEFORE
            # ranking/truncation inside select_nearest -- a candidate excluded
            # by the outgroup pool is backfilled by the next-best candidate,
            # not just dropped from an already-truncated top-n list.
            candidates = select_nearest(
                focal.species, focal.lineage, pool, n=n, scope_rank=scope_rank, excluded=outgroup_set,
            )
            companions = [
                {'species': c.species, 'taxon_group': by_species[c.species].lineage[RANK_NAMES.index(c.rank_name)], 'reason': ''}
                for c in candidates
            ]
        elif mode == 'trait':
            candidates = select_trait(
                focal.species, focal.lineage, pool, extra['trait'], extra['value'], n=n,
                traits_by_species=traits_by_species, scope_rank=scope_rank, excluded=outgroup_set,
            )
            companions = [
                {'species': c.species, 'taxon_group': by_species[c.species].lineage[RANK_NAMES.index(c.rank_name)], 'reason': ''}
                for c in candidates
            ]
        elif mode == 'explicit':
            members = _resolve_members(extra['members'], by_species, f"study for {focal_name!r}")
            overlap = set(members) & outgroup_set
            if overlap:
                raise SystemExit(f"study for {focal_name!r}: explicit member(s) {overlap} also in outgroup pool {pool_name!r}")
            reason = extra.get('reason', '')
            companions = [
                {'species': m, 'taxon_group': [t for t in by_species[m].lineage if t][-1] if any(by_species[m].lineage) else '', 'reason': reason}
                for m in members
            ]
        else:
            raise SystemExit(f"study for {focal_name!r}: unknown mode {mode!r} (must be nearest/trait/explicit)")

        rows = [_row('IN', focal, short_map[focal.species], '', link_dir=link_dir, source_db_map=source_db_map)]
        for c in companions:
            rows.append(_row(
                'IN', by_species[c['species']], short_map[c['species']], c['taxon_group'],
                link_dir=link_dir, source_db_map=source_db_map,
            ))
        for m in outgroup_members:
            s = by_species[m]
            taxon_group = [t for t in s.lineage if t][-1] if any(s.lineage) else ''
            rows.append(_row('OUT', s, short_map[m], taxon_group, link_dir=link_dir, source_db_map=source_db_map))

        config_path = outdir / f"{short_map[focal.species]}_{batch_name}.csv"
        with open(config_path, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=CONFIG_FIELDS)
            writer.writeheader()
            writer.writerows(rows)

        map_path = outdir / f"{short_map[focal.species]}_{batch_name}.map.tsv"
        with open(map_path, 'w') as fh:
            for r in rows:
                fh.write(f"{r['Species']}\t{r['Short']}\n")

        summaries.append({
            'batch': batch_name,
            'focal': focal_name,
            'config_path': str(config_path),
            'map_path': str(map_path),
            'companions': companions,
            'outgroup_pool': pool_name,
            'outgroup_size': len(outgroup_members),
        })

    return summaries


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--master-pool', required=True)
    p.add_argument('--trait-definitions', required=True)
    p.add_argument('--traits', required=True)
    p.add_argument('--batch-spec', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument(
        '--link-dir',
        help=(
            'If given, symlink each selected species\' Protein/DNA files into '
            '<link-dir>/pep/ and <link-dir>/dna/ and write basenames (not absolute '
            'paths) into the rendered config\'s Protein/DNA columns, so the config is '
            'directly usable as --data_dir <link-dir> with main.nf. Omit to keep '
            'writing absolute master-pool paths (not main.nf-ready without this flag).'
        ),
    )
    p.add_argument(
        '--source-db',
        help=(
            'Optional path to a Species,SourceDB seed CSV (e.g. config_support/source_db.csv) '
            'populating the rendered config\'s optional SourceDB column -- see '
            'lib/report_common.py\'s genomeDbLink() for the accepted value forms. '
            'A species not listed gets an empty SourceDB (falls back to the report\'s '
            'default NCBI Protein search / remote-homology cluster linkout).'
        ),
    )
    args = p.parse_args()

    summaries = render_batch(
        args.master_pool, args.trait_definitions, args.traits, args.batch_spec, args.outdir,
        link_dir=args.link_dir, source_db_path=args.source_db,
    )

    for s in summaries:
        print(f"\n{s['focal']}  ({s['config_path']})", file=sys.stderr)
        for c in s['companions']:
            reason = f" -- {c['reason']}" if c['reason'] else ''
            print(f"  + {c['species']}  [{c['taxon_group']}]{reason}", file=sys.stderr)
        print(f"  outgroup: {s['outgroup_pool']} ({s['outgroup_size']} species)", file=sys.stderr)


if __name__ == '__main__':
    main()
