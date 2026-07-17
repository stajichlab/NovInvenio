#!/usr/bin/env python3
"""
Generate a per-run analysis config CSV (GROUP,Species,Strain,Protein,DNA,Short,
TaxonGroup) by selecting an ingroup clade and a sampled outgroup from a master
samples file.

The master samples file (configs/samples.csv) lists every available proteome
once, with a semicolon-delimited Lineage string per species instead of a
single flat TaxonGroup (Species,Strain,Protein,DNA,Short,Lineage). This script
slices out a run-specific config from that pool by matching lineage tokens —
it is taxonomy-string matching, not phylogenetic inference. If a real species
tree becomes available, a --tree mode (MRCA/sister-clade selection via ete3 or
dendropy) would be the natural extension, but is not implemented here.

Ingroup selection:
  --ingroup-taxon NAME   All samples whose Lineage contains this token.
  --ingroup-short SHORT  Pin exact Short IDs into the ingroup (repeatable),
                         for focal-species configs. At least one of
                         --ingroup-taxon / --ingroup-short is required.

Outgroup selection:
  --outgroup-taxon NAME  Candidate outgroup pool (repeatable). If omitted,
                         defaults to "siblings of --ingroup-taxon's parent
                         clade" — every other child of the ingroup's parent
                         lineage segment.
  --exclude-taxon NAME   Drop any sample matching this token (repeatable).
  --max-per-outgroup-taxon N
                         Stratified cap: keep at most N species per matched
                         outgroup clade, so no single lineage dominates.
                         Deterministic (sorted by Short) unless --random.
"""
import argparse
import csv
import random
import sys
from dataclasses import dataclass
from pathlib import Path

PROTEIN_SUBDIRS = ['pep', 'proteins']
DNA_SUBDIRS = ['dna', 'genome', 'scaffolds']


@dataclass
class MasterSample:
    species: str
    strain: str
    protein: str
    dna: str
    short: str
    lineage: list[str]


def load_samples(path):
    samples = []
    with open(path, newline='') as fh:
        for row in csv.DictReader(fh):
            samples.append(MasterSample(
                species=row['Species'].strip(),
                strain=row.get('Strain', '').strip(),
                protein=row['Protein'].strip(),
                dna=row.get('DNA', '').strip(),
                short=row['Short'].strip(),
                lineage=[tok.strip() for tok in row['Lineage'].split(';') if tok.strip()],
            ))

    seen = {}
    for s in samples:
        if s.short in seen:
            raise SystemExit(
                f"Duplicate Short '{s.short}' in {path} "
                f"(rows: {seen[s.short].species!r}, {s.species!r})"
            )
        seen[s.short] = s
    return samples


def matches_taxon(sample, token):
    return token.lower() in {t.lower() for t in sample.lineage}


def parent_segment(sample, token):
    """Lineage segment immediately before `token`, or None if not found / first."""
    lowered = [t.lower() for t in sample.lineage]
    tok = token.lower()
    if tok not in lowered:
        return None
    idx = lowered.index(tok)
    return sample.lineage[idx - 1] if idx > 0 else None


def child_segment(sample, parent_token):
    """Lineage segment immediately after `parent_token` in this sample's own lineage."""
    lowered = [t.lower() for t in sample.lineage]
    tok = parent_token.lower()
    if tok not in lowered:
        return None
    idx = lowered.index(tok)
    if idx + 1 < len(sample.lineage):
        return sample.lineage[idx + 1]
    return sample.lineage[idx]


def fallback_taxon_group(sample):
    """Deepest-but-one lineage segment, used when a sample has no taxon match
    (e.g. an --ingroup-short pin with no --ingroup-taxon)."""
    if len(sample.lineage) >= 2:
        return sample.lineage[-2]
    if sample.lineage:
        return sample.lineage[-1]
    return ''


def select_ingroup(samples, ingroup_taxon, ingroup_shorts):
    by_short = {s.short: s for s in samples}
    for short in ingroup_shorts:
        if short not in by_short:
            raise SystemExit(f"--ingroup-short '{short}' not found in samples file")

    picked = {}
    if ingroup_taxon:
        for s in samples:
            if matches_taxon(s, ingroup_taxon):
                picked[s.short] = (s, ingroup_taxon)
    for short in ingroup_shorts:
        s = by_short[short]
        picked.setdefault(s.short, (s, fallback_taxon_group(s)))
    return picked


def select_outgroup_pool(samples, ingroup_taxon, outgroup_taxa, exclude_taxa, ingroup_picked):
    pool = {}
    if outgroup_taxa:
        for token in outgroup_taxa:
            for s in samples:
                if s.short in ingroup_picked or s.short in pool:
                    continue
                if matches_taxon(s, token):
                    pool[s.short] = (s, token)
    else:
        if not ingroup_taxon:
            raise SystemExit(
                "no --outgroup-taxon given and no --ingroup-taxon to derive a "
                "sibling default from — pass --outgroup-taxon explicitly"
            )
        ref = next((s for s in ingroup_picked.values() if matches_taxon(s[0], ingroup_taxon)), None)
        if ref is None:
            raise SystemExit(f"no ingroup sample matched --ingroup-taxon '{ingroup_taxon}'")
        parent = parent_segment(ref[0], ingroup_taxon)
        if parent is None:
            raise SystemExit(
                f"--ingroup-taxon '{ingroup_taxon}' has no parent lineage segment to "
                "derive sibling outgroups from — pass --outgroup-taxon explicitly"
            )
        for s in samples:
            if s.short in ingroup_picked:
                continue
            if matches_taxon(s, ingroup_taxon):
                continue
            if parent.lower() in {t.lower() for t in s.lineage}:
                sibling = child_segment(s, parent)
                pool[s.short] = (s, sibling or parent)

    for token in exclude_taxa:
        pool = {short: v for short, v in pool.items() if not matches_taxon(v[0], token)}

    return pool


def apply_stratified_cap(pool, max_per_taxon, use_random, seed):
    if max_per_taxon is None:
        return pool

    by_group: dict[str, list] = {}
    for short, (sample, group) in pool.items():
        by_group.setdefault(group, []).append((short, sample, group))

    rng = random.Random(seed) if use_random else None
    kept = {}
    for group, entries in by_group.items():
        entries = sorted(entries, key=lambda e: e[0])
        if use_random:
            rng.shuffle(entries)
        for short, sample, grp in entries[:max_per_taxon]:
            kept[short] = (sample, grp)
    return kept


def resolve_under(data_dir, filename, subdirs):
    if not filename:
        return True
    candidates = [Path(data_dir) / filename] + [Path(data_dir) / d / filename for d in subdirs]
    return any(c.exists() for c in candidates)


def validate_data_dir(rows, data_dir):
    missing = []
    for row in rows:
        if not resolve_under(data_dir, row['Protein'], PROTEIN_SUBDIRS):
            missing.append(f"{row['Short']}: Protein '{row['Protein']}' not found under {data_dir}")
        if row['DNA'] and not resolve_under(data_dir, row['DNA'], DNA_SUBDIRS):
            missing.append(f"{row['Short']}: DNA '{row['DNA']}' not found under {data_dir}")
    if missing:
        raise SystemExit("Missing files under --data-dir:\n  " + "\n  ".join(missing))


def build_rows(ingroup_picked, outgroup_picked):
    rows = []
    seen_short = set()
    for short, (s, taxon_group) in sorted(ingroup_picked.items()):
        seen_short.add(short)
        rows.append({
            'GROUP': 'IN', 'Species': s.species, 'Strain': s.strain,
            'Protein': s.protein, 'DNA': s.dna, 'Short': s.short,
            'TaxonGroup': taxon_group,
        })
    for short, (s, taxon_group) in sorted(outgroup_picked.items()):
        if short in seen_short:
            raise SystemExit(f"Short '{short}' selected in both ingroup and outgroup")
        rows.append({
            'GROUP': 'OUT', 'Species': s.species, 'Strain': s.strain,
            'Protein': s.protein, 'DNA': s.dna, 'Short': s.short,
            'TaxonGroup': taxon_group,
        })
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--samples', required=True, help='Master samples CSV (Species,Strain,Protein,DNA,Short,Lineage)')
    p.add_argument('--ingroup-taxon', help='Lineage token defining the ingroup clade')
    p.add_argument('--ingroup-short', action='append', default=[], help='Pin exact Short ID into the ingroup (repeatable)')
    p.add_argument('--outgroup-taxon', action='append', default=[], help='Lineage token for the outgroup pool (repeatable)')
    p.add_argument('--exclude-taxon', action='append', default=[], help='Lineage token to drop from consideration (repeatable)')
    p.add_argument('--max-per-outgroup-taxon', type=int, default=None, help='Stratified cap: max outgroup species per matched clade')
    p.add_argument('--random', action='store_true', help='Randomly subsample within each stratified group instead of deterministic sort')
    p.add_argument('--seed', type=int, default=None, help='RNG seed (only used with --random)')
    p.add_argument('--data-dir', help='If given, verify Protein/DNA files resolve under this directory')
    p.add_argument('--output', required=True, help='Destination config CSV')
    args = p.parse_args()

    if not args.ingroup_taxon and not args.ingroup_short:
        raise SystemExit("must pass at least one of --ingroup-taxon or --ingroup-short")

    samples = load_samples(args.samples)

    ingroup_picked = select_ingroup(samples, args.ingroup_taxon, args.ingroup_short)
    if not ingroup_picked:
        raise SystemExit(f"no species matched --ingroup-taxon '{args.ingroup_taxon}' or --ingroup-short — check spelling against {args.samples}")

    outgroup_pool = select_outgroup_pool(
        samples, args.ingroup_taxon, args.outgroup_taxon, args.exclude_taxon, ingroup_picked
    )
    if not outgroup_pool:
        raise SystemExit("no species matched the outgroup selection — check --outgroup-taxon/--exclude-taxon against " + args.samples)

    outgroup_picked = apply_stratified_cap(outgroup_pool, args.max_per_outgroup_taxon, args.random, args.seed)

    rows = build_rows(ingroup_picked, outgroup_picked)

    if args.data_dir:
        validate_data_dir(rows, args.data_dir)

    with open(args.output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['GROUP', 'Species', 'Strain', 'Protein', 'DNA', 'Short', 'TaxonGroup'])
        writer.writeheader()
        writer.writerows(rows)

    n_in = sum(1 for r in rows if r['GROUP'] == 'IN')
    n_out = sum(1 for r in rows if r['GROUP'] == 'OUT')
    print(f"Wrote {args.output}: {n_in} ingroup, {n_out} outgroup", file=sys.stderr)


if __name__ == '__main__':
    main()
