# Targeted Config-Builder Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the targeted config-builder tool designed in `docs/superpowers/specs/2026-09-05-config-builder-design.md`: a master species pool (BFD samples + representative-assembly join), a hand-curated trait layer, a lineage-proximity/trait-based ingroup-companion picker, and a renderer that turns a small YAML batch spec into one or more NovInvenio-ready config CSVs.

**Architecture:** Five small `lib/` modules, each independently testable, composed by one `bin/` CLI script per pipeline stage (pool building, then rendering). No changes to `main.nf` or any Nextflow workflow — the renderer only ever emits the same `GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup` CSV contract `main.nf` already consumes.

**Tech Stack:** Python 3.12, `pyyaml` (already a pixi dependency), `pytest`. No new dependencies.

**Spec:** `docs/superpowers/specs/2026-09-05-config-builder-design.md` — read it in full before touching any task below. Its design decisions (Species-name keying, prefix-based rank-aligned lineage proximity, the `repr_assignments.tsv` join semantics, the `nearest`/`trait`/`explicit` mode algorithm, disjointness rules, IN/OUT-only GROUP scheme, YAML batch format) are settled; this plan implements them and does not revisit them.

## Global Constraints

- Every `bin/` script is `#!/usr/bin/env python3`, executable (`chmod +x`), uses `argparse` (no positional args), and imports shared logic via `sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))` then bare `import <module>` — matches `bin/make_config.py`'s existing pattern.
- Every new shared module lives in `lib/` as a bare-importable module (no package-relative imports) — matches `lib/config_parser.py`.
- Tests live in `tests/test_<module>.py`; `tests/conftest.py` already puts `lib/` on `sys.path`, so test files import shared modules the same bare way (`from lineage import ...`), no per-file `sys.path` hacking needed.
- The master pool's `Lineage` column is **fixed-width**: exactly 7 semicolon-separated fields in the order `PHYLUM;SUBPHYLUM;CLASS;SUBCLASS;ORDER;FAMILY;GENUS`, with an **empty string kept in place** for a rank that isn't recorded (never dropped) — this is what makes position-based rank alignment safe, and is a deliberate departure from `bin/convert_bfd_samples.py`'s existing behavior (which drops empty ranks and is not touched by this plan).
- Every hard-error path in the spec (undeclared trait/value, `none` coexisting with another value, zero/multiple `is_representative` rows, focal-in-own-outgroup-pool, empty `mode: trait` filter result) raises `SystemExit` with a message naming the offending species/row — never silently skips or falls back.

---

### Task 1: Lineage proximity ranking

**Files:**
- Create: `lib/lineage.py`
- Test: `tests/test_lineage.py`

**Interfaces:**
- Produces: `RANK_NAMES: list[str]` (7 names, `PHYLUM` first, `GENUS` last), `lineage_match(lineage_a: list[str], lineage_b: list[str]) -> str` (deepest shared rank *name*, `''` if none share).

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_lineage.py
from lineage import RANK_NAMES, lineage_match


def test_rank_names_fixed_order():
    assert RANK_NAMES == [
        'PHYLUM', 'SUBPHYLUM', 'CLASS', 'SUBCLASS', 'ORDER', 'FAMILY', 'GENUS',
    ]


def test_exact_match_through_genus():
    a = ['Mucoromycota', 'Mucoromycotina', 'Mucoromycetes', '', 'Mucorales', 'Mucoraceae', 'Mucor']
    b = ['Mucoromycota', 'Mucoromycotina', 'Mucoromycetes', '', 'Mucorales', 'Mucoraceae', 'Mucor']
    assert lineage_match(a, b) == 'GENUS'


def test_missing_rank_does_not_break_a_deeper_match():
    # a has no SUBPHYLUM recorded, b does -- must not be misjudged as
    # diverging at SUBPHYLUM; they still agree all the way to FAMILY.
    a = ['Mucoromycota', '', 'Mucoromycetes', '', 'Mucorales', 'Mucoraceae', 'Mucor']
    b = ['Mucoromycota', 'Mucoromycotina', 'Mucoromycetes', '', 'Mucorales', 'Mucoraceae', 'Rhizopus']
    assert lineage_match(a, b) == 'FAMILY'


def test_true_divergence_stops_at_last_agreeing_rank():
    a = ['Mucoromycota', '', 'Mucoromycetes', '', 'Mucorales', 'Mucoraceae', 'Mucor']
    b = ['Mucoromycota', '', 'Mucoromycetes', '', 'Entomophthorales', 'Ancylistaceae', 'Conidiobolus']
    assert lineage_match(a, b) == 'CLASS'


def test_no_shared_rank():
    a = ['Mucoromycota', '', 'Mucoromycetes', '', 'Mucorales', 'Mucoraceae', 'Mucor']
    b = ['Chytridiomycota', '', 'Chytridiomycetes', '', 'Spizellomycetales', 'Spizellomycetaceae', 'Spizellomyces']
    assert lineage_match(a, b) == ''


def test_case_insensitive():
    a = ['mucoromycota', '', '', '', '', '', '']
    b = ['MUCOROMYCOTA', '', '', '', '', '', '']
    assert lineage_match(a, b) == 'PHYLUM'
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest tests/test_lineage.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'lineage'`

- [ ] **Step 3: Write the implementation**

```python
# lib/lineage.py
"""Lineage-proximity ranking shared by the master-pool builder and the
targeted config renderer. See docs/superpowers/specs/2026-09-05-config-builder-design.md
("Lineage-proximity ranking") for the full rationale.
"""

RANK_NAMES = ['PHYLUM', 'SUBPHYLUM', 'CLASS', 'SUBCLASS', 'ORDER', 'FAMILY', 'GENUS']


def lineage_match(lineage_a: list[str], lineage_b: list[str]) -> str:
    """Deepest rank NAME both lineages still agree on.

    Both lineages must be exactly len(RANK_NAMES) long (one slot per rank,
    '' for a rank that isn't recorded). A rank where either side is '' is
    skipped -- it neither extends nor breaks the match -- so a missing
    SUBPHYLUM/SUBCLASS never causes a same-genus pair to be misjudged as
    merely same-order. Comparison is case-insensitive. Returns '' if no
    rank matches (including if they diverge at PHYLUM itself).
    """
    if len(lineage_a) != len(RANK_NAMES) or len(lineage_b) != len(RANK_NAMES):
        raise ValueError(
            f"lineage must have exactly {len(RANK_NAMES)} fields "
            f"(got {len(lineage_a)} and {len(lineage_b)})"
        )
    deepest = ''
    for name, a, b in zip(RANK_NAMES, lineage_a, lineage_b):
        if not a or not b:
            continue
        if a.lower() != b.lower():
            break
        deepest = name
    return deepest
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi run pytest tests/test_lineage.py -v`
Expected: PASS (6 tests)

- [ ] **Step 5: Commit**

```bash
git add lib/lineage.py tests/test_lineage.py
git commit -m "feat: add prefix-based, rank-name-aligned lineage proximity matching"
```

---

### Task 2: Master pool data structure and loader

**Files:**
- Create: `lib/master_pool.py`
- Test: `tests/test_master_pool.py`

**Interfaces:**
- Consumes: `RANK_NAMES` from `lineage` (Task 1).
- Produces: `MasterSample` dataclass (`species: str`, `strain: str`, `protein_path: str`, `dna_path: str`, `lineage: list[str]`, `ncbi_taxid: str`), `MASTER_POOL_FIELDS: list[str]` (`['Species', 'Strain', 'ProteinPath', 'DNAPath', 'Lineage', 'NCBI_TaxID']`), `load_master_pool(path) -> list[MasterSample]`, `make_short(species: str, used: set[str]) -> str`, `assign_shorts(samples: list[MasterSample]) -> dict[str, str]` (Species → Short).

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_master_pool.py
import pytest
from master_pool import (
    MASTER_POOL_FIELDS,
    MasterSample,
    assign_shorts,
    load_master_pool,
    make_short,
)

POOL_CSV = """\
Species,Strain,ProteinPath,DNAPath,Lineage,NCBI_TaxID
Mucor circinelloides,1006PhL,/data/Mucor_circinelloides.pep.fa,/data/Mucor_circinelloides.dna.fa,Mucoromycota;Mucoromycotina;Mucoromycetes;;Mucorales;Mucoraceae;Mucor,36698
Rhizopus arrhizus,,/data/Rhizopus_arrhizus.pep.fa,/data/Rhizopus_arrhizus.dna.fa,Mucoromycota;Mucoromycotina;Mucoromycetes;;Mucorales;Rhizopodaceae;Rhizopus,64495
"""


def test_load_master_pool_parses_fixed_width_lineage(tmp_path):
    p = tmp_path / 'pool.csv'
    p.write_text(POOL_CSV)
    samples = load_master_pool(p)
    assert len(samples) == 2
    mucor = samples[0]
    assert mucor.species == 'Mucor circinelloides'
    assert mucor.strain == '1006PhL'
    assert mucor.protein_path == '/data/Mucor_circinelloides.pep.fa'
    assert mucor.lineage == [
        'Mucoromycota', 'Mucoromycotina', 'Mucoromycetes', '', 'Mucorales', 'Mucoraceae', 'Mucor',
    ]
    assert mucor.ncbi_taxid == '36698'


def test_load_master_pool_rejects_wrong_lineage_width(tmp_path):
    p = tmp_path / 'bad_pool.csv'
    p.write_text(
        "Species,Strain,ProteinPath,DNAPath,Lineage,NCBI_TaxID\n"
        "Mucor circinelloides,,x,y,Mucoromycota;Mucorales;Mucor,36698\n"
    )
    with pytest.raises(SystemExit, match='expected 7'):
        load_master_pool(p)


def test_load_master_pool_rejects_duplicate_species(tmp_path):
    p = tmp_path / 'dup_pool.csv'
    p.write_text(POOL_CSV + POOL_CSV.splitlines()[1] + "\n")
    with pytest.raises(SystemExit, match='Duplicate Species'):
        load_master_pool(p)


def test_make_short_disambiguates_collisions():
    used: set[str] = set()
    # genus[:4] upper-first + lower + epithet[:4], concatenated and truncated to 8:
    # 'Muco' + 'circ' -> 'Mucocirc' (already 8 chars, no truncation)
    assert make_short('Mucor circinelloides', used) == 'Mucocirc'
    assert make_short('Mucor mucedo', used) == 'Mucomuce'
    # third species collides on the same 8-char base and gets a numeric suffix
    assert make_short('Mucor circinans', used) == 'Mucocir2'


def test_assign_shorts_is_deterministic_regardless_of_input_order(tmp_path):
    p = tmp_path / 'pool.csv'
    p.write_text(POOL_CSV)
    samples = load_master_pool(p)
    forward = assign_shorts(samples)
    backward = assign_shorts(list(reversed(samples)))
    assert forward == backward == {
        'Mucor circinelloides': 'Mucocirc',
        'Rhizopus arrhizus': 'Rhizarrh',
    }
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest tests/test_master_pool.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'master_pool'`

- [ ] **Step 3: Write the implementation**

```python
# lib/master_pool.py
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi run pytest tests/test_master_pool.py -v`
Expected: PASS (5 tests)

- [ ] **Step 5: Commit**

```bash
git add lib/master_pool.py tests/test_master_pool.py
git commit -m "feat: add master pool data structure, loader, and deterministic Short assignment"
```

---

### Task 3: Representative-assembly join

**Files:**
- Modify: `lib/master_pool.py` (add one function)
- Modify: `tests/test_master_pool.py` (add tests)

**Interfaces:**
- Produces: `load_representative_picks(path) -> dict[str, str]` (Species → the representative assembly's `out` dirname, from `repr_assignments.tsv`).

- [ ] **Step 1: Write the failing tests**

```python
# append to tests/test_master_pool.py
from master_pool import load_representative_picks

REPR_TSV = (
    "out\tspecies\tis_representative\trepresentative_out\tani_to_representative\treuse_eligible\n"
    "Mucor_circinelloides\tMucor circinelloides\tFalse\tMucor_circinelloides_1006PhL\t99.46\tTrue\n"
    "Mucor_circinelloides_1006PhL\tMucor circinelloides\tTrue\tMucor_circinelloides_1006PhL\t100.0\tTrue\n"
    "Rhizopus_arrhizus\tRhizopus arrhizus\tTrue\tRhizopus_arrhizus\t100.0\tFalse\n"
)

REPR_TSV_ZERO_TRUE = (
    "out\tspecies\tis_representative\trepresentative_out\tani_to_representative\treuse_eligible\n"
    "Phycomyces_blakesleeanus\tPhycomyces blakesleeanus\tFalse\tPhycomyces_blakesleeanus_NRRL1555\t98.0\tTrue\n"
)

REPR_TSV_MULTIPLE_TRUE = (
    "out\tspecies\tis_representative\trepresentative_out\tani_to_representative\treuse_eligible\n"
    "Basidiobolus_A\tBasidiobolus meristosporus\tTrue\tBasidiobolus_A\t100.0\tFalse\n"
    "Basidiobolus_B\tBasidiobolus meristosporus\tTrue\tBasidiobolus_A\t99.9\tFalse\n"
)


def test_load_representative_picks(tmp_path):
    p = tmp_path / 'repr.tsv'
    p.write_text(REPR_TSV)
    picks = load_representative_picks(p)
    assert picks == {
        'Mucor circinelloides': 'Mucor_circinelloides_1006PhL',
        'Rhizopus arrhizus': 'Rhizopus_arrhizus',
    }


def test_load_representative_picks_rejects_zero_true_rows(tmp_path):
    p = tmp_path / 'repr_zero.tsv'
    p.write_text(REPR_TSV_ZERO_TRUE)
    with pytest.raises(SystemExit, match='Phycomyces blakesleeanus'):
        load_representative_picks(p)


def test_load_representative_picks_rejects_multiple_true_rows(tmp_path):
    p = tmp_path / 'repr_multi.tsv'
    p.write_text(REPR_TSV_MULTIPLE_TRUE)
    with pytest.raises(SystemExit, match='Basidiobolus meristosporus'):
        load_representative_picks(p)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest tests/test_master_pool.py -v -k representative`
Expected: FAIL with `ImportError: cannot import name 'load_representative_picks'`

- [ ] **Step 3: Write the implementation**

Append to `lib/master_pool.py`:

```python
def load_representative_picks(path, expected_species: set[str] | None = None) -> dict[str, str]:
    """Species -> representative assembly's `out` dirname, from
    repr_assignments.tsv (columns: out, species, is_representative, ...).
    Hard-errors if a species has zero or more than one is_representative
    == 'True' row. If expected_species is given, also hard-errors on any
    species in it missing from the result entirely (zero rows at all,
    not just zero True rows).
    """
    import csv as _csv

    true_rows: dict[str, list[str]] = {}
    seen_species: set[str] = set()
    with open(path, newline='') as fh:
        for row in _csv.DictReader(fh, delimiter='\t'):
            seen_species.add(row['species'])
            if row['is_representative'] == 'True':
                true_rows.setdefault(row['species'], []).append(row['out'])

    if expected_species is not None:
        never_seen = expected_species - seen_species
        for sp in never_seen:
            true_rows.setdefault(sp, [])

    bad = {sp: outs for sp, outs in true_rows.items() if len(outs) != 1}
    if bad:
        lines = [
            f"{sp}: {len(outs)} is_representative=True row(s) ({', '.join(outs) or 'none'})"
            for sp, outs in sorted(bad.items())
        ]
        raise SystemExit(f"{path}: species without exactly one representative row:\n  " + "\n  ".join(lines))

    return {sp: outs[0] for sp, outs in true_rows.items()}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi run pytest tests/test_master_pool.py -v`
Expected: PASS (8 tests total)

- [ ] **Step 5: Commit**

```bash
git add lib/master_pool.py tests/test_master_pool.py
git commit -m "feat: join repr_assignments.tsv representative-assembly picks, hard-error on 0/>1 True rows"
```

---

### Task 4: `bin/build_master_pool.py` — the pool-building CLI

**Files:**
- Create: `bin/build_master_pool.py`
- Test: `tests/test_build_master_pool.py`

**Interfaces:**
- Consumes: `RANK_NAMES` (lineage), `load_representative_picks` (master_pool, Task 3), `MASTER_POOL_FIELDS` (master_pool, Task 2).
- Produces: `render_master_pool(bfd_samples_path, annotation_dir, repr_assignments_path) -> list[dict]` (rows in `MASTER_POOL_FIELDS` shape, ready to `csv.DictWriter`), importable from the CLI's `main()` and directly from tests.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_build_master_pool.py
import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))
from build_master_pool import render_master_pool  # noqa: E402

BFD_SAMPLES_CSV = (
    "ASMID,SPECIES_IN,STRAIN,BIOPROJECT,NCBI_TAXONID,BUSCO_LINEAGE,PHYLUM,SUBPHYLUM,CLASS,SUBCLASS,"
    "ORDER,FAMILY,GENUS,SPECIES,TRANSL_TABLE,LOCUSTAG\n"
    "GCA_1,Mucor circinelloides,1006PhL,PRJ1,36698,fungi_odb12,Mucoromycota,Mucoromycotina,"
    "Mucoromycetes,,Mucorales,Mucoraceae,Mucor,Mucor circinelloides,1,Mucci\n"
)

REPR_TSV = (
    "out\tspecies\tis_representative\trepresentative_out\tani_to_representative\treuse_eligible\n"
    "Mucor_circinelloides_1006PhL\tMucor circinelloides\tTrue\tMucor_circinelloides_1006PhL\t100.0\tTrue\n"
)


def _make_annotation_dir(tmp_path):
    d = tmp_path / 'genome_annotation' / 'Mucor_circinelloides_1006PhL' / 'predict_results'
    d.mkdir(parents=True)
    (d / 'Mucor_circinelloides_1006PhL.proteins.fa').write_text('>x\nMKV\n')
    (d / 'Mucor_circinelloides_1006PhL.scaffolds.fa').write_text('>x\nACGT\n')
    return tmp_path / 'genome_annotation'


def test_render_master_pool_joins_representative_pick_and_keeps_absolute_paths(tmp_path):
    bfd = tmp_path / 'bfd_samples.csv'
    bfd.write_text(BFD_SAMPLES_CSV)
    repr_tsv = tmp_path / 'repr_assignments.tsv'
    repr_tsv.write_text(REPR_TSV)
    annotation_dir = _make_annotation_dir(tmp_path)

    rows = render_master_pool(bfd, annotation_dir, repr_tsv)

    assert len(rows) == 1
    row = rows[0]
    assert row['Species'] == 'Mucor circinelloides'
    assert row['Strain'] == '1006PhL'
    assert row['NCBI_TaxID'] == '36698'
    assert row['Lineage'] == 'Mucoromycota;Mucoromycotina;Mucoromycetes;;Mucorales;Mucoraceae;Mucor'
    assert row['ProteinPath'] == str(annotation_dir / 'Mucor_circinelloides_1006PhL' / 'predict_results' / 'Mucor_circinelloides_1006PhL.proteins.fa')
    assert Path(row['ProteinPath']).is_absolute()


def test_render_master_pool_errors_on_missing_annotation_dir(tmp_path):
    bfd = tmp_path / 'bfd_samples.csv'
    bfd.write_text(BFD_SAMPLES_CSV)
    repr_tsv = tmp_path / 'repr_assignments.tsv'
    repr_tsv.write_text(REPR_TSV)
    empty_annotation_dir = tmp_path / 'no_such_dir'

    with pytest.raises(SystemExit, match='Mucor circinelloides'):
        render_master_pool(bfd, empty_annotation_dir, repr_tsv)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run pytest tests/test_build_master_pool.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'build_master_pool'`

- [ ] **Step 3: Write the implementation**

```python
#!/usr/bin/env python3
"""
Build the master species pool consumed by bin/build_targeted_configs.py,
from a BFD samples.csv (ASMID,SPECIES_IN,STRAIN,BIOPROJECT,NCBI_TAXONID,
BUSCO_LINEAGE,PHYLUM,SUBPHYLUM,CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIES,
TRANSL_TABLE,LOCUSTAG) plus its matching repr_assignments.tsv
(_reuse_assignments/repr_assignments.tsv: out,species,is_representative,
representative_out,ani_to_representative,reuse_eligible).

One row per species (never per assembly): the assembly whose dirname
matches that species' repr_assignments.tsv is_representative=True `out`
value is the only one kept. Unlike bin/convert_bfd_samples.py, Lineage
keeps a fixed 7-field width (PHYLUM;SUBPHYLUM;CLASS;SUBCLASS;ORDER;
FAMILY;GENUS, '' for an unrecorded rank, never dropped) and
Protein/DNA are resolved to absolute paths, not basenames -- see
docs/superpowers/specs/2026-09-05-config-builder-design.md ("Master pool").

Usage:
    bin/build_master_pool.py \\
        --bfd-samples /path/to/Fungi_BFD/samples.csv \\
        --annotation-dir /path/to/Fungi_BFD/genome_annotation \\
        --repr-assignments /path/to/Fungi_BFD_runs/genome_annotation/_reuse_assignments/repr_assignments.tsv \\
        --output config_support/master_pool.csv
"""
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from lineage import RANK_NAMES  # noqa: E402
from master_pool import MASTER_POOL_FIELDS, load_representative_picks  # noqa: E402

LINEAGE_FIELDS = ['PHYLUM', 'SUBPHYLUM', 'CLASS', 'SUBCLASS', 'ORDER', 'FAMILY', 'GENUS']
assert LINEAGE_FIELDS == RANK_NAMES  # single source of truth check


def find_annotation(annotation_dir, dirname):
    d = Path(annotation_dir) / dirname / 'predict_results'
    protein = d / f"{dirname}.proteins.fa"
    scaffolds = d / f"{dirname}.scaffolds.fa"
    if protein.exists() and scaffolds.exists():
        return protein.resolve(), scaffolds.resolve()
    return None, None


def render_master_pool(bfd_samples_path, annotation_dir, repr_assignments_path) -> list[dict]:
    with open(bfd_samples_path, newline='') as fh:
        bfd_rows = list(csv.DictReader(fh))

    species_set = {row['SPECIES'].strip() for row in bfd_rows}
    picks = load_representative_picks(repr_assignments_path, expected_species=species_set)

    by_dirname = {
        f"{row['SPECIES'].strip()}_{row['STRAIN'].strip()}".rstrip('_').replace(' ', '_'): row
        for row in bfd_rows
    }

    out_rows = []
    missing = []
    for species, dirname in sorted(picks.items()):
        row = by_dirname.get(dirname)
        if row is None:
            missing.append(f"{species}: representative dirname {dirname!r} not found in {bfd_samples_path}")
            continue
        protein, scaffolds = find_annotation(annotation_dir, dirname)
        if protein is None:
            missing.append(f"{species}: no predict_results/{dirname}.proteins.fa+scaffolds.fa under {annotation_dir}")
            continue
        lineage = [row[f].strip() for f in LINEAGE_FIELDS]
        out_rows.append({
            'Species': species,
            'Strain': row['STRAIN'].strip(),
            'ProteinPath': str(protein),
            'DNAPath': str(scaffolds),
            'Lineage': ';'.join(lineage),
            'NCBI_TaxID': row['NCBI_TAXONID'].strip(),
        })

    if missing:
        raise SystemExit("Could not build master pool rows:\n  " + "\n  ".join(missing))
    return out_rows


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--bfd-samples', required=True, help='Fungi_BFD/samples.csv')
    p.add_argument('--annotation-dir', required=True, help='Fungi_BFD/genome_annotation directory')
    p.add_argument('--repr-assignments', required=True, help='_reuse_assignments/repr_assignments.tsv')
    p.add_argument('--output', required=True, help='Destination master pool CSV')
    args = p.parse_args()

    rows = render_master_pool(args.bfd_samples, args.annotation_dir, args.repr_assignments)

    with open(args.output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=MASTER_POOL_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {args.output}: {len(rows)} species", file=sys.stderr)


if __name__ == '__main__':
    main()
```

- [ ] **Step 4: Make it executable and run tests to verify they pass**

Run: `chmod +x bin/build_master_pool.py && pixi run pytest tests/test_build_master_pool.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add bin/build_master_pool.py tests/test_build_master_pool.py
git commit -m "feat: add bin/build_master_pool.py, joining BFD samples with representative-assembly picks"
```

---

### Task 5: Trait loader

**Files:**
- Create: `lib/trait_data.py`
- Test: `tests/test_trait_data.py`

**Interfaces:**
- Produces: `TraitRow` dataclass (`species`, `trait`, `value`, `source`, `notes`), `load_trait_definitions(path) -> dict[str, set[str]]` (trait name → legal values), `load_traits(path, definitions) -> dict[str, list[TraitRow]]` (Species → its rows), `has_trait(traits_by_species, species, trait, value) -> bool`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_trait_data.py
import pytest
from trait_data import has_trait, load_trait_definitions, load_traits

DEFS_YAML = """\
traits:
  spore_motility:
    description: motility
    values:
      motile:
        description: flagellated
        ontology_term:
      nonmotile:
        description: not flagellated
        ontology_term:
  animal_association:
    description: host association
    values:
      none:
        description: none
        ontology_term:
      pathogen:
        description: pathogen
        ontology_term:
"""

TRAITS_CSV = """\
Species,trait,value,source,notes
Mucor circinelloides,spore_motility,nonmotile,,sporangiospores
Mucor circinelloides,animal_association,pathogen,,opportunistic
"""


def test_load_trait_definitions():
    defs = load_trait_definitions_from_text(DEFS_YAML)
    assert defs == {
        'spore_motility': {'motile', 'nonmotile'},
        'animal_association': {'none', 'pathogen'},
    }


def load_trait_definitions_from_text(text, tmp_path=None):
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / 'defs.yaml'
        p.write_text(text)
        return load_trait_definitions(p)


def test_load_traits_and_has_trait(tmp_path):
    defs_path = tmp_path / 'defs.yaml'
    defs_path.write_text(DEFS_YAML)
    traits_path = tmp_path / 'traits.csv'
    traits_path.write_text(TRAITS_CSV)

    defs = load_trait_definitions(defs_path)
    by_species = load_traits(traits_path, defs)

    assert has_trait(by_species, 'Mucor circinelloides', 'animal_association', 'pathogen')
    assert not has_trait(by_species, 'Mucor circinelloides', 'animal_association', 'none')
    assert not has_trait(by_species, 'Phycomyces blakesleeanus', 'spore_motility', 'motile')


def test_load_traits_rejects_undeclared_trait(tmp_path):
    defs_path = tmp_path / 'defs.yaml'
    defs_path.write_text(DEFS_YAML)
    traits_path = tmp_path / 'traits.csv'
    traits_path.write_text("Species,trait,value,source,notes\nMucor circinelloides,made_up_trait,x,,\n")

    with pytest.raises(SystemExit, match='undeclared trait'):
        load_traits(traits_path, load_trait_definitions(defs_path))


def test_load_traits_rejects_undeclared_value(tmp_path):
    defs_path = tmp_path / 'defs.yaml'
    defs_path.write_text(DEFS_YAML)
    traits_path = tmp_path / 'traits.csv'
    traits_path.write_text("Species,trait,value,source,notes\nMucor circinelloides,spore_motility,flying,,\n")

    with pytest.raises(SystemExit, match='undeclared value'):
        load_traits(traits_path, load_trait_definitions(defs_path))


def test_load_traits_rejects_none_coexisting_with_another_value(tmp_path):
    defs_path = tmp_path / 'defs.yaml'
    defs_path.write_text(DEFS_YAML)
    traits_path = tmp_path / 'traits.csv'
    traits_path.write_text(
        "Species,trait,value,source,notes\n"
        "Mucor circinelloides,animal_association,none,,\n"
        "Mucor circinelloides,animal_association,pathogen,,\n"
    )

    with pytest.raises(SystemExit, match="'none' must never coexist"):
        load_traits(traits_path, load_trait_definitions(defs_path))
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest tests/test_trait_data.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'trait_data'`

- [ ] **Step 3: Write the implementation**

```python
# lib/trait_data.py
"""Species trait loader: config_support/traits/trait_definitions.yaml (the
controlled vocabulary) + config_support/traits/traits.csv (the actual
species-to-trait rows). See config_support/traits/README.md and
docs/superpowers/specs/2026-09-05-config-builder-design.md ("Trait data").
"""
import csv
from dataclasses import dataclass

import yaml


@dataclass
class TraitRow:
    species: str
    trait: str
    value: str
    source: str
    notes: str


def load_trait_definitions(path) -> dict[str, set[str]]:
    with open(path) as fh:
        doc = yaml.safe_load(fh)
    return {trait: set(spec['values'].keys()) for trait, spec in doc['traits'].items()}


def load_traits(path, definitions: dict[str, set[str]]) -> dict[str, list[TraitRow]]:
    by_species: dict[str, list[TraitRow]] = {}
    with open(path, newline='') as fh:
        for row in csv.DictReader(fh):
            trait = row['trait'].strip()
            value = row['value'].strip()
            species = row['Species'].strip()
            if trait not in definitions:
                raise SystemExit(f"{path}: undeclared trait {trait!r} (species: {species!r})")
            if value not in definitions[trait]:
                raise SystemExit(
                    f"{path}: undeclared value {value!r} for trait {trait!r} (species: {species!r})"
                )
            by_species.setdefault(species, []).append(TraitRow(
                species=species, trait=trait, value=value,
                source=(row.get('source') or '').strip(),
                notes=(row.get('notes') or '').strip(),
            ))

    for species, rows in by_species.items():
        values_by_trait: dict[str, set[str]] = {}
        for r in rows:
            values_by_trait.setdefault(r.trait, set()).add(r.value)
        for trait, values in values_by_trait.items():
            if 'none' in values and len(values) > 1:
                raise SystemExit(
                    f"{path}: {species!r} trait {trait!r} has 'none' coexisting with "
                    f"{sorted(values - {'none'})} -- 'none' must never coexist with another value"
                )

    return by_species


def has_trait(traits_by_species: dict[str, list[TraitRow]], species: str, trait: str, value: str) -> bool:
    return any(r.trait == trait and r.value == value for r in traits_by_species.get(species, []))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi run pytest tests/test_trait_data.py -v`
Expected: PASS (5 tests)

- [ ] **Step 5: Commit**

```bash
git add lib/trait_data.py tests/test_trait_data.py
git commit -m "feat: add trait_definitions.yaml/traits.csv loader with hard-error validation"
```

---

### Task 6: Candidate selection (`nearest` / `trait` / disjointness)

**Files:**
- Create: `lib/targeted_selection.py`
- Test: `tests/test_targeted_selection.py`

**Interfaces:**
- Consumes: `RANK_NAMES`, `lineage_match` (lineage, Task 1); `MasterSample` (master_pool, Task 2); `has_trait` (trait_data, Task 5).
- Produces: `Candidate` dataclass (`species: str`, `rank_name: str`), `DEFAULT_SCOPE_RANK = 'ORDER'`, `candidate_pool(focal_species, focal_lineage, pool, scope_rank=DEFAULT_SCOPE_RANK) -> list[Candidate]`, `rank_candidates(candidates) -> list[Candidate]`, `exclude_species(candidates, excluded: set[str]) -> list[Candidate]`, `select_nearest(focal_species, focal_lineage, pool, n, scope_rank=DEFAULT_SCOPE_RANK, excluded: set[str] = frozenset()) -> list[Candidate]`, `select_trait(focal_species, focal_lineage, pool, trait, value, n, traits_by_species, scope_rank=DEFAULT_SCOPE_RANK, excluded: set[str] = frozenset()) -> list[Candidate]`.
  **Note:** `excluded` is applied *before* ranking/truncation inside `select_nearest`/`select_trait`, not after — this is deliberate: a caller that truncates to `n` and only then excludes outgroup-pool members would silently under-fill a study whenever a top candidate happens to be excluded, instead of backfilling from the next-best candidate. See the new `test_select_nearest_backfills_past_an_excluded_top_candidate` test below, which exists specifically to catch a regression to the wrong order.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_targeted_selection.py
import pytest
from master_pool import MasterSample
from targeted_selection import (
    candidate_pool,
    exclude_species,
    rank_candidates,
    select_nearest,
    select_trait,
)


def _sample(species, phylum, klass, order, family, genus):
    return MasterSample(
        species=species, strain='', protein_path='', dna_path='',
        lineage=[phylum, '', klass, '', order, family, genus], ncbi_taxid='',
    )


FOCAL = _sample('Mucor circinelloides', 'Mucoromycota', 'Mucoromycetes', 'Mucorales', 'Mucoraceae', 'Mucor')
POOL = [
    FOCAL,
    _sample('Rhizopus arrhizus', 'Mucoromycota', 'Mucoromycetes', 'Mucorales', 'Rhizopodaceae', 'Rhizopus'),
    _sample('Lichtheimia corymbifera', 'Mucoromycota', 'Mucoromycetes', 'Mucorales', 'Lichtheimiaceae', 'Lichtheimia'),
    _sample('Basidiobolus meristosporus', 'Basidiobolomycota', 'Basidiobolomycetes', 'Basidiobolales', 'Basidiobolaceae', 'Basidiobolus'),
]


def test_candidate_pool_excludes_focal_and_out_of_scope_species():
    candidates = candidate_pool(FOCAL.species, FOCAL.lineage, POOL, scope_rank='ORDER')
    assert {c.species for c in candidates} == {'Rhizopus arrhizus', 'Lichtheimia corymbifera'}


def test_rank_candidates_orders_deepest_match_first_then_alphabetically():
    same_family = _sample('Mucor mucedo', 'Mucoromycota', 'Mucoromycetes', 'Mucorales', 'Mucoraceae', 'Mucor')
    candidates = candidate_pool(FOCAL.species, FOCAL.lineage, POOL + [same_family], scope_rank='ORDER')
    ranked = rank_candidates(candidates)
    assert [c.species for c in ranked] == ['Mucor mucedo', 'Lichtheimia corymbifera', 'Rhizopus arrhizus']


def test_select_nearest_takes_top_n():
    picked = select_nearest(FOCAL.species, FOCAL.lineage, POOL, n=1, scope_rank='ORDER')
    assert len(picked) == 1
    # deterministic tiebreak: both tie at rank ORDER, 'Lichtheimia' sorts before 'Rhizopus'
    assert picked[0].species == 'Lichtheimia corymbifera'


def test_select_nearest_backfills_past_an_excluded_top_candidate():
    # Mucor mucedo shares GENUS with the focal (closer than Rhizopus/Lichtheimia,
    # which only share ORDER) -- excluding it must NOT just shrink the result to
    # fewer than n; it must promote the next-best candidate instead.
    same_genus = _sample('Mucor mucedo', 'Mucoromycota', 'Mucoromycetes', 'Mucorales', 'Mucoraceae', 'Mucor')
    picked = select_nearest(
        FOCAL.species, FOCAL.lineage, POOL + [same_genus], n=1,
        scope_rank='ORDER', excluded={'Mucor mucedo'},
    )
    assert [c.species for c in picked] == ['Lichtheimia corymbifera']


def test_select_trait_filters_then_ranks():
    traits_by_species = {
        'Rhizopus arrhizus': [_TraitRow('Rhizopus arrhizus', 'thermotolerance', 'high')],
    }
    picked = select_trait(
        FOCAL.species, FOCAL.lineage, POOL, 'thermotolerance', 'high', n=2,
        traits_by_species=traits_by_species, scope_rank='ORDER',
    )
    assert [c.species for c in picked] == ['Rhizopus arrhizus']


def test_select_trait_errors_on_empty_filtered_pool():
    with pytest.raises(SystemExit, match='no candidate'):
        select_trait(
            FOCAL.species, FOCAL.lineage, POOL, 'thermotolerance', 'high', n=2,
            traits_by_species={}, scope_rank='ORDER',
        )


def test_exclude_species_removes_outgroup_pool_members():
    candidates = rank_candidates(candidate_pool(FOCAL.species, FOCAL.lineage, POOL, scope_rank='ORDER'))
    remaining = exclude_species(candidates, {'Rhizopus arrhizus'})
    assert [c.species for c in remaining] == ['Lichtheimia corymbifera']


def _TraitRow(species, trait, value):
    from trait_data import TraitRow
    return TraitRow(species=species, trait=trait, value=value, source='', notes='')
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest tests/test_targeted_selection.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'targeted_selection'`

- [ ] **Step 3: Write the implementation**

```python
# lib/targeted_selection.py
"""mode: nearest / trait / explicit candidate selection, and the
ingroup/outgroup disjointness rule. See
docs/superpowers/specs/2026-09-05-config-builder-design.md ("Trait mode",
"Ingroup/outgroup disjointness").
"""
from dataclasses import dataclass

from lineage import RANK_NAMES, lineage_match
from trait_data import has_trait

DEFAULT_SCOPE_RANK = 'ORDER'


@dataclass
class Candidate:
    species: str
    rank_name: str  # deepest shared rank name with the focal, e.g. 'FAMILY'


def _depth(rank_name: str) -> int:
    """Sort key: deeper rank = smaller index = closer. '' (no match) sorts last."""
    return RANK_NAMES.index(rank_name) if rank_name else -1


def candidate_pool(focal_species: str, focal_lineage: list[str], pool: list, scope_rank: str = DEFAULT_SCOPE_RANK) -> list[Candidate]:
    """Every species in `pool` (excluding the focal itself) sharing at
    least rank `scope_rank` with the focal."""
    if scope_rank not in RANK_NAMES:
        raise ValueError(f"scope_rank {scope_rank!r} must be one of {RANK_NAMES}")
    scope_depth = RANK_NAMES.index(scope_rank)
    out = []
    for s in pool:
        if s.species == focal_species:
            continue
        rank_name = lineage_match(focal_lineage, s.lineage)
        if _depth(rank_name) >= scope_depth:
            out.append(Candidate(species=s.species, rank_name=rank_name))
    return out


def rank_candidates(candidates: list[Candidate]) -> list[Candidate]:
    """Deepest match first; ties broken alphabetically by species name."""
    return sorted(candidates, key=lambda c: (-_depth(c.rank_name), c.species))


def exclude_species(candidates: list[Candidate], excluded: set[str]) -> list[Candidate]:
    return [c for c in candidates if c.species not in excluded]


def select_nearest(focal_species, focal_lineage, pool, n, scope_rank: str = DEFAULT_SCOPE_RANK, excluded: set[str] = frozenset()) -> list[Candidate]:
    # excluded is applied BEFORE ranking/truncation -- see this task's Interfaces
    # note: excluding after truncation would silently under-fill a study instead
    # of backfilling from the next-best candidate.
    candidates = exclude_species(candidate_pool(focal_species, focal_lineage, pool, scope_rank), excluded)
    return rank_candidates(candidates)[:n]


def select_trait(focal_species, focal_lineage, pool, trait, value, n, traits_by_species, scope_rank: str = DEFAULT_SCOPE_RANK, excluded: set[str] = frozenset()) -> list[Candidate]:
    candidates = exclude_species(candidate_pool(focal_species, focal_lineage, pool, scope_rank), excluded)
    ranked = rank_candidates(candidates)
    filtered = [c for c in ranked if has_trait(traits_by_species, c.species, trait, value)]
    if not filtered:
        raise SystemExit(
            f"mode: trait -- no candidate for focal {focal_species!r} in the "
            f"{scope_rank}-scoped pool has {trait}={value!r} (after excluding the outgroup pool)"
        )
    return filtered[:n]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi run pytest tests/test_targeted_selection.py -v`
Expected: PASS (6 tests)

- [ ] **Step 5: Commit**

```bash
git add lib/targeted_selection.py tests/test_targeted_selection.py
git commit -m "feat: add nearest/trait candidate selection with lineage-scoped filtering and disjointness helper"
```

---

### Task 7: `bin/build_targeted_configs.py` — the renderer

**Files:**
- Create: `bin/build_targeted_configs.py`
- Test: `tests/test_build_targeted_configs.py`

**Interfaces:**
- Consumes: everything from Tasks 1-6.
- Produces: `render_batch(master_pool_path, trait_definitions_path, traits_path, batch_spec_path, outdir) -> list[dict]` (one summary dict per study: `{'batch': str, 'focal': str, 'config_path': str, 'map_path': str, 'companions': list[dict], 'outgroup_pool': str, 'outgroup_size': int}`), and `main()` wiring it to argparse.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_build_targeted_configs.py
import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))
from build_targeted_configs import render_batch  # noqa: E402

MASTER_POOL_CSV = """\
Species,Strain,ProteinPath,DNAPath,Lineage,NCBI_TaxID
Mucor circinelloides,1006PhL,/data/Mucci.pep.fa,/data/Mucci.dna.fa,Mucoromycota;Mucoromycotina;Mucoromycetes;;Mucorales;Mucoraceae;Mucor,36698
Rhizopus arrhizus,,/data/Rhiar.pep.fa,/data/Rhiar.dna.fa,Mucoromycota;Mucoromycotina;Mucoromycetes;;Mucorales;Rhizopodaceae;Rhizopus,64495
Lichtheimia corymbifera,,/data/Licor.pep.fa,/data/Licor.dna.fa,Mucoromycota;Mucoromycotina;Mucoromycetes;;Mucorales;Lichtheimiaceae;Lichtheimia,64752
Basidiobolus meristosporus,,/data/Bamer.pep.fa,/data/Bamer.dna.fa,Basidiobolomycota;;Basidiobolomycetes;;Basidiobolales;Basidiobolaceae;Basidiobolus,,
Neurospora crassa,OR74A,/data/Ncra.pep.fa,/data/Ncra.dna.fa,Ascomycota;Pezizomycotina;Sordariomycetes;;Sordariales;Sordariaceae;Neurospora,367110
Aspergillus nidulans,,/data/Anid.pep.fa,/data/Anid.dna.fa,Ascomycota;Pezizomycotina;Eurotiomycetes;;Eurotiales;Aspergillaceae;Aspergillus,227321
"""

TRAIT_DEFINITIONS_YAML = """\
traits:
  thermotolerance:
    description: x
    values:
      high: {description: x, ontology_term:}
      low: {description: x, ontology_term:}
"""

TRAITS_CSV = """\
Species,trait,value,source,notes
Rhizopus arrhizus,thermotolerance,high,,
"""

BATCH_YAML = """\
batch: test_batch_v1
outgroup_pools:
  dikarya_v1:
    members: ["Neurospora crassa", "Aspergillus nidulans"]
studies:
  - focal: Mucor circinelloides
    ingroup_extra: {mode: nearest, n: 2}
    outgroup_pool: dikarya_v1
"""

BATCH_YAML_TRAIT_MODE = """\
batch: test_batch_v2
outgroup_pools:
  dikarya_v1:
    members: ["Neurospora crassa", "Aspergillus nidulans"]
studies:
  - focal: Mucor circinelloides
    ingroup_extra: {mode: trait, trait: thermotolerance, value: high, n: 1}
    outgroup_pool: dikarya_v1
"""


def _write_fixtures(tmp_path, batch_yaml=BATCH_YAML):
    pool = tmp_path / 'pool.csv'
    pool.write_text(MASTER_POOL_CSV)
    defs = tmp_path / 'trait_definitions.yaml'
    defs.write_text(TRAIT_DEFINITIONS_YAML)
    traits = tmp_path / 'traits.csv'
    traits.write_text(TRAITS_CSV)
    batch = tmp_path / 'batch.yaml'
    batch.write_text(batch_yaml)
    outdir = tmp_path / 'out'
    outdir.mkdir()
    return pool, defs, traits, batch, outdir


def test_render_batch_nearest_mode_produces_config_and_map(tmp_path):
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path)
    summaries = render_batch(pool, defs, traits, batch, outdir)

    assert len(summaries) == 1
    summary = summaries[0]
    assert summary['focal'] == 'Mucor circinelloides'
    assert {c['species'] for c in summary['companions']} == {'Rhizopus arrhizus', 'Lichtheimia corymbifera'}

    config_path = Path(summary['config_path'])
    assert config_path.exists()
    with open(config_path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    groups = {r['Species']: r['GROUP'] for r in rows}
    assert groups == {
        'Mucor circinelloides': 'IN',
        'Rhizopus arrhizus': 'IN',
        'Lichtheimia corymbifera': 'IN',
        'Neurospora crassa': 'OUT',
        'Aspergillus nidulans': 'OUT',
    }

    map_path = Path(summary['map_path'])
    assert map_path.exists()
    map_text = map_path.read_text()
    assert 'Mucor circinelloides\t' in map_text


def test_render_batch_trait_mode_filters_to_matching_candidate(tmp_path):
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path, batch_yaml=BATCH_YAML_TRAIT_MODE)
    summaries = render_batch(pool, defs, traits, batch, outdir)
    assert [c['species'] for c in summaries[0]['companions']] == ['Rhizopus arrhizus']


def test_render_batch_errors_when_focal_is_inside_its_own_outgroup_pool(tmp_path):
    bad_batch = """\
batch: bad
outgroup_pools:
  bad_pool:
    members: ["Mucor circinelloides"]
studies:
  - focal: Mucor circinelloides
    ingroup_extra: {mode: nearest, n: 1}
    outgroup_pool: bad_pool
"""
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path, batch_yaml=bad_batch)
    with pytest.raises(SystemExit, match='inside its own outgroup pool'):
        render_batch(pool, defs, traits, batch, outdir)


def test_render_batch_excludes_outgroup_members_from_ingroup_candidates(tmp_path):
    # Neurospora crassa can never be picked as an ORDER-scoped companion of
    # Mucor here anyway (different phylum), so this test uses a pool member
    # that WOULD otherwise rank -- Rhizopus -- placed into the outgroup pool.
    batch_yaml = """\
batch: disjoint
outgroup_pools:
  contains_rhizopus:
    members: ["Rhizopus arrhizus", "Neurospora crassa"]
studies:
  - focal: Mucor circinelloides
    ingroup_extra: {mode: nearest, n: 2}
    outgroup_pool: contains_rhizopus
"""
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path, batch_yaml=batch_yaml)
    summaries = render_batch(pool, defs, traits, batch, outdir)
    assert [c['species'] for c in summaries[0]['companions']] == ['Lichtheimia corymbifera']
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run pytest tests/test_build_targeted_configs.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'build_targeted_configs'`

- [ ] **Step 3: Write the implementation**

```python
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
        --outdir configs/
"""
import argparse
import csv
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from master_pool import assign_shorts, load_master_pool  # noqa: E402
from targeted_selection import (  # noqa: E402
    DEFAULT_SCOPE_RANK,
    select_nearest,
    select_trait,
)
from trait_data import load_trait_definitions, load_traits  # noqa: E402

CONFIG_FIELDS = ['GROUP', 'Species', 'Strain', 'Protein', 'DNA', 'Short', 'TaxonGroup']


def _by_species(pool):
    return {s.species: s for s in pool}


def _resolve_members(names, by_species, context):
    missing = [n for n in names if n not in by_species]
    if missing:
        raise SystemExit(f"{context}: species not found in master pool: {missing}")
    return names


def _row(group, sample, short, taxon_group):
    return {
        'GROUP': group, 'Species': sample.species, 'Strain': sample.strain,
        'Protein': sample.protein_path, 'DNA': sample.dna_path, 'Short': short,
        'TaxonGroup': taxon_group,
    }


def render_batch(master_pool_path, trait_definitions_path, traits_path, batch_spec_path, outdir) -> list[dict]:
    pool = load_master_pool(master_pool_path)
    by_species = _by_species(pool)
    short_map = assign_shorts(pool)

    definitions = load_trait_definitions(trait_definitions_path)
    traits_by_species = load_traits(traits_path, definitions)

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
            companions = [{'species': c.species, 'taxon_group': c.rank_name, 'reason': ''} for c in candidates]
        elif mode == 'trait':
            candidates = select_trait(
                focal.species, focal.lineage, pool, extra['trait'], extra['value'], n=n,
                traits_by_species=traits_by_species, scope_rank=scope_rank, excluded=outgroup_set,
            )
            companions = [{'species': c.species, 'taxon_group': c.rank_name, 'reason': ''} for c in candidates]
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

        rows = [_row('IN', focal, short_map[focal.species], '')]
        for c in companions:
            rows.append(_row('IN', by_species[c['species']], short_map[c['species']], c['taxon_group']))
        for m in outgroup_members:
            s = by_species[m]
            taxon_group = [t for t in s.lineage if t][-1] if any(s.lineage) else ''
            rows.append(_row('OUT', s, short_map[m], taxon_group))

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
    args = p.parse_args()

    summaries = render_batch(args.master_pool, args.trait_definitions, args.traits, args.batch_spec, args.outdir)

    for s in summaries:
        print(f"\n{s['focal']}  ({s['config_path']})", file=sys.stderr)
        for c in s['companions']:
            reason = f" -- {c['reason']}" if c['reason'] else ''
            print(f"  + {c['species']}  [{c['taxon_group']}]{reason}", file=sys.stderr)
        print(f"  outgroup: {s['outgroup_pool']} ({s['outgroup_size']} species)", file=sys.stderr)


if __name__ == '__main__':
    main()
```

- [ ] **Step 4: Make it executable and run tests to verify they pass**

Run: `chmod +x bin/build_targeted_configs.py && pixi run pytest tests/test_build_targeted_configs.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add bin/build_targeted_configs.py tests/test_build_targeted_configs.py
git commit -m "feat: add bin/build_targeted_configs.py, rendering a YAML batch spec into config CSV(s)"
```

---

### Task 8: Full test suite and lint pass

**Files:**
- None created; verification only.

- [ ] **Step 1: Run the full test suite**

Run: `pixi run pytest tests/ -v`
Expected: PASS — every test from Tasks 1-7 plus the existing suite (`test_config_parser.py`, `test_hits.py`, etc.) green.

- [ ] **Step 2: Run lint**

Run: `pixi run lint`
Expected: no new violations in `bin/build_master_pool.py`, `bin/build_targeted_configs.py`, `lib/lineage.py`, `lib/master_pool.py`, `lib/trait_data.py`, `lib/targeted_selection.py`.

- [ ] **Step 3: Fix any lint violations and re-run**

If Step 2 reports issues, fix them in the relevant file(s) and re-run `pixi run lint` and `pixi run pytest tests/ -v` until both are clean.

- [ ] **Step 4: Commit (only if Step 3 required changes)**

```bash
git add -A
git commit -m "chore: fix lint violations in targeted config-builder modules"
```

---

## Spec coverage check

- Master pool (Species-keyed, absolute paths, NCBI_TaxID, Short-at-render-time + `.map.tsv`) → Tasks 2, 4.
- `repr_assignments.tsv` join, hard-error on 0/>1 True rows → Task 3.
- Trait data (`trait_definitions.yaml`/`traits.csv`, hard-error rules) → Task 5.
- Lineage-proximity ranking (prefix-based, rank-name-aligned, missing-rank case) → Task 1.
- `mode: nearest`/`trait`/`explicit`, trait-mode's filter-then-rank, empty-pool error → Tasks 6, 7.
- Ingroup/outgroup disjointness, focal-in-own-outgroup-pool error → Tasks 6, 7 (tests in Task 7).
- Batch YAML format, `batch:` key, per-study summary output, `.map.tsv` → Task 7.
- Workload estimate: the spec's v1 decision is "no per-study estimator needed at this scale" — no task required.
- Not in this plan (explicitly deferred by the spec): the `textual` TUI, browser GUI, `bin/import_funguild_traits.py`, whole-pool (non-lineage-scoped) `mode: trait`, tree-backed ranking, timing-calibrated estimation.
