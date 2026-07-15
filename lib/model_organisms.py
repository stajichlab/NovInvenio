"""
Model organism gene annotation support.

Loads a YAML config that lists reference organisms.  For each organism,
resolves protein IDs from the analysis to gene names and product descriptions.

Three ID-transform strategies are supported:

  direct        protein_id is used as the gene lookup key unchanged
  strip         a regex is stripped from protein_id to get the gene key
                  e.g.  Afu1g00220-T-p1  →  Afu1g00220   (strip_pattern: "-T-p\\d+$")
  diamond_fasta protein_id → diamond hit → FungiDB protein_id → FASTA header gene field
                  e.g.  FC69C3D3_000001-T1  (local Ncra)
                        → NCU10129-t26_1-p1  (FungiDB protein, via diamond_hits TSV)
                        → NCU10129           (via fasta_gene_field: "gene" in header)

YAML format::

  model_organisms:
    - short: Ncra
      gene_names_csv: db/modelorgs/Neurospora_crassa_gene_names_FungiDB.csv
      id_transform: diamond_fasta
      diamond_hits: db/modelorgs/Ncra_vs_FungiDB_Ncra.diamond.tsv
      protein_fasta: db/modelorgs/FungiDB-68_NcrassaOR74A_AnnotatedProteins.fasta
      fasta_gene_field: gene        # parse "gene=..." from FASTA header

    - short: Afum
      gene_names_csv: db/modelorgs/Aspergillus_fumigatus_gene_names_FungiDB.csv
      id_transform: strip
      strip_pattern: "-T-p\\d+$"

Optional column overrides (defaults match FungiDB CSV export format):
      gene_id_col:   "Gene ID"
      gene_name_col: "Gene Name or Symbol"
      product_col:   "Product Description"
      csv_delimiter: ","              # default comma; set to "\\t" for TSV
"""

import csv
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml


# ── Data classes ─────────────────────────────────────────────────────────────

@dataclass
class ModelOrgConfig:
    short: str
    gene_names_csv: str
    id_transform: str = 'direct'          # direct | strip | diamond_fasta

    # strip strategy
    strip_pattern: Optional[str] = None

    # diamond_fasta strategy
    diamond_hits: Optional[str] = None
    protein_fasta: Optional[str] = None
    fasta_gene_field: str = 'gene'

    # gene_names_csv column names (FungiDB defaults)
    gene_id_col: str = 'Gene ID'
    gene_name_col: str = 'Gene Name or Symbol'
    product_col: str = 'Product Description'
    csv_delimiter: str = ','


# ── Loader ────────────────────────────────────────────────────────────────────

def load_model_org_configs(yaml_path: str,
                           launch_dir: Optional[str] = None) -> dict[str, ModelOrgConfig]:
    """
    Parse a model organisms YAML file.
    Returns {short: ModelOrgConfig}.

    Relative paths in the YAML are resolved against *launch_dir* when provided
    (the project root, i.e. the directory where ``nextflow run`` was called),
    falling back to the YAML file's own directory when running outside Nextflow.
    """
    yaml_path = Path(yaml_path).resolve()
    base_dir = Path(launch_dir).resolve() if launch_dir else yaml_path.parent

    with open(yaml_path) as fh:
        data = yaml.safe_load(fh)

    entries = data.get('model_organisms', [])
    configs: dict[str, ModelOrgConfig] = {}

    for entry in entries:
        short = entry['short']
        mo = ModelOrgConfig(
            short=short,
            gene_names_csv=_resolve(entry['gene_names_csv'], base_dir),
            id_transform=entry.get('id_transform', 'direct'),
            strip_pattern=entry.get('strip_pattern'),
            diamond_hits=_resolve(entry.get('diamond_hits'), base_dir),
            protein_fasta=_resolve(entry.get('protein_fasta'), base_dir),
            fasta_gene_field=entry.get('fasta_gene_field', 'gene'),
            gene_id_col=entry.get('gene_id_col', 'Gene ID'),
            gene_name_col=entry.get('gene_name_col', 'Gene Name or Symbol'),
            product_col=entry.get('product_col', 'Product Description'),
            csv_delimiter=entry.get('csv_delimiter', ','),
        )

        _validate(mo)
        configs[short] = mo

    return configs


def _resolve(path_str: Optional[str], base_dir: Path) -> Optional[str]:
    if path_str is None:
        return None
    p = Path(path_str)
    if not p.is_absolute():
        p = base_dir / p
    return str(p)


def _validate(mo: ModelOrgConfig) -> None:
    if mo.id_transform not in ('direct', 'strip', 'diamond_fasta'):
        sys.exit(
            f"[model_organisms] Unknown id_transform '{mo.id_transform}' for {mo.short}. "
            "Choose: direct | strip | diamond_fasta"
        )
    if mo.id_transform == 'strip' and not mo.strip_pattern:
        sys.exit(f"[model_organisms] strip_pattern required when id_transform=strip for {mo.short}")
    if mo.id_transform == 'diamond_fasta':
        if not mo.diamond_hits:
            sys.exit(f"[model_organisms] diamond_hits required when id_transform=diamond_fasta for {mo.short}")
        if not mo.protein_fasta:
            sys.exit(f"[model_organisms] protein_fasta required when id_transform=diamond_fasta for {mo.short}")


# ── Annotator ─────────────────────────────────────────────────────────────────

class ModelOrgAnnotator:
    """
    Loads all reference data for a set of model organisms and provides
    fast lookup of gene name and product description.

    Usage::

        annotator = ModelOrgAnnotator.from_yaml("configs/modelorgs.yaml")
        gene_name, product = annotator.annotate("FC69C3D3_000001-T1", "Ncra")
    """

    def __init__(self, configs: dict[str, ModelOrgConfig]) -> None:
        self._configs = configs
        # Per-short lookup tables, loaded lazily on first use
        self._gene_info: dict[str, dict[str, tuple[str, str]]] = {}
        self._diamond: dict[str, dict[str, str]] = {}
        self._fasta_gene: dict[str, dict[str, str]] = {}

    @classmethod
    def from_yaml(cls, yaml_path: str,
                  launch_dir: Optional[str] = None) -> 'ModelOrgAnnotator':
        return cls(load_model_org_configs(yaml_path, launch_dir=launch_dir))

    # -- public --

    def annotate(self, protein_id: str, source_short: str) -> tuple[str, str]:
        """
        Return (gene_name, product_description) for the given protein.
        Returns ('', '') if no annotation is found or the source is not
        in the configured model organisms.
        """
        mo = self._configs.get(source_short)
        if mo is None:
            return '', ''

        gene_key = self._resolve_gene_key(protein_id, mo)
        if not gene_key:
            return '', ''

        info = self._gene_info_for(mo)
        result = info.get(gene_key, ('', ''))
        return result

    def shorts(self) -> list[str]:
        return list(self._configs.keys())

    # -- private --

    def _resolve_gene_key(self, protein_id: str, mo: ModelOrgConfig) -> Optional[str]:
        if mo.id_transform == 'direct':
            return protein_id

        if mo.id_transform == 'strip':
            return re.sub(mo.strip_pattern, '', protein_id)

        if mo.id_transform == 'diamond_fasta':
            diamond = self._diamond_map_for(mo)
            ref_protein = diamond.get(protein_id)
            if not ref_protein:
                return None
            fasta_genes = self._fasta_gene_map_for(mo)
            return fasta_genes.get(ref_protein)

        return None

    def _gene_info_for(self, mo: ModelOrgConfig) -> dict[str, tuple[str, str]]:
        if mo.short not in self._gene_info:
            self._gene_info[mo.short] = _load_gene_names_csv(mo)
        return self._gene_info[mo.short]

    def _diamond_map_for(self, mo: ModelOrgConfig) -> dict[str, str]:
        if mo.short not in self._diamond:
            self._diamond[mo.short] = _load_diamond_hits(mo.diamond_hits)
        return self._diamond[mo.short]

    def _fasta_gene_map_for(self, mo: ModelOrgConfig) -> dict[str, str]:
        if mo.short not in self._fasta_gene:
            self._fasta_gene[mo.short] = _load_fasta_gene_map(
                mo.protein_fasta, mo.fasta_gene_field
            )
        return self._fasta_gene[mo.short]


# ── File parsers ──────────────────────────────────────────────────────────────

def _load_gene_names_csv(mo: ModelOrgConfig) -> dict[str, tuple[str, str]]:
    """
    Return {gene_id: (gene_name, product_description)}.
    Handles quoted fields (FungiDB export wraps values in extra quotes).
    """
    info: dict[str, tuple[str, str]] = {}
    path = mo.gene_names_csv
    if not path or not Path(path).exists():
        print(
            f"[model_organisms] WARNING: gene_names_csv not found for {mo.short}: {path}",
            file=sys.stderr,
        )
        return info

    delim = mo.csv_delimiter.replace('\\t', '\t')

    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        for row in reader:
            gid = row.get(mo.gene_id_col, '').strip('"').strip()
            name = row.get(mo.gene_name_col, '').strip('"').strip()
            prod = row.get(mo.product_col, '').strip('"').strip()
            if gid:
                info[gid] = (
                    '' if name == 'N/A' else name,
                    '' if prod == 'N/A' else prod,
                )
    return info


def _load_diamond_hits(tsv_path: Optional[str]) -> dict[str, str]:
    """Return {query_id: subject_id} keeping only the first (best) hit per query."""
    hits: dict[str, str] = {}
    if not tsv_path or not Path(tsv_path).exists():
        return hits
    with open(tsv_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            qid, sid = parts[0], parts[1]
            hits.setdefault(qid, sid)
    return hits


def _load_fasta_gene_map(fasta_path: Optional[str], gene_field: str) -> dict[str, str]:
    """
    Return {protein_id: gene_id} parsed from FASTA header lines.
    The gene_field parameter specifies which key=value field to extract,
    e.g. gene_field='gene' matches 'gene=NCU10129' → 'NCU10129'.
    """
    gene_map: dict[str, str] = {}
    if not fasta_path or not Path(fasta_path).exists():
        return gene_map
    pattern = re.compile(rf'{re.escape(gene_field)}=(\S+)')
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            pid = line[1:].split()[0]
            m = pattern.search(line)
            if m:
                gene_map[pid] = m.group(1)
    return gene_map
