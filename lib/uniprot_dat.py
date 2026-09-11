"""Parse a UniProt SwissProt/TrEMBL flat file (.dat.gz) into per-protein annotation.

Ported from NovInvenio_Investigations' bin/extract_dat_annotations.py (not imported
from it -- this repo must not have a runtime dependency on a private sibling repo).
Kept in lock-step with that script's DR-line extraction rules; see
docs/superpowers/specs/2026-09-09-uniprot-xref-linkout-design.md for the full design
rationale (allow-listed DB set, packing rule, per-DB field quirks).

Deliberately not using Biopython's SwissProt parser: that parses the entire record
(references, comments, sequence, everything) to get at a handful of line types. This
reads only AC/GN/DE/OX/DR lines and resets on the `//` record separator.

Yielded dict per record: accession, taxon_id, gene_name, description, go_ids,
pfam_ids, pfam_names, interpro_ids, ec_numbers, alphafold_id, xrefs.
  - go_ids / pfam_ids / pfam_names / interpro_ids / ec_numbers are "|"-separated.
  - xrefs is a "|"-separated set of allow-listed cross-reference databases
    (VEuPathDB/GeneID/RefSeq/KEGG/EnsemblFungi/EnsemblBacteria), packed as `DB:id`
    pairs (first colon only separates DB from id -- id itself may contain colons,
    e.g. "KEGG:ncr:NCU10683"). RefSeq is deduped to its first DR line per record
    (protein-accession field only, never the paired transcript accession).
"""
import gzip
import re
from pathlib import Path

AC_RE = re.compile(r"^AC\s+([A-Z0-9]+)")
OX_RE = re.compile(r"NCBI_TaxID=(\d+)")
GN_RE = re.compile(r"(?:Name|ORFNames)=([^;{]+)")
DE_RECNAME_RE = re.compile(r"^DE\s+RecName:\s*Full=([^;{]+)")
DE_SUBNAME_RE = re.compile(r"^DE\s+SubName:\s*Full=([^;{]+)")
DR_GO_RE = re.compile(r"^DR\s+GO;\s*(GO:\d+);\s*[A-Z]:[^;]*;\s*([A-Za-z0-9_]+):")
DR_PFAM_RE = re.compile(r"^DR\s+Pfam;\s*(PF\d+);\s*([^;]+);")
DR_INTERPRO_RE = re.compile(r"^DR\s+InterPro;\s*(IPR\d+);")
DR_ALPHAFOLD_RE = re.compile(r"^DR\s+AlphaFoldDB;\s*([A-Z0-9]+);")
DE_EC_RE = re.compile(r"EC=([\d.]+(?:-)?)")

# Only RefSeq is deduped to its first DR line per record: RefSeq's protein-vs-
# transcript field ambiguity is a real correctness hazard. The other DBs here can
# legitimately repeat per record (different loci/paralogs sharing an accession) and
# all such entries are kept.
XREF_FIRST_FIELD_DBS = {"VEuPathDB", "GeneID", "KEGG", "RefSeq"}
XREF_LAST_FIELD_DBS = {"EnsemblFungi", "EnsemblBacteria"}
# Some Ensembl* DR lines carry a trailing isoform bracket tag after the final period,
# e.g. "DR   EnsemblFungi; YML032C_mRNA; YML032C; YML032C. [P06778-1]" -- strip that too.
DR_LINE_RE = re.compile(r"^DR\s+(\w+);\s*(.*?)\.?(?:\s*\[[^\]]*\])?\s*$")


def _first_xref_field(fields: list[str]) -> str:
    return fields[0].strip() if fields else ""


def _last_xref_field(fields: list[str]) -> str:
    # Field position/count for Ensembl* DR lines is not stable across organisms
    # (Ncra: protein;protein;gene -- Scer: transcript;protein;gene), but the last
    # non-placeholder field is consistently the stable gene ID.
    cleaned = [f.strip() for f in fields if f.strip() and f.strip() != "-"]
    return cleaned[-1] if cleaned else ""


def parse_dat_gz(path):
    """Yield one dict per protein entry in a UniProt .dat(.gz) flat file."""
    path = Path(path)
    accession = None
    taxon_id = None
    gene_name = None
    rec_description = None
    sub_description = None
    go_ids = []
    pfam_ids = []
    pfam_names = []
    interpro_ids = []
    ec_numbers = []
    alphafold_id = None
    xrefs = []
    xref_seen_dbs = set()

    def reset():
        nonlocal accession, taxon_id, gene_name, rec_description, sub_description
        nonlocal go_ids, pfam_ids, pfam_names, interpro_ids, ec_numbers, alphafold_id
        nonlocal xrefs, xref_seen_dbs
        accession = taxon_id = gene_name = rec_description = sub_description = None
        go_ids, pfam_ids, pfam_names, interpro_ids, ec_numbers = [], [], [], [], []
        alphafold_id = None
        xrefs = []
        xref_seen_dbs = set()

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("//"):
                if accession:
                    yield {
                        "accession": accession,
                        "taxon_id": taxon_id or "",
                        "gene_name": gene_name or "",
                        "description": (rec_description or sub_description or "").strip(),
                        "go_ids": "|".join(go_ids),
                        "pfam_ids": "|".join(pfam_ids),
                        "pfam_names": "|".join(pfam_names),
                        "interpro_ids": "|".join(interpro_ids),
                        "ec_numbers": "|".join(ec_numbers),
                        "alphafold_id": alphafold_id or "",
                        "xrefs": "|".join(xrefs),
                    }
                reset()
                continue
            if accession is None and line.startswith("AC"):
                m = AC_RE.match(line)
                if m:
                    accession = m.group(1)
                continue
            if gene_name is None and line.startswith("GN"):
                m = GN_RE.search(line)
                if m:
                    gene_name = m.group(1).strip()
                continue
            if line.startswith("DE"):
                if rec_description is None:
                    m = DE_RECNAME_RE.match(line)
                    if m:
                        rec_description = m.group(1).strip()
                if sub_description is None:
                    m = DE_SUBNAME_RE.match(line)
                    if m:
                        sub_description = m.group(1).strip()
                m = DE_EC_RE.search(line)
                if m:
                    ec_numbers.append(m.group(1))
                continue
            if taxon_id is None and line.startswith("OX"):
                m = OX_RE.search(line)
                if m:
                    taxon_id = m.group(1)
                continue
            if line.startswith("DR"):
                m = DR_GO_RE.match(line)
                if m:
                    go_ids.append(f"{m.group(1)}:{m.group(2)}")
                    continue
                m = DR_PFAM_RE.match(line)
                if m:
                    pfam_ids.append(m.group(1))
                    pfam_names.append(m.group(2).strip())
                    continue
                m = DR_INTERPRO_RE.match(line)
                if m:
                    interpro_ids.append(m.group(1))
                    continue
                if alphafold_id is None:
                    m = DR_ALPHAFOLD_RE.match(line)
                    if m:
                        alphafold_id = m.group(1)
                        continue
                m = DR_LINE_RE.match(line)
                if m:
                    db, rest = m.group(1), m.group(2)
                    if db in XREF_FIRST_FIELD_DBS:
                        if db == "RefSeq" and "RefSeq" in xref_seen_dbs:
                            continue
                        val = _first_xref_field(rest.split(";"))
                        if val:
                            xrefs.append(f"{db}:{val}")
                            xref_seen_dbs.add(db)
                        continue
                    if db in XREF_LAST_FIELD_DBS:
                        val = _last_xref_field(rest.split(";"))
                        if val:
                            xrefs.append(f"{db}:{val}")
                        continue
