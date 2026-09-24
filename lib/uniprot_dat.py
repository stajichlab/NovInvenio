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
import hashlib
import re
from pathlib import Path

from compressed_io import open_maybe_compressed

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
ID_RE = re.compile(r"^ID\s+(\S+)\s+(Reviewed|Unreviewed);")
PE_RE = re.compile(r"^PE\s+(\d)")
RX_PMID_RE = re.compile(r"PubMed=(\d+)")
RX_DOI_RE = re.compile(r"DOI=([^;\s]+)")
GENOME_RP = "LARGE SCALE GENOMIC DNA"
REF_LINE_TYPES = ("RP", "RX", "RT", "RA", "RG", "RL", "RC")

RECORD_COLUMNS = (
    "accession", "entry_name", "reviewed", "protein_existence", "taxon_id",
    "gene_name", "description", "ec_numbers", "go_ids", "pfam_ids", "pfam_names",
    "interpro_ids", "alphafold_id", "xrefs", "pubs", "seq_md5", "seq_len",
)

# Only RefSeq is deduped to its first DR line per record: RefSeq's protein-vs-
# transcript field ambiguity is a real correctness hazard. The other DBs here can
# legitimately repeat per record (different loci/paralogs sharing an accession) and
# all such entries are kept.
XREF_FIRST_FIELD_DBS = {
    "VEuPathDB", "GeneID", "KEGG", "RefSeq",
    # Added 2026-09-23 (UniProt library index spec): each of these carries its
    # stable ID in the first field, checked against UP000001805 / UP000002311.
    "PDB", "PANTHER", "OrthoDB", "STRING", "eggNOG", "Gene3D", "SUPFAM", "PROSITE",
    "SMART", "CDD", "PRINTS", "PIRSF", "HAMAP", "NCBIfam", "FunFam", "MEROPS",
    "CAZy", "ESTHER", "TCDB", "BRENDA", "UniPathway",
}
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


def normalize_seq(seq: str) -> str:
    """Upper case, no whitespace, trailing stop '*' removed -- the one form both
    the library index and the per-run lookup hash."""
    return "".join(seq.split()).upper().rstrip("*")


def seq_md5(seq: str) -> str:
    return hashlib.md5(normalize_seq(seq).encode("ascii", "replace")).hexdigest()


def _pack_ref(ref: dict) -> str:
    """One publication as "PMID;DOI;scope;title" (title last, may contain ';';
    '|' is the list separator so it becomes '/')."""
    title = " ".join(ref["rt"]).strip().strip('"').rstrip(";").rstrip('"').strip()
    title = title.replace("|", "/")
    scope = "proteome" if GENOME_RP in ref["rp"] else "protein"
    return f'{ref["pmid"]};{ref["doi"]};{scope};{title}'


def parse_dat_gz(path, stats=None):
    """Yield one dict per protein entry (keys = RECORD_COLUMNS) in a UniProt
    .dat(.gz/.zst) flat file. Entries with no AC or no sequence are skipped; when
    `stats` is a dict it receives n_records / n_skipped_no_ac / n_skipped_no_seq."""
    path = Path(path)
    if stats is not None:
        stats.update(n_records=0, n_skipped_no_ac=0, n_skipped_no_seq=0)
    st = {}

    def reset():
        st.update(accession=None, entry_name=None, reviewed=None, pe=None, taxon_id=None,
                  gene_name=None, rec_description=None, sub_description=None,
                  go_ids=[], pfam_ids=[], pfam_names=[], interpro_ids=[], ec_numbers=[],
                  alphafold_id=None, xrefs=[], xref_seen_dbs=set(), refs=[], cur_ref=None,
                  in_seq=False, seq_chunks=[])

    reset()
    with open_maybe_compressed(path) as fh:
        for line in fh:
            if line.startswith("//"):
                seq = normalize_seq("".join(st["seq_chunks"]))
                if not st["accession"]:
                    if stats is not None:
                        stats["n_skipped_no_ac"] += 1
                elif not seq:
                    if stats is not None:
                        stats["n_skipped_no_seq"] += 1
                else:
                    if stats is not None:
                        stats["n_records"] += 1
                    pubs = [_pack_ref(r) for r in st["refs"] if r["pmid"] or r["doi"]]
                    yield {
                        "accession": st["accession"],
                        "entry_name": st["entry_name"] or "",
                        "reviewed": st["reviewed"] or "0",
                        "protein_existence": st["pe"] or "",
                        "taxon_id": st["taxon_id"] or "",
                        "gene_name": st["gene_name"] or "",
                        "description": (st["rec_description"] or st["sub_description"] or "").strip(),
                        "ec_numbers": "|".join(st["ec_numbers"]),
                        "go_ids": "|".join(st["go_ids"]),
                        "pfam_ids": "|".join(st["pfam_ids"]),
                        "pfam_names": "|".join(st["pfam_names"]),
                        "interpro_ids": "|".join(st["interpro_ids"]),
                        "alphafold_id": st["alphafold_id"] or "",
                        "xrefs": "|".join(st["xrefs"]),
                        "pubs": "|".join(pubs),
                        "seq_md5": hashlib.md5(seq.encode("ascii", "replace")).hexdigest(),
                        "seq_len": str(len(seq)),
                    }
                reset()
                continue
            if st["in_seq"]:
                if line.startswith("     "):
                    st["seq_chunks"].append(line)
                continue
            code = line[:2]
            if code == "ID":
                m = ID_RE.match(line)
                if m:
                    st["entry_name"] = m.group(1)
                    st["reviewed"] = "1" if m.group(2) == "Reviewed" else "0"
                continue
            if code == "RN":
                st["cur_ref"] = {"pmid": "", "doi": "", "rp": "", "rt": []}
                st["refs"].append(st["cur_ref"])
                continue
            if code in REF_LINE_TYPES:
                ref = st["cur_ref"]
                if ref is not None:
                    body = line[5:].rstrip("\n")
                    if code == "RP":
                        ref["rp"] += " " + body
                    elif code == "RX":
                        m = RX_PMID_RE.search(body)
                        if m:
                            ref["pmid"] = m.group(1)
                        m = RX_DOI_RE.search(body)
                        if m:
                            ref["doi"] = m.group(1)
                    elif code == "RT":
                        ref["rt"].append(body.strip())
                continue
            st["cur_ref"] = None
            if code == "SQ":
                st["in_seq"] = True
                continue
            if code == "PE":
                m = PE_RE.match(line)
                if m:
                    st["pe"] = m.group(1)
                continue
            if st["accession"] is None and code == "AC":
                m = AC_RE.match(line)
                if m:
                    st["accession"] = m.group(1)
                continue
            if st["gene_name"] is None and code == "GN":
                m = GN_RE.search(line)
                if m:
                    st["gene_name"] = m.group(1).strip()
                continue
            if code == "DE":
                if st["rec_description"] is None:
                    m = DE_RECNAME_RE.match(line)
                    if m:
                        st["rec_description"] = m.group(1).strip()
                if st["sub_description"] is None:
                    m = DE_SUBNAME_RE.match(line)
                    if m:
                        st["sub_description"] = m.group(1).strip()
                m = DE_EC_RE.search(line)
                if m:
                    st["ec_numbers"].append(m.group(1))
                continue
            if st["taxon_id"] is None and code == "OX":
                m = OX_RE.search(line)
                if m:
                    st["taxon_id"] = m.group(1)
                continue
            if code == "DR":
                m = DR_GO_RE.match(line)
                if m:
                    st["go_ids"].append(f"{m.group(1)}:{m.group(2)}")
                    continue
                m = DR_PFAM_RE.match(line)
                if m:
                    st["pfam_ids"].append(m.group(1))
                    st["pfam_names"].append(m.group(2).strip())
                    continue
                m = DR_INTERPRO_RE.match(line)
                if m:
                    st["interpro_ids"].append(m.group(1))
                    continue
                if st["alphafold_id"] is None:
                    m = DR_ALPHAFOLD_RE.match(line)
                    if m:
                        st["alphafold_id"] = m.group(1)
                        continue
                m = DR_LINE_RE.match(line)
                if m:
                    db, rest = m.group(1), m.group(2)
                    if db in XREF_FIRST_FIELD_DBS:
                        if db == "RefSeq" and "RefSeq" in st["xref_seen_dbs"]:
                            continue
                        val = _first_xref_field(rest.split(";"))
                        if val:
                            st["xrefs"].append(f"{db}:{val}")
                            st["xref_seen_dbs"].add(db)
                        continue
                    if db in XREF_LAST_FIELD_DBS:
                        val = _last_xref_field(rest.split(";"))
                        if val:
                            st["xrefs"].append(f"{db}:{val}")
                        continue
