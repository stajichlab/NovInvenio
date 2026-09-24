"""Reader for a UniProt library index built by bin/uniprot_build_index.py
(docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md).

Index layout (format 1): manifest.json, seq_index.sqlite (tables seq, acc, refseq,
proteome), records/<proteome_id>.records.tsv.zst (columns = uniprot_dat.RECORD_COLUMNS).
"""
import csv
import json
import sqlite3
from pathlib import Path

from compressed_io import open_maybe_compressed

FORMAT_VERSION = 1


class IndexFormatError(Exception):
    """The index directory is missing a file or has an unsupported format."""


def _binomial(name):
    """First two words, lower case -- 'Neurospora crassa (strain ...)' and
    'Neurospora crassa OR74A' both give 'neurospora crassa'."""
    return " ".join((name or "").lower().replace("(", " ").split()[:2])


class UniProtIndex:
    def __init__(self, index_dir):
        self.dir = Path(index_dir)
        for need in ("manifest.json", "seq_index.sqlite", "records"):
            if not (self.dir / need).exists():
                raise IndexFormatError(f"UniProt index {self.dir} is incomplete: {need} missing")
        self.manifest = json.loads((self.dir / "manifest.json").read_text())
        v = self.manifest.get("format_version")
        if v != FORMAT_VERSION:
            raise IndexFormatError(f"UniProt index {self.dir} has format_version {v}; "
                                   f"this code reads {FORMAT_VERSION}")
        self.con = sqlite3.connect(f"file:{self.dir / 'seq_index.sqlite'}?mode=ro", uri=True)
        self._species = {pid: (tax, name) for pid, tax, name in self.con.execute(
            "select proteome_id, taxid, species_name from proteome")}

    def by_accession(self, acc):
        r = self.con.execute("select accession, proteome_id, taxid from acc where accession=?",
                             (acc,)).fetchone()
        return dict(zip(("accession", "proteome_id", "taxid"), r)) if r else None

    def by_refseq(self, refseq_id):
        return [r[0] for r in self.con.execute(
            "select accession from refseq where refseq_id=? order by accession", (refseq_id,))]

    def by_md5(self, md5):
        return [dict(zip(("accession", "proteome_id", "taxid", "reviewed"), r))
                for r in self.con.execute(
                    "select accession, proteome_id, taxid, reviewed from seq where md5=?", (md5,))]

    def taxid_for_species(self, name):
        """Taxid of the library proteome(s) whose species binomial matches, when exactly
        one taxid matches; else None."""
        want = _binomial(name)
        taxids = {tax for tax, sp in self._species.values() if _binomial(sp) == want}
        return taxids.pop() if len(taxids) == 1 else None

    def taxids_for_species(self, name):
        """Every library taxid whose species binomial matches (strain-level taxids of the
        same species all count), as a set; empty when none match."""
        want = _binomial(name)
        return {tax for tax, sp in self._species.values() if _binomial(sp) == want}

    def species_name(self, proteome_id):
        return self._species.get(proteome_id, (None, ""))[1]

    def records(self, proteome_id, accessions):
        """accession -> record dict, for the requested accessions of one proteome."""
        out = {}
        path = self.dir / "records" / f"{proteome_id}.records.tsv.zst"
        with open_maybe_compressed(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["accession"] in accessions:
                    out[row["accession"]] = row
        return out
