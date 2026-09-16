"""Shared data structures for the pangenome cluster-profile analysis
(notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md).

PresenceMatrix holds one row per gene family (a tier-1 mmseqs/diamond
cluster representative) and one column per strain, with a three-state
call per cell -- see the Global Constraints in this plan's header for
what each state means.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from compressed_io import open_maybe_compressed, open_maybe_compressed_write

PRESENT = "present"
GENOME_ONLY = "genome_only"
ABSENT = "absent"
STATES = (PRESENT, GENOME_ONLY, ABSENT)


def read_cluster_tsv(path: str | Path) -> dict[str, str]:
    """Parse an mmseqs/diamond cluster TSV (rep\tmember per line) into
    {member_id: rep_id}."""
    member_to_rep: dict[str, str] = {}
    with open_maybe_compressed(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            member_to_rep[member] = rep
    return member_to_rep


def build_families(member_to_rep: dict[str, str]) -> dict[str, list[str]]:
    """Return {rep_id: sorted [member_ids]} for every cluster, including
    singletons -- unlike NovInvenio's lib/clusters.py::build_families,
    which drops singleton clusters (irrelevant there; here a singleton
    family is a real, reportable frequency bin)."""
    families: dict[str, list[str]] = {}
    for member, rep in member_to_rep.items():
        families.setdefault(rep, []).append(member)
    return {rep: sorted(members) for rep, members in families.items()}


@dataclass
class PresenceMatrix:
    families: list[str]
    strains: list[str]
    calls: dict[tuple[str, str], str] = field(default_factory=dict)
    copy_number: dict[tuple[str, str], int] = field(default_factory=dict)

    def set_call(self, family: str, strain: str, state: str, copies: int = 0) -> None:
        if state not in STATES:
            raise ValueError(f"invalid state {state!r}, must be one of {STATES}")
        self.calls[(family, strain)] = state
        if copies:
            self.copy_number[(family, strain)] = copies

    def call(self, family: str, strain: str) -> str:
        return self.calls.get((family, strain), ABSENT)

    def is_present(self, family: str, strain: str) -> bool:
        return self.call(family, strain) in (PRESENT, GENOME_ONLY)

    def presence_vector(self, family: str) -> list[bool]:
        return [self.is_present(family, s) for s in self.strains]

    def strain_count(self, family: str) -> int:
        return sum(self.presence_vector(family))

    def frequency(self, family: str) -> float:
        if not self.strains:
            return 0.0
        return self.strain_count(family) / len(self.strains)

    def to_tsv(self, path: str | Path, write_copy_number: bool = True) -> None:
        """Write the family x strain state matrix.

        The main TSV holds ONLY the three-state call string per cell -- that
        format is the stable contract every downstream script reads, so copy
        numbers are persisted in a SIDECAR file (see `copy_number_path`)
        rather than by changing what a cell looks like. The sidecar is written
        only when this matrix actually carries copy numbers, so a state-only
        matrix still produces exactly one output file, as before.

        `path` may end in `.gz`/`.zst` to write compressed output -- the
        sidecar path (`copy_number_path`) preserves whatever suffix `path`
        ends in, so a compressed matrix gets a compressed sidecar too, per
        this workspace's storage-compression convention for large TSV
        intermediates (a full family x strain matrix at pangenome scale is
        exactly the kind of file that convention targets).
        """
        with open_maybe_compressed_write(path) as fh:
            fh.write("family\t" + "\t".join(self.strains) + "\n")
            for fam in self.families:
                row = [self.call(fam, s) for s in self.strains]
                fh.write(fam + "\t" + "\t".join(row) + "\n")
        if write_copy_number and self.copy_number:
            with open_maybe_compressed_write(copy_number_path(path)) as fh:
                fh.write("family\tstrain\tcopies\n")
                for (fam, strain), copies in sorted(self.copy_number.items()):
                    fh.write(f"{fam}\t{strain}\t{copies}\n")

    @classmethod
    def from_tsv(cls, path: str | Path, read_copy_number: bool = True) -> "PresenceMatrix":
        """Load a matrix written by `to_tsv`.

        The copy-number sidecar is loaded when it exists; its absence is never
        an error (a matrix with no copy numbers, or one written before sidecars
        existed, simply loads with every copy number 0).

        Raises:
            ValueError: if any cell holds a value that is not one of `STATES`
                -- a corrupted or hand-edited matrix, which would otherwise
                load silently and evaluate as "not present" everywhere.
        """
        with open_maybe_compressed(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            strains = header[1:]
            families: list[str] = []
            calls: dict[tuple[str, str], str] = {}
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                fam = parts[0]
                families.append(fam)
                for strain, state in zip(strains, parts[1:]):
                    if state not in STATES:
                        raise ValueError(
                            f"{path}: invalid presence state {state!r} for family "
                            f"{fam!r} / strain {strain!r}; must be one of {STATES}"
                        )
                    calls[(fam, strain)] = state
        pm = cls(families=families, strains=strains)
        pm.calls = calls
        if read_copy_number:
            pm.copy_number = read_copy_number_sidecar(path)
        return pm


_COMPRESSION_SUFFIXES = (".gz", ".zst")


def copy_number_path(matrix_path: str | Path) -> Path:
    """The copy-number sidecar path for a given presence-matrix path.

    A compressed matrix path (`matrix.tsv.zst`) gets a compressed sidecar
    path too (`matrix.copy_number.tsv.zst`), not `matrix.tsv.zst.copy_number.tsv`
    -- the compression suffix must stay LAST so `open_maybe_compressed`/
    `open_maybe_compressed_write`'s extension check still recognizes it.
    """
    matrix_path = str(matrix_path)
    for suffix in _COMPRESSION_SUFFIXES:
        if matrix_path.endswith(suffix):
            return Path(matrix_path[: -len(suffix)] + ".copy_number.tsv" + suffix)
    return Path(matrix_path + ".copy_number.tsv")


def read_copy_number_sidecar(matrix_path: str | Path) -> dict[tuple[str, str], int]:
    """Load the copy-number sidecar (see `copy_number_path`) if it exists;
    return {} if it does not.

    Tries `copy_number_path(matrix_path)` itself first (already
    compression-suffixed to match `matrix_path`, e.g. `.copy_number.tsv.zst`
    for a `.zst` matrix -- what `PresenceMatrix.to_tsv` actually writes),
    then that same path with an extra `.gz`/`.zst` appended -- a defensive
    fallback for a sidecar that was compressed after the fact without its
    name being recomputed, not the shape `to_tsv` itself produces."""
    sidecar = copy_number_path(matrix_path)
    candidates = [sidecar, Path(str(sidecar) + ".gz"), Path(str(sidecar) + ".zst")]
    sidecar = next((c for c in candidates if c.exists()), None)
    if sidecar is None:
        return {}
    copies: dict[tuple[str, str], int] = {}
    with open_maybe_compressed(sidecar) as fh:
        for i, line in enumerate(fh):
            line = line.rstrip("\n")
            if not line:
                continue
            if i == 0 and line.startswith("family\t"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            copies[(parts[0], parts[1])] = int(parts[2])
    return copies
