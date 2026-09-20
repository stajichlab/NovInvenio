"""Map a Pfam domain NAME to a display class for the island synteny view.

This is data, not logic: one curated keyword -> class table, matched
case-insensitively as a substring of the domain name, because Pfam names
carry version/variant suffixes (Ank_2, Ank_4, ...) and casing is not stable
across Pfam releases.

A family whose domains fall in two classes takes the first match in
CLASS_ORDER, which is a documented precedence, not an accident of dict
iteration order. The order is most-specific-biology first: an island
carrying both an NLR domain and a transporter is far more interesting as an
NLR island, and transporter domains are the more promiscuous of the two.
"""
from __future__ import annotations

UNANNOTATED = "unannotated"
OTHER = "other"

# Precedence, most specific first. OTHER and UNANNOTATED are terminal.
CLASS_ORDER: tuple[str, ...] = (
    "nlr",
    "secondary_metabolite",
    "transporter",
    OTHER,
    UNANNOTATED,
)

CLASS_LABELS: dict[str, str] = {
    "nlr": "NLR / incompatibility",
    "secondary_metabolite": "Secondary metabolite",
    "transporter": "Transporter",
    OTHER: "Other",
    UNANNOTATED: "No Pfam hit",
}

# keyword (lowercase, substring-matched) -> class
_KEYWORDS: dict[str, str] = {
    "nacht": "nlr",
    "het": "nlr",
    "ank": "nlr",
    "tpr": "nlr",
    "ketoacyl-synt": "secondary_metabolite",
    "amp-binding": "secondary_metabolite",
    "pp-binding": "secondary_metabolite",
    "thioesterase": "secondary_metabolite",
    "acp_": "secondary_metabolite",
    "ps-dh": "secondary_metabolite",
    "mfs_": "transporter",
    "sugar_tr": "transporter",
    "abc_tran": "transporter",
    "mfs_1": "transporter",
}

_MISSING = {"", "-", "na", "none"}


def classify_domain(name: str) -> str:
    """Class for one Pfam domain name. Missing/sentinel -> UNANNOTATED."""
    if name is None:
        return UNANNOTATED
    key = name.strip().lower()
    if key in _MISSING:
        return UNANNOTATED
    matches = {cls for kw, cls in _KEYWORDS.items() if kw in key}
    if not matches:
        return OTHER
    for cls in CLASS_ORDER:
        if cls in matches:
            return cls
    return OTHER


def dominant_class(names: list[str]) -> str:
    """Most frequent class among `names`, ties broken by CLASS_ORDER.

    OTHER and UNANNOTATED never win over a real class: they are what you get
    when there is nothing better, so they are only returned if no name in the
    list classified into a real class.
    """
    classes = [classify_domain(n) for n in (names or [])]
    real = [c for c in classes if c not in (OTHER, UNANNOTATED)]
    pool = real or classes
    if not pool:
        return UNANNOTATED
    counts: dict[str, int] = {}
    for c in pool:
        counts[c] = counts.get(c, 0) + 1
    best = max(counts.values())
    tied = [c for c, n in counts.items() if n == best]
    for cls in CLASS_ORDER:
        if cls in tied:
            return cls
    return UNANNOTATED
