import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pfam_classes import CLASS_ORDER, classify_domain, dominant_class


def test_nlr_domains_classify_as_nlr():
    assert classify_domain("NACHT") == "nlr"
    assert classify_domain("HET") == "nlr"
    assert classify_domain("Ank_2") == "nlr"


def test_secondary_metabolite_domains():
    assert classify_domain("ketoacyl-synt") == "secondary_metabolite"
    assert classify_domain("AMP-binding") == "secondary_metabolite"
    assert classify_domain("Thioesterase") == "secondary_metabolite"


def test_transporter_domains():
    assert classify_domain("MFS_1") == "transporter"
    assert classify_domain("ABC_tran") == "transporter"


def test_matching_is_case_insensitive():
    # Pfam name casing is not stable across releases.
    assert classify_domain("nacht") == "nlr"
    assert classify_domain("Ketoacyl-Synt") == "secondary_metabolite"


def test_unrecognised_domain_is_other():
    assert classify_domain("DUF1653") == "other"


def test_empty_or_missing_domain_is_unannotated():
    assert classify_domain("") == "unannotated"
    assert classify_domain("-") == "unannotated"


def test_dominant_class_is_most_frequent():
    assert dominant_class(["MFS_1", "ABC_tran", "NACHT"]) == "transporter"


def test_dominant_class_breaks_ties_by_precedence():
    # One each: precedence decides, and it must not depend on input order.
    assert dominant_class(["MFS_1", "NACHT"]) == "nlr"
    assert dominant_class(["NACHT", "MFS_1"]) == "nlr"
    assert CLASS_ORDER.index("nlr") < CLASS_ORDER.index("transporter")


def test_dominant_class_of_nothing_is_unannotated():
    assert dominant_class([]) == "unannotated"
    assert dominant_class(["-"]) == "unannotated"
