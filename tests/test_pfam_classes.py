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


def test_het_false_positives_reject_thet_phet():
    # het must be anchored to prefix to avoid matching "thet" and "phet" substrings
    # (dna_pol3_theta, phetrs_b1, rbd_dgktheta are false positives under substring match)
    assert classify_domain("dna_pol3_theta") == "other"
    assert classify_domain("phetrs_b1") == "other"
    assert classify_domain("rbd_dgktheta") == "other"


def test_het_true_positives_still_classify_as_nlr():
    # Genuine het-family domains must remain classified as nlr
    assert classify_domain("het") == "nlr"
    assert classify_domain("het-s") == "nlr"
    assert classify_domain("het-c") == "nlr"
    assert classify_domain("hetr_c") == "nlr"


def test_tpr_false_positives_no_longer_nlr():
    # tpr (tetratricopeptide repeat) is a generic protein-protein interaction
    # scaffold found in hundreds of unrelated proteins; it is NOT NLR-specific.
    # These must classify as "other", not "nlr".
    assert classify_domain("apc5_tpr") == "other"
    assert classify_domain("cnot10_tpr") == "other"


def test_mfs_covers_all_variants():
    # "mfs_1" is redundant; "mfs_" already matches mfs_1, mfs_2..5, mfs_mot1, etc.
    assert classify_domain("mfs_mot1") == "transporter"
    assert classify_domain("mfs_mycoplasma") == "transporter"


def test_ank_still_classifies_as_nlr():
    # ank (ankyrin repeat) is genuinely NLR-related; verified against 25 Pfam names.
    assert classify_domain("ank_2") == "nlr"
    assert classify_domain("ankrd13_c") == "nlr"
