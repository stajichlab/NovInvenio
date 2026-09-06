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
