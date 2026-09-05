"""
Guards for the shared colour-skin registry (lib/skins.py).

The point of these tests is that adding a skin -- especially a high-chroma one
like ``neuromancer`` -- cannot quietly regress legibility or drop a token that
some report happens to use. Without them "every skin stays readable" is an
aspiration; with them it is a build failure.
"""
import re

import pytest
from skins import (
    DEFAULT_DARK_SKIN,
    DEFAULT_SKIN,
    MIN_ON_SERIES_CONTRAST,
    MIN_SERIES_CONTRAST,
    MIN_TEXT_CONTRAST,
    REQUIRED_TOKENS,
    SKINS,
    contrast_ratio,
    relative_luminance,
    skin_boot_js,
    skin_picker_html,
    skin_picker_js,
    skins_css,
)

SKIN_IDS = sorted(SKINS)


def test_defaults_are_real_skins():
    assert DEFAULT_SKIN in SKINS
    assert DEFAULT_DARK_SKIN in SKINS


@pytest.mark.parametrize('skin_id', SKIN_IDS)
def test_skin_defines_every_token(skin_id):
    """A missing token would silently inherit from the previously applied skin."""
    assert set(SKINS[skin_id]['tokens']) == set(REQUIRED_TOKENS)


@pytest.mark.parametrize('skin_id', SKIN_IDS)
def test_skin_has_label_and_scheme(skin_id):
    skin = SKINS[skin_id]
    assert skin['label']
    assert skin['scheme'] in ('light', 'dark')


@pytest.mark.parametrize('skin_id', SKIN_IDS)
def test_body_text_meets_wcag_aa(skin_id):
    tokens = SKINS[skin_id]['tokens']
    for name in ('--text-primary', '--text-secondary'):
        ratio = contrast_ratio(tokens[name], tokens['--page'])
        assert ratio >= MIN_TEXT_CONTRAST, (
            f'{skin_id} {name} on --page is {ratio:.2f}:1, '
            f'below the {MIN_TEXT_CONTRAST}:1 floor'
        )


@pytest.mark.parametrize('skin_id', SKIN_IDS)
def test_data_marks_meet_wcag_aa(skin_id):
    """--series-1/-2 and --warn are graphics, so the 3:1 UI-component floor applies."""
    tokens = SKINS[skin_id]['tokens']
    for name in ('--series-1', '--series-2', '--warn', '--muted'):
        ratio = contrast_ratio(tokens[name], tokens['--page'])
        assert ratio >= MIN_SERIES_CONTRAST, (
            f'{skin_id} {name} on --page is {ratio:.2f}:1, '
            f'below the {MIN_SERIES_CONTRAST}:1 floor'
        )


@pytest.mark.parametrize('skin_id', SKIN_IDS)
def test_series_pair_is_luminance_separated(skin_id):
    """The two evidence colours must not collapse into each other.

    What actually carries this pair under red-green CVD is *hue* separation on
    the blue-yellow axis (see lib/skins.py's colour-role contract), and no
    luminance test can check that. This is the weaker companion guard: it
    catches a pair that has collapsed outright, e.g. two neons a skin author
    picked for looks without noticing they render as the same grey.

    The floor is 1.10 because that is where the *shipped* palettes already
    sit -- paper's validated blue/green pair is 1.12 and dark's is 1.36 -- so
    the assertion is "no worse than what already ships", not an ideal. Raising
    it would mean recolouring paper, which CLAUDE.md requires re-validating
    for colour-blind separation first.
    """
    tokens = SKINS[skin_id]['tokens']
    l1 = relative_luminance(tokens['--series-1'])
    l2 = relative_luminance(tokens['--series-2'])
    lo, hi = sorted((l1, l2))
    ratio = (hi + 0.05) / (lo + 0.05)
    assert ratio >= 1.10, (
        f'{skin_id} --series-1/--series-2 differ by only {ratio:.2f}:1 in '
        'luminance; they will merge in greyscale'
    )


@pytest.mark.parametrize('skin_id', SKIN_IDS)
def test_skin_emits_a_css_block(skin_id):
    css = skins_css()
    assert f':root[data-skin="{skin_id}"]' in css
    for token in REQUIRED_TOKENS:
        assert f'{token}:' in css


def test_css_has_system_default_and_print_override():
    css = skins_css()
    assert '@media (prefers-color-scheme: dark)' in css
    assert ':root:where(:not([data-skin]))' in css
    assert '@media print' in css
    # The print block must pin the light palette whatever is on screen.
    print_block = css.split('@media print', 1)[1]
    assert SKINS[DEFAULT_SKIN]['tokens']['--page'] in print_block


def test_picker_lists_every_skin_plus_follow_system():
    html = skin_picker_html()
    assert 'value=""' in html          # follow-system option
    for skin_id, skin in SKINS.items():
        assert f'value="{skin_id}"' in html
        assert skin['label'] in html
    # One <option> per skin, plus follow-system.
    assert html.count('<option') == len(SKINS) + 1


def test_boot_and_picker_js_agree_on_the_valid_skin_list():
    boot = skin_boot_js()
    for skin_id in SKINS:
        assert f'"{skin_id}"' in boot
    # Both halves must key off the same storage slot.
    key = re.search(r'localStorage\.getItem\("([^"]+)"\)', boot).group(1)
    assert f'localStorage.setItem("{key}"' in skin_picker_js()


def test_storage_access_is_guarded():
    """A file:// document can have an opaque origin where localStorage throws."""
    for js in (skin_boot_js(), skin_picker_js()):
        assert 'try {' in js and 'catch' in js


def test_contrast_ratio_reference_values():
    assert contrast_ratio('#000000', '#ffffff') == pytest.approx(21.0, abs=0.01)
    assert contrast_ratio('#ffffff', '#ffffff') == pytest.approx(1.0, abs=0.01)
    # Order must not matter.
    assert contrast_ratio('#2a78d6', '#f9f9f7') == pytest.approx(
        contrast_ratio('#f9f9f7', '#2a78d6')
    )


@pytest.mark.parametrize('skin_id', SKIN_IDS)
def test_on_series_label_is_readable(skin_id):
    """--on-series is drawn on top of a --series-* fill (the presence chips).

    This is the check that caught white-on-neon: a hardcoded #fff label sat at
    ~1.5:1 on the bright cyan and amber of the neon skin.
    """
    tokens = SKINS[skin_id]['tokens']
    for series in ('--series-1', '--series-2'):
        ratio = contrast_ratio(tokens['--on-series'], tokens[series])
        assert ratio >= MIN_ON_SERIES_CONTRAST, (
            f'{skin_id} --on-series on {series} is {ratio:.2f}:1, '
            f'below the {MIN_ON_SERIES_CONTRAST}:1 floor'
        )
