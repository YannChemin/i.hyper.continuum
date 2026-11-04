"""Regression tests for the convex-hull continuum of i.hyper.continuum.

The module is a GRASS script: it imports ``grass.script`` at module scope and
carries a parser header, so it cannot be imported outside a session. The two
continuum functions are pure numpy, so they are lifted out of the source here
and tested on their own.

The property that matters is that the hull can descend. An envelope built
only from rising points stops at the reflectance maximum and extrapolates a
continuum far above the spectrum over the whole falling limb, which on a
VNIR-SWIR cube is the entire SWIR.

SPDX-License-Identifier: GPL-2.0-or-later
"""

import os

import numpy as np
import pytest

MODULE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "i.hyper.continuum.py")


def _load_continuum_functions():
    """Return the module's pure continuum functions, without GRASS."""
    with open(MODULE, encoding="utf-8") as handle:
        source = handle.read()
    namespace = {"np": np}
    for name in ("calculate_linear_continuum",
                 "calculate_convex_hull_continuum"):
        start = source.index(f"def {name}(")
        end = source.index("\ndef ", start + 5)
        # exec of this repository's own source, sliced to two pure functions:
        # the alternative is importing a GRASS script outside a session.
        exec(compile(source[start:end], MODULE, "exec"), namespace)  # noqa: S102
    return namespace["calculate_convex_hull_continuum"]


convex_hull_continuum = _load_continuum_functions()


@pytest.fixture
def wavelengths():
    """A VNIR-SWIR axis of the shape this module is used on."""
    return np.linspace(400.0, 2500.0, 200)


def _peaked_spectrum(wl):
    """Reflectance that rises to a NIR peak and falls away through the SWIR."""
    return 0.1 + 0.5 * np.exp(-((wl - 900.0) / 350.0) ** 2)


def test_continuum_never_below_the_spectrum(wavelengths):
    """The defining property: a continuum-removed spectrum cannot exceed 1."""
    refl = _peaked_spectrum(wavelengths)
    continuum = convex_hull_continuum(wavelengths, refl)
    assert np.all(continuum >= refl - 1e-9)
    assert np.all(refl / continuum <= 1.0 + 1e-9)


def test_hull_follows_the_falling_limb(wavelengths):
    """The hull must descend after the peak rather than extrapolate flat.

    This is the regression: an envelope of monotonically rising points held
    the continuum at the 900 nm peak across the whole SWIR, which divided
    every SWIR band by a value far above it.
    """
    refl = _peaked_spectrum(wavelengths)
    continuum = convex_hull_continuum(wavelengths, refl)
    peak = float(refl.max())
    assert continuum[-1] < 0.5 * peak, (
        "continuum did not descend after the peak: "
        f"{continuum[-1]:.4f} against a peak of {peak:.4f}")
    # It tracks the spectrum closely where the spectrum is itself convex.
    assert continuum[-1] == pytest.approx(refl[-1], rel=1e-6)


def test_absorption_depth_is_measured_against_the_hull(wavelengths):
    """A known absorption on a straight baseline returns its own depth.

    The baseline is linear, so it is its own upper hull and the only thing
    the continuum removes is the dip. On a curved spectrum the hull is a
    chord above the curve, and the measured depth then includes the
    spectrum's own departure from that chord -- which is the definition, not
    an error, and is why the baseline here is straight.
    """
    refl = np.linspace(0.60, 0.20, wavelengths.size)
    centre = int(np.argmin(np.abs(wavelengths - 2200.0)))
    refl[centre - 2:centre + 3] *= 0.70
    continuum = convex_hull_continuum(wavelengths, refl)
    depth = 1.0 - float(refl[centre] / continuum[centre])
    assert depth == pytest.approx(0.30, abs=0.01)


def test_endpoints_are_hull_vertices(wavelengths):
    """Nothing is extrapolated: the continuum meets the spectrum at both ends."""
    refl = _peaked_spectrum(wavelengths)
    continuum = convex_hull_continuum(wavelengths, refl)
    assert continuum[0] == pytest.approx(refl[0], rel=1e-9)
    assert continuum[-1] == pytest.approx(refl[-1], rel=1e-9)


@pytest.mark.parametrize("shape", ["flat", "rising", "falling"])
def test_monotonic_and_flat_spectra(wavelengths, shape):
    """A spectrum with no absorption has a continuum equal to itself."""
    refl = {"flat": np.ones_like(wavelengths),
            "rising": np.linspace(0.1, 0.6, wavelengths.size),
            "falling": np.linspace(0.6, 0.1, wavelengths.size)}[shape]
    continuum = convex_hull_continuum(wavelengths, refl)
    assert continuum == pytest.approx(refl, abs=1e-9)


def test_non_finite_bands_are_skipped(wavelengths):
    """A null band is not a hull vertex and does not void the continuum."""
    refl = _peaked_spectrum(wavelengths)
    refl[::13] = np.nan
    continuum = convex_hull_continuum(wavelengths, refl)
    good = np.isfinite(refl)
    assert np.all(np.isfinite(continuum))
    assert np.all(continuum[good] >= refl[good] - 1e-9)


def test_degenerate_input_falls_back():
    """Fewer than two usable bands cannot define a hull."""
    out = convex_hull_continuum(np.array([400.0]), np.array([0.2]))
    assert np.isfinite(out).all()
    out = convex_hull_continuum(np.array([400.0, 2500.0]),
                                np.array([np.nan, np.nan]))
    assert out.size == 2
