"""Tests for the Mapping1D family exposed to Python (Phase B).

Each mapping is a single object that BOTH draws a value (Forward) and
reports its own normalized density (Density), so a model/channel routing
both through one shared instance cannot let the sampler and the density
drift apart (Contract C1).  These tests verify the pybind bindings work,
that every concrete map is a polymorphic Mapping1D, that Forward/Inverse
round-trip, and that the density is positive and normalized to 1.
"""

import numpy as np
import pytest

from siren import injection


# (name, factory, lo, hi)
MAPS = [
    ("BreitWigner",
     lambda: injection.BreitWignerMapping(1.0, 0.1, 0.5, 2.0), 0.5, 2.0),
    ("PowerLaw",
     lambda: injection.PowerLawMapping(0.8, 0.0, 0.01, 100.0), 0.01, 100.0),
    ("Tabulated",
     lambda: injection.TabulatedMapping(
         [0.0, 1.0, 2.0, 3.0, 4.0], [0.0, 1.0, 4.0, 9.0, 16.0], 0.0, 4.0),
     0.0, 4.0),
    ("Propagator",
     lambda: injection.PropagatorMapping(0.01, 0.0, 10.0), 0.0, 10.0),
    ("Uniform",
     lambda: injection.UniformMapping(1.0, 10.0), 1.0, 10.0),
]
_IDS = [m[0] for m in MAPS]


@pytest.mark.parametrize("name,factory,lo,hi", MAPS, ids=_IDS)
def test_mapping_is_bound_and_polymorphic(name, factory, lo, hi):
    m = factory()
    assert isinstance(m, injection.Mapping1D)


@pytest.mark.parametrize("name,factory,lo,hi", MAPS, ids=_IDS)
def test_mapping_forward_in_range_and_roundtrips(name, factory, lo, hi):
    m = factory()
    span = hi - lo
    for r in np.linspace(0.0, 1.0, 21):
        x = m.Forward(r)
        assert lo - 1e-9 * span <= x <= hi + 1e-9 * span, \
            f"{name}: Forward({r}) = {x} outside [{lo}, {hi}]"
        r_back = m.Inverse(x)
        assert r_back == pytest.approx(r, abs=1e-6), \
            f"{name}: round-trip {r} -> {x} -> {r_back}"


@pytest.mark.parametrize("name,factory,lo,hi", MAPS, ids=_IDS)
def test_mapping_density_positive_and_normalized(name, factory, lo, hi):
    m = factory()
    # The normalized CDF runs 0 -> 1 across the range.
    assert m.Inverse(lo) == pytest.approx(0.0, abs=1e-9)
    assert m.Inverse(hi) == pytest.approx(1.0, abs=1e-9)

    xs = np.linspace(lo, hi, 2001)[1:-1]
    dens = np.array([m.Density(x) for x in xs])
    assert np.all(dens > 0.0), f"{name}: non-positive density encountered"

    grid = np.linspace(lo, hi, 40001)
    integral = np.trapz([m.Density(x) for x in grid], grid)
    assert integral == pytest.approx(1.0, abs=2e-2), \
        f"{name}: density integral = {integral}, expected 1"


def test_fixed_map_adaptive_hooks_are_noops():
    """The base's Accumulate/Refine hooks exist and are no-ops on fixed maps."""
    m = injection.UniformMapping(0.0, 1.0)
    d_before = m.Density(0.5)
    m.Accumulate(0.5, 1.0)
    m.Refine()
    assert m.Density(0.5) == d_before
