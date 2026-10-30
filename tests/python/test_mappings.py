"""Tests for the Mapping1D family exposed to Python.

Each mapping is a single object that BOTH draws a value (Forward) and
reports its own normalized density (Density), so a model/channel routing
both through one shared instance cannot let the sampler and the density
drift apart (Contract C1).  These tests verify the pybind bindings work,
that every concrete map is a polymorphic Mapping1D, that Forward/Inverse
round-trip, and that the density is positive and normalized to 1.
"""

import numpy as np
import pytest

# numpy.trapz was renamed to numpy.trapezoid in NumPy 2.0 and removed.
try:
    from numpy import trapezoid as _trapz   # NumPy >= 2.0
except ImportError:                          # NumPy < 2.0
    from numpy import trapz as _trapz

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
    ("Log",
     lambda: injection.LogMapping(0.1, 100.0), 0.1, 100.0),
    ("Exponential",
     lambda: injection.ExponentialMapping(2.0, 0.0, 10.0), 0.0, 10.0),
    ("Gaussian",
     lambda: injection.GaussianMapping(5.0, 2.0, 0.0, 10.0), 0.0, 10.0),
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
    integral = _trapz([m.Density(x) for x in grid], grid)
    assert integral == pytest.approx(1.0, abs=2e-2), \
        f"{name}: density integral = {integral}, expected 1"




def test_log_mapping_matches_powerlaw_limit():
    """LogMapping is the nu->1 case PowerLaw cannot represent (1/(1-nu) blows
    up); check it reproduces a near-1 PowerLaw within the shared range."""
    lo, hi = 0.5, 50.0
    log_map = injection.LogMapping(lo, hi)
    pl = injection.PowerLawMapping(0.999, 0.0, lo, hi)  # nu just below 1
    for r in np.linspace(0.05, 0.95, 19):
        assert log_map.Forward(r) == pytest.approx(pl.Forward(r), rel=2e-2)


