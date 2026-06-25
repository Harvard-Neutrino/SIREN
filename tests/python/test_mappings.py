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
    ("Log",
     lambda: injection.LogMapping(0.1, 100.0), 0.1, 100.0),
    ("Exponential",
     lambda: injection.ExponentialMapping(2.0, 0.0, 10.0), 0.0, 10.0),
    ("Gaussian",
     lambda: injection.GaussianMapping(5.0, 2.0, 0.0, 10.0), 0.0, 10.0),
    # Adaptive starts uniform (a no-op until refined), so it satisfies the
    # generic round-trip / normalization checks in this fresh state.
    ("Adaptive",
     lambda: injection.AdaptiveMapping(0.0, 1.0, 32), 0.0, 1.0),
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


def test_log_mapping_matches_powerlaw_limit():
    """LogMapping is the nu->1 case PowerLaw cannot represent (1/(1-nu) blows
    up); check it reproduces a near-1 PowerLaw within the shared range."""
    lo, hi = 0.5, 50.0
    log_map = injection.LogMapping(lo, hi)
    pl = injection.PowerLawMapping(0.999, 0.0, lo, hi)  # nu just below 1
    for r in np.linspace(0.05, 0.95, 19):
        assert log_map.Forward(r) == pytest.approx(pl.Forward(r), rel=2e-2)


def test_adaptive_mapping_converges_to_target():
    """The adaptive map drives its proposal toward a target density via the
    Accumulate/Refine hooks: the importance-weight variance collapses and the
    learned density tracks the (rising) target shape."""
    rng = np.random.default_rng(7)

    def f(x):  # target on [0,1]: integral 1, strongly peaked at x=1
        return 3.0 * x * x

    m = injection.AdaptiveMapping(0.0, 1.0, 24, 0.6, 1e-3)

    def coeff_of_variation(mapping, n=20000):
        x = np.array([mapping.Forward(r) for r in rng.random(n)])
        g = np.array([mapping.Density(xi) for xi in x])
        w = f(x) / g
        return float(w.std() / w.mean())

    cv_start = coeff_of_variation(m)          # uniform proposal
    for _ in range(8):
        for r in rng.random(4000):            # accumulate a batch ...
            x = m.Forward(r)
            g = m.Density(x)
            if g > 0:
                m.Accumulate(x, f(x) / g)
        m.Refine()                            # ... then refine once
    cv_end = coeff_of_variation(m)

    assert cv_end < 0.3 * cv_start, f"CV {cv_start:.3f} -> {cv_end:.3f}"
    # learned the increasing shape: density much higher near the peak than the tail
    assert m.Density(0.9) > 10.0 * m.Density(0.1)
