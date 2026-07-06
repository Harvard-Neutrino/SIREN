"""Weighter.explain() over the engine breakdown, weight_all as ndarray, and
the removal of private-injector reach-throughs.

Data-free where possible; the assembled-chain tests skip cleanly when the CCM
detector data files are absent.
"""

import math
import os
import warnings

import numpy as np
import pytest

import siren
from siren import dataclasses as dc
from siren import injection
from siren import interactions
from siren import distributions
from siren import detector
from siren import math as smath
from siren import utilities
from siren import _util

_NuMu = dc.Particle.ParticleType.NuMu


def _skip_unless_ccm_data():
    try:
        det_dir = _util.get_detector_model_path("CCM")
    except ValueError as e:
        pytest.skip("CCM detector model path unavailable: {}".format(e))
    for name in ("materials.dat", "densities.dat"):
        if not os.path.exists(os.path.join(det_dir, name)):
            pytest.skip("CCM detector data missing: " + name)


def _load_ccm_detector():
    dm = detector.DetectorModel()
    det_dir = _util.get_detector_model_path("CCM")
    dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
    dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
    return dm


def _build_pair(dm, seed=4242):
    xs = interactions.DummyCrossSection()
    inj_dists = [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, 0.5, 5.0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
    ]
    inj = injection.Injector(
        number_of_events=50, detector_model=dm, seed=seed,
        primary_type=_NuMu, primary_interactions=[xs],
        primary_injection_distributions=inj_dists)

    weighter = injection.Weighter()
    weighter.injectors = [inj]
    weighter.detector_model = dm
    weighter.primary_type = _NuMu
    weighter.primary_interactions = [xs]
    weighter.primary_physical_distributions = [
        distributions.PrimaryMass(0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
    ]
    return inj, weighter


# ------------------------------------------------------------------ #
#  explain(): breakdown reconstructs the weight                        #
# ------------------------------------------------------------------ #

def test_explain_total_matches_scalar_weight():
    """explain(tree).total equals the scalar event weight to 1e-12."""
    _skip_unless_ccm_data()
    dm = _load_ccm_detector()
    inj, weighter = _build_pair(dm)

    trees = inj.generate(5, on_shortfall="ignore")
    if not trees:
        pytest.skip("no events generated for this data-free chain")
    for tree in trees:
        breakdown = weighter.explain(tree)
        scalar = weighter(tree)
        if math.isfinite(scalar) and math.isfinite(breakdown.total):
            assert abs(breakdown.total - scalar) <= 1e-12 * max(1.0, abs(scalar))
        else:
            assert math.isnan(breakdown.total) == math.isnan(scalar)


def test_explain_returns_weight_breakdown_with_vertices():
    """explain(tree) yields a report.WeightBreakdown with per-vertex lines."""
    from siren.report import WeightBreakdown
    _skip_unless_ccm_data()
    dm = _load_ccm_detector()
    inj, weighter = _build_pair(dm)
    trees = inj.generate(3, on_shortfall="ignore")
    if not trees:
        pytest.skip("no events generated")
    bd = weighter.explain(trees[0])
    assert isinstance(bd, WeightBreakdown)
    assert len(bd.vertices) >= 1
    assert "WeightBreakdown" in str(bd)


# ------------------------------------------------------------------ #
#  breakdown() is a deprecated alias                                  #
# ------------------------------------------------------------------ #

def test_breakdown_is_deprecated_alias():
    """breakdown(tree) emits a DeprecationWarning and matches explain()."""
    _skip_unless_ccm_data()
    dm = _load_ccm_detector()
    inj, weighter = _build_pair(dm)
    trees = inj.generate(2, on_shortfall="ignore")
    if not trees:
        pytest.skip("no events generated")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        bd = weighter.breakdown(trees[0])
    assert any(issubclass(x.category, DeprecationWarning) for x in w)
    assert bd.total == weighter.explain(trees[0]).total


# ------------------------------------------------------------------ #
#  weight_all returns a numpy array                                   #
# ------------------------------------------------------------------ #

def test_weight_all_returns_ndarray():
    """weight_all(trees) returns a numpy.ndarray matching per-event weights."""
    _skip_unless_ccm_data()
    dm = _load_ccm_detector()
    inj, weighter = _build_pair(dm)
    trees = inj.generate(5, on_shortfall="ignore")
    if not trees:
        pytest.skip("no events generated")
    batch = weighter.weight_all(trees)
    assert isinstance(batch, np.ndarray)
    individual = np.array([weighter(t) for t in trees])
    assert np.array_equal(batch, individual, equal_nan=True)


# ------------------------------------------------------------------ #
#  the private reach-throughs are gone                                #
# ------------------------------------------------------------------ #

def test_no_private_injector_reach_through_in_source():
    """Weighter.py reaches the raw injector via the public engine property."""
    src_path = os.path.join(os.path.dirname(_util.__file__), "Weighter.py")
    with open(src_path) as f:
        source = f.read()
    assert "_Injector__injector" not in source


# ------------------------------------------------------------------ #
#  spec-form constructor accepts injectors                            #
# ------------------------------------------------------------------ #

def test_spec_form_ctor_accepts_injectors():
    """Weighter(injector, primary_physical=[...]) builds a working weighter."""
    _skip_unless_ccm_data()
    dm = _load_ccm_detector()
    xs = interactions.DummyCrossSection()
    inj_dists = [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, 0.5, 5.0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
    ]
    inj = injection.Injector(
        number_of_events=10, detector_model=dm, seed=7,
        primary_type=_NuMu, primary_interactions=[xs],
        primary_injection_distributions=inj_dists)

    weighter = injection.Weighter(
        inj,
        primary_physical=[
            distributions.PrimaryMass(0),
            distributions.PrimaryNeutrinoHelicityDistribution(),
            distributions.IsotropicDirection(),
        ])
    trees = inj.generate(3, on_shortfall="ignore")
    if not trees:
        pytest.skip("no events generated")
    weights = weighter.weight_all(trees)
    assert isinstance(weights, np.ndarray)
    assert len(weights) == len(trees)
