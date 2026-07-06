"""Legacy keyword and Vertex-spec entry points compile to one engine config.

The two authoring surfaces are a translation shim over a single build path, so
at a fixed seed they must produce identical event streams. These tests
parametrize a shared scenario over both surfaces and assert stream identity of
the generated trees and their weights.
"""

import os

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


def _distributions():
    return [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, 0.5, 5.0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(
            smath.Vector3D(0, 0, 0), 25.0),
    ]


SEED = 20240601
N = 25


def _build_legacy(dm):
    """The scenario spelled with legacy keyword arguments."""
    return injection.Injector(
        number_of_events=N,
        detector_model=dm,
        seed=SEED,
        primary_type=_NuMu,
        primary_interactions=[interactions.DummyCrossSection()],
        primary_injection_distributions=_distributions(),
    )


def _build_spec(dm):
    """The same scenario spelled with a Vertex spec."""
    v = siren.Vertex(
        _NuMu,
        interactions.DummyCrossSection(),
        distributions=_distributions(),
    )
    return injection.Injector(detector=dm, primary=v, events=N, seed=SEED)


def _stream(injector):
    """Realized (energy, vertex, depth) stream from a fixed-seed run."""
    trees = injector.generate(N, on_shortfall="raise")
    rows = []
    for t in trees:
        root = t.tree[0].record
        rows.append((
            tuple(root.primary_momentum),
            tuple(root.interaction_vertex),
        ))
    return rows


@pytest.fixture
def ccm_detector():
    _skip_unless_ccm_data()
    return _load_ccm_detector()


def test_legacy_and_spec_produce_identical_streams(ccm_detector):
    """Legacy kwargs and a Vertex spec yield byte-identical event streams."""
    legacy = _build_legacy(ccm_detector)
    spec = _build_spec(ccm_detector)

    legacy_stream = _stream(legacy)
    spec_stream = _stream(spec)

    assert len(legacy_stream) == len(spec_stream) == N
    assert legacy_stream == spec_stream


def test_legacy_and_spec_produce_identical_weights(ccm_detector):
    """The two surfaces weight their (identical) events identically."""
    legacy = _build_legacy(ccm_detector)
    spec = _build_spec(ccm_detector)

    legacy_trees = legacy.generate(N, on_shortfall="raise")
    spec_trees = spec.generate(N, on_shortfall="raise")

    phys = [
        distributions.PrimaryMass(0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
    ]

    def _weighter_for(inj):
        w = injection.Weighter()
        w.injectors = [inj]
        w.detector_model = ccm_detector
        w.primary_type = _NuMu
        w.primary_interactions = inj.primary_interactions
        w.primary_physical_distributions = list(phys)
        return w

    lw = _weighter_for(legacy).weight_all(legacy_trees)
    sw = _weighter_for(spec).weight_all(spec_trees)

    assert list(lw) == list(sw)


@pytest.mark.parametrize("build", [_build_legacy, _build_spec],
                         ids=["legacy", "spec"])
def test_both_surfaces_generate_requested_count(build, ccm_detector):
    """Each surface generates exactly the requested number of events."""
    inj = build(ccm_detector)
    trees = inj.generate(N, on_shortfall="raise")
    assert len(trees) == N
