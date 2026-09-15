"""TrivialCrossSection driven through the injector and weighter.

The class is a pass-through cross section whose physics content is the total
cross section alone. Unit tests cover the tabulated total (interpolation,
threshold, clamping) and the signature conventions. The end-to-end tests
inject against the CCM detector model with a constant stand-in cross section
and then weight the same events under two physical hypotheses, checking that
the weight ratio is exactly the ratio of the physical total cross sections:
the mechanism that defers the physical total to an external generator.
"""
import math as pymath
import os

import pytest

siren = pytest.importorskip("siren")

from siren import dataclasses as dc
from siren import distributions
from siren import injection
from siren import interactions
from siren import math as smath
from siren import utilities
from siren import _util

NuMu = dc.Particle.ParticleType.NuMu
NuE = dc.Particle.ParticleType.NuE
Nucleon = dc.Particle.ParticleType.Nucleon
Ar40Nucleus = dc.Particle.ParticleType.Ar40Nucleus

SIGMA_CONST = 1.0e-38  # cm^2
ENERGY_MIN = 0.5
ENERGY_MAX = 5.0


def _table_sigma(energy):
    """Linear stand-in for a generator total: 0.7e-38 cm^2 per GeV."""
    return 0.7e-38 * energy


def _table_cross_section():
    # Two knots suffice: linear interpolation reproduces _table_sigma exactly
    # across the injected energy range.
    return interactions.TrivialCrossSection(
        {NuMu: ([0.1, 10.0], [_table_sigma(0.1), _table_sigma(10.0)])},
        [Nucleon],
    )


def _constant_cross_section():
    return interactions.TrivialCrossSection(SIGMA_CONST, [NuMu], [Nucleon])


def _detector_model():
    try:
        dm = siren.detector.DetectorModel()
        det_dir = _util.get_detector_model_path("CCM")
        dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
        dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
    except Exception as e:  # detector model files not available in this env
        pytest.skip(f"CCM detector model unavailable: {e}")
    return dm


@pytest.fixture(scope="module")
def detector_model():
    return _detector_model()


# --- unit-level behavior ---


def test_constant_cross_section_is_flat():
    xs = _constant_cross_section()
    assert xs.TotalCrossSection(NuMu, 0.5) == pytest.approx(SIGMA_CONST)
    assert xs.TotalCrossSection(NuMu, 50.0) == pytest.approx(SIGMA_CONST)
    assert xs.TotalCrossSection(NuE, 1.0) == 0.0


def test_tabulated_cross_section_interpolates():
    xs = _table_cross_section()
    for energy in (0.1, 0.5, 1.0, 2.5, 10.0):
        assert xs.TotalCrossSection(NuMu, energy) == pytest.approx(_table_sigma(energy))
    # Below the first knot the total vanishes; above the last it clamps.
    assert xs.TotalCrossSection(NuMu, 0.05) == 0.0
    assert xs.TotalCrossSection(NuMu, 100.0) == pytest.approx(_table_sigma(10.0))


def test_signatures_are_pass_through():
    xs = _constant_cross_section()
    signatures = xs.GetPossibleSignaturesFromParents(NuMu, Nucleon)
    assert len(signatures) == 1
    sig = signatures[0]
    assert sig.primary_type == NuMu
    assert sig.target_type == Nucleon
    assert list(sig.secondary_types) == [NuMu, Nucleon]
    assert xs.GetPossibleSignaturesFromParents(NuE, Nucleon) == []


def test_nucleus_target_convention():
    xs = interactions.TrivialCrossSection(SIGMA_CONST, [NuMu], [Ar40Nucleus])
    signatures = xs.GetPossibleSignaturesFromParents(NuMu, Ar40Nucleus)
    assert len(signatures) == 1
    assert signatures[0].target_type == Ar40Nucleus


def test_invalid_construction_raises():
    with pytest.raises(Exception):
        interactions.TrivialCrossSection(0.0, [NuMu], [Nucleon])
    with pytest.raises(Exception):
        interactions.TrivialCrossSection({NuMu: ([3.0, 1.0], [1e-38, 2e-38])}, [Nucleon])


# --- end-to-end: injection with a stand-in, physical totals deferred ---


def _make_injector(dm, standin, seed=31, n_inject=1000):
    dists = [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, ENERGY_MIN, ENERGY_MAX),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
    ]
    primary_inj = injection.PrimaryInjectionProcess()
    primary_inj.primary_type = NuMu
    primary_inj.interactions = interactions.InteractionCollection(NuMu, [standin])
    primary_inj.distributions = dists

    rand = utilities.SIREN_random(seed)
    return injection._Injector(n_inject, dm, primary_inj, rand), dists


def _physical_process(cross_section, dists):
    primary_phys = injection.PhysicalProcess()
    primary_phys.primary_type = NuMu
    primary_phys.interactions = interactions.InteractionCollection(NuMu, [cross_section])
    # Share the injection distribution instances so every distribution factor
    # cancels in the weight; only the interaction factors differ between
    # physical hypotheses.
    primary_phys.distributions = [d for d in dists
                                  if not isinstance(d, distributions.VertexPositionDistribution)]
    return primary_phys


def _generate(inj, n):
    events = []
    for _ in range(n):
        try:
            ev = inj.GenerateEvent()
        except RuntimeError:
            break
        if len(ev.tree) > 0:
            events.append(ev)
    return events


def test_final_state_is_pass_through(detector_model):
    standin = _constant_cross_section()
    inj, _ = _make_injector(detector_model, standin, seed=32)
    events = _generate(inj, 10)
    assert events
    for ev in events:
        record = ev.tree[0].record
        assert list(record.signature.secondary_types) == [NuMu, Nucleon]
        primary = record.primary_momentum
        neutrino = record.secondary_momenta[0]
        for a, b in zip(primary, neutrino):
            assert b == pytest.approx(a, rel=1e-12, abs=1e-15)


def test_weight_ratio_defers_total_cross_section(detector_model):
    standin = _constant_cross_section()
    inj, dists = _make_injector(detector_model, standin, seed=33)
    events = _generate(inj, 25)
    assert events

    # Physical hypothesis A: the same constant total used for generation.
    weighter_const = injection._Weighter(
        [inj], detector_model, _physical_process(standin, dists))
    weights_const = [weighter_const.EventWeight(ev) for ev in events]
    assert all(pymath.isfinite(w) and w > 0.0 for w in weights_const)

    # Physical hypothesis B: a tabulated (generator-provided) total.
    weighter_table = injection._Weighter(
        [inj], detector_model, _physical_process(_table_cross_section(), dists))
    weights_table = [weighter_table.EventWeight(ev) for ev in events]
    assert all(pymath.isfinite(w) and w > 0.0 for w in weights_table)

    # The weights differ exactly by the ratio of physical totals at the event
    # energy: the interaction probability is linear in the total cross section
    # for a thin target, the normalized position density and the final-state
    # probability are cross-section independent, and the generation side is
    # common to both weighters.
    for ev, w_const, w_table in zip(events, weights_const, weights_table):
        energy = ev.tree[0].record.primary_momentum[0]
        expected = _table_sigma(energy) / SIGMA_CONST
        assert w_table / w_const == pytest.approx(expected, rel=1e-6)


def test_collection_retrieval_preserves_type():
    # The collection hands back the concrete type, so the tabulated totals
    # remain reachable after registration in an InteractionCollection.
    standin = _constant_cross_section()
    collection = interactions.InteractionCollection(NuMu, [standin])
    assert collection.HasCrossSections()
    xs = collection.GetCrossSectionsForTarget(Nucleon)[0]
    assert isinstance(xs, interactions.TrivialCrossSection)
    assert xs.TotalCrossSection(NuMu, 1.0) == pytest.approx(SIGMA_CONST)
