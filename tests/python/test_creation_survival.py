"""PropagatedFromCreation charges the survival from creation to a bounded volume.

A particle of known decay length is created on the beam axis at z0 and decays
inside a sphere whose near edge is at z = DISTANCE. Propagated() weights the
decay from that edge: weight x attempts = 1 - exp(-2 R / lambda) for this
single-row source. PropagatedFromCreation() also multiplies by the survival up
to the edge, exp(-(DISTANCE - z0) / lambda), computed here in closed form.
"""

import math
import pickle

import pytest

import siren
from siren import distributions as d

PT = siren.dataclasses.ParticleType
HBARC = siren.utilities.Constants.hbarc
MASS = 1.0
MOMENTUM = 3.0
ENERGY = math.hypot(MASS, MOMENTUM)
DISTANCE = 10.0  # m, from the origin to the near edge of the sphere
RADIUS = 1.0  # m


class FixedLengthDecay(siren.DecayModel):
    """N4 -> e- e+ whose lab decay length is 1 / inverse_length metres."""

    parent = PT.N4
    measure = siren.Measure.SolidAngleRest()

    def __init__(self, inverse_length):
        super().__init__()
        # lambda = (p / m) * hbarc / width
        self.width = inverse_length * MOMENTUM / MASS * HBARC
        self.daughters = (PT.EMinus, PT.EPlus)

    def total_width(self):
        return self.width

    def differential_width(self, record):
        return self.width / (4 * math.pi)

    def sample(self, record, random):
        self.sample_isotropic(record, random)


def native_decay():
    """A native N4 decay, so that archives may hold it (Python models may not)."""
    return siren.interactions.HNLDipoleDecay(
        MASS, [1e-7, 2e-7, 3e-7], siren.interactions.HNLDipoleDecay.Majorana)


def vertex(mode, inverse_length, z0=0.0, decay=None):
    source = d.PrimaryExternalDistribution(
        ["E", "px", "py", "pz", "x0", "y0", "z0", "m"],
        [[ENERGY, 0.0, 0.0, MOMENTUM, 0.0, 0.0, z0, MASS]])
    centre = siren.math.Vector3D(0, 0, DISTANCE + RADIUS)
    sphere = siren.geometry.Sphere(siren.geometry.Placement(centre), RADIUS, 0)
    decay = FixedLengthDecay(inverse_length) if decay is None else decay
    return siren.Vertex(PT.N4, [decay], distributions=[source],
                        position=siren.dist.BoundedPrimaryVertex(sphere, 100.0),
                        weighting=mode)


def simulation(mode, inverse_length, z0=0.0, events=16, seed=31, decay=None):
    return siren.Simulation(detector=siren.detector.DetectorModel(),
                            primary=vertex(mode, inverse_length, z0, decay),
                            events=events, seed=seed)


def run(mode, inverse_length, z0=0.0, **kwargs):
    return simulation(mode, inverse_length, z0, **kwargs).run(
        on_failure="raise", on_shortfall="raise")


# Survival from 0.61 to 1e-304. One minus an interaction probability would
# round everything from about 37 decay lengths on to zero.
DEPTHS = [0.5, 5.0, 50.0, 300.0, 700.0]


@pytest.mark.parametrize("depth", DEPTHS)
def test_weight_includes_closed_form_survival(depth):
    inverse_length = depth / DISTANCE
    plain = run(siren.Propagated(), inverse_length)
    full = run(siren.PropagatedFromCreation(), inverse_length)
    assert plain.attempts == full.attempts
    survival = math.exp(-depth)
    inside = -math.expm1(-2 * RADIUS * inverse_length)
    # SIREN forms beta*gamma from the four-momentum and the closed form from
    # p/m, so depths agree to a few ulp and survivals to depth * ulp.
    tolerance = 1e-14 * max(1.0, depth)
    for a, b, wa, wb in zip(plain.events, full.events, plain.weights, full.weights):
        assert list(a.tree[0].record.interaction_vertex) == list(b.tree[0].record.interaction_vertex)
        assert wa * plain.attempts == pytest.approx(inside, rel=1e-12, abs=0)
        assert wb / wa == pytest.approx(survival, rel=tolerance, abs=0)
        assert wb * full.attempts == pytest.approx(survival * inside, rel=tolerance, abs=0)


def test_creation_inside_the_volume_has_unit_survival():
    plain = run(siren.Propagated(), 2.0, z0=DISTANCE + RADIUS)
    full = run(siren.PropagatedFromCreation(), 2.0, z0=DISTANCE + RADIUS)
    assert list(full.weights) == list(plain.weights)


def test_breakdown_reports_survival():
    inverse_length = 3.0
    result = run(siren.PropagatedFromCreation(), inverse_length, events=4)
    plain = run(siren.Propagated(), inverse_length, events=4)
    for i in range(len(result)):
        line = result.explain(i).vertices[0]
        assert line.survival == pytest.approx(math.exp(-DISTANCE * inverse_length), rel=1e-13)
        assert plain.explain(i).vertices[0].survival == 1.0
        assert line.physical == pytest.approx(plain.explain(i).vertices[0].physical * line.survival,
                                              rel=1e-13)
    assert "survival=" in str(result.explain(0))


def test_binary_archives_keep_the_mode_and_weights(tmp_path):
    sim = simulation(siren.PropagatedFromCreation(), None, events=4, decay=native_decay())
    result = sim.run(on_failure="raise", on_shortfall="raise")
    engine = sim.weighter.engine
    assert engine.GetPrimaryPhysicalProcess().GetWeightingMode() == siren.PropagatedFromCreation()
    weights = [engine.EventWeight(tree) for tree in result.events]

    path = str(tmp_path / "creation")
    engine.SaveWeighter(path)
    loaded = siren.injection._Weighter([sim.injector.engine], path)
    assert loaded.GetPrimaryPhysicalProcess().GetWeightingMode() == siren.PropagatedFromCreation()
    assert [loaded.EventWeight(tree) for tree in result.events] == weights

    sim.injector.save(str(tmp_path / "injector"))
    restored = siren.injection.Injector.load(str(tmp_path / "injector"))
    assert restored.engine.GetPrimaryProcess().GetWeightingMode() == siren.PropagatedFromCreation()


def test_pickle_keeps_the_mode():
    # The external table above has no portable-binary pickle registration
    # (SIR-20); these native distributions do. Configuration round trip only.
    vertex = siren.Vertex(PT.N4, [native_decay()], distributions=[
        d.PrimaryMass(MASS), d.Monoenergetic(ENERGY),
        d.FixedDirection(siren.math.Vector3D(0, 0, 1)),
        d.PointSourcePositionDistribution(siren.math.Vector3D(), DISTANCE)],
        weighting=siren.PropagatedFromCreation())
    injector = siren.injection.Injector(detector=siren.detector.DetectorModel(),
                                        primary=vertex, events=1, seed=3)
    restored = pickle.loads(pickle.dumps(injector))
    mode = restored.engine.GetPrimaryProcess().GetWeightingMode()
    assert mode == siren.PropagatedFromCreation() and mode.survival_from_creation


@pytest.mark.parametrize("preset", [siren.Fixed, siren.ExternalBounds])
def test_survival_requires_geometry_propagation(preset):
    mode = preset()
    mode.survival_from_creation = True
    with pytest.raises(RuntimeError, match="survival_from_creation"):
        mode.Validate()
    process = siren.injection.PhysicalProcess()
    with pytest.raises(RuntimeError, match="survival_from_creation"):
        process.weighting_mode = mode
    assert process.weighting_mode == siren.Propagated()


def test_unrepresentable_survival_raises():
    # exp(-800) is below the smallest double.
    with pytest.raises(siren.errors.WeightCalculationError, match="CreationSurvival"):
        run(siren.PropagatedFromCreation(), 80.0, events=1)


def test_pooled_injectors_must_agree():
    a = simulation(siren.Propagated(), 1.0, events=2, seed=1)
    b = simulation(siren.PropagatedFromCreation(), 1.0, events=2, seed=2)
    weighter = siren.injection.Weighter(a.injector, b.injector)
    tree = a.run(on_failure="raise", on_shortfall="raise").events[0]
    with pytest.raises(siren.errors.ConfigurationError, match="survival_from_creation"):
        weighter.event_weight(tree)
    same = siren.injection.Weighter(a.injector, a.injector)
    assert same.event_weight(tree) > 0


def test_results_with_different_modes_do_not_merge():
    plain = run(siren.Propagated(), 1.0, events=2)
    full = run(siren.PropagatedFromCreation(), 1.0, events=2)
    with pytest.raises(siren.errors.ConfigurationError):
        siren.Results.merge([plain, full])
