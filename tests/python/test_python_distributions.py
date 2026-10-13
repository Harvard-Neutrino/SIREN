"""Python-subclassed distributions driven through the C++ injector and weighter.

Unit tests call bound base-class methods so the C++ virtual dispatch crosses
back into the python overrides; end-to-end tests run python distributions
through the Injector and Weighter against the CCM detector model with a
dummy cross section. The distribution classes live at module scope so that
pickle can locate them by qualified name.
"""
import os
import math as pymath

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

TIME_CONST = 41.25
BOX_HALF_WIDTH = 0.5
BOX_CENTER = (0.0, 0.0, 1.0)


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


def _interaction_collection():
    return interactions.InteractionCollection(NuMu, [interactions.DummyCrossSection()])


class UniformEnergy(distributions.PrimaryEnergyDistribution):
    """Uniform energy on [emin, emax] via the SampleEnergy hook."""

    def __init__(self, emin, emax):
        super().__init__()
        self.emin = emin
        self.emax = emax

    def SampleEnergy(self, rand, detector_model, interactions, record):
        return rand.Uniform(self.emin, self.emax)

    def GenerationProbability(self, detector_model, interactions, record):
        energy = record.primary_momentum[0]
        if self.emin <= energy <= self.emax:
            return 1.0 / (self.emax - self.emin)
        return 0.0


class NamedEnergy(UniformEnergy):
    """Uniform energy with an explicit Name override."""

    def Name(self):
        return "not-the-class-name"


class CloneableEnergy(UniformEnergy):
    """Uniform energy with an explicit clone override."""

    def clone(self):
        return CloneableEnergy(self.emin, self.emax)


class FixedZDirection(distributions.PrimaryDirectionDistribution):
    """Always points along +z via the SampleDirection hook."""

    def __init__(self):
        super().__init__()

    def SampleDirection(self, rand, detector_model, interactions, record):
        return smath.Vector3D(0.0, 0.0, 1.0)

    def GenerationProbability(self, detector_model, interactions, record):
        return 1.0


class PairEnergyDirection(distributions.PrimaryEnergyDirectionDistribution):
    """Fixed (energy, direction) pair via the SampleEnergyAndDirection hook."""

    def __init__(self, energy):
        super().__init__()
        self.fixed_energy = energy

    def SampleEnergyAndDirection(self, rand, detector_model, interactions, record):
        return (self.fixed_energy, smath.Vector3D(0.0, 0.0, 1.0))

    def GenerationProbability(self, detector_model, interactions, record):
        return 1.0


class BoxVertex(distributions.VertexPositionDistribution):
    """Uniform vertex position inside an axis-aligned box."""

    def __init__(self):
        super().__init__()

    def SamplePosition(self, rand, detector_model, interactions, record):
        cx, cy, cz = BOX_CENTER
        h = BOX_HALF_WIDTH
        x = rand.Uniform(cx - h, cx + h)
        y = rand.Uniform(cy - h, cy + h)
        z = rand.Uniform(cz - h, cz + h)
        initial = smath.Vector3D(x, y, cz - 20.0)
        vertex = smath.Vector3D(x, y, z)
        return (initial, vertex)

    def GenerationProbability(self, detector_model, interactions, record):
        cx, cy, cz = BOX_CENTER
        h = BOX_HALF_WIDTH
        x, y, z = record.interaction_vertex
        inside = (abs(x - cx) <= h) and (abs(y - cy) <= h) and (abs(z - cz) <= h)
        if not inside:
            return 0.0
        return 1.0 / (2.0 * h) ** 3

    def InjectionBounds(self, detector_model, interactions, record):
        x, y, z = record.interaction_vertex
        cz = BOX_CENTER[2]
        h = BOX_HALF_WIDTH
        first = smath.Vector3D(x, y, cz - h)
        last = smath.Vector3D(x, y, cz + h)
        return (first, last)


class InitialTimeShift(distributions.PrimaryInjectionDistribution):
    """Generic injection distribution that stamps a fixed initial time."""

    def __init__(self):
        super().__init__()

    def Sample(self, rand, detector_model, interactions, record):
        record.initial_time = TIME_CONST

    def GenerationProbability(self, detector_model, interactions, record):
        return 1.0

    def SetVariables(self):
        # The declaration vocabulary has no member for the initial time;
        # an empty set is the honest declaration for a time-only stamp.
        # Subclassing the generic base makes SetVariables mandatory (it is
        # pure there); the typed bases inherit their C++ declarations.
        return set()


class UnitWeightable(distributions.WeightableDistribution):
    """Physical-side distribution contributing a constant factor."""

    def __init__(self, factor=1.0):
        super().__init__()
        self.factor = factor

    def GenerationProbability(self, detector_model, interactions, record):
        return self.factor


class LengthSettingSecondaryVertex(distributions.SecondaryVertexPositionDistribution):
    """Secondary vertex distribution that pins the decay length."""

    def __init__(self, length):
        super().__init__()
        self.length = length
        self.sample_calls = 0

    def SampleVertex(self, rand, detector_model, interactions, record):
        self.sample_calls += 1
        record.length = self.length

    def GenerationProbability(self, detector_model, interactions, record):
        return 1.0

    def InjectionBounds(self, detector_model, interactions, record):
        x, y, z = record.interaction_vertex
        return (smath.Vector3D(x, y, z - 1.0), smath.Vector3D(x, y, z + 1.0))


class FlagAreaDistribution(distributions.PrimaryAreaDistribution):
    """Area distribution that records that it was called."""

    def __init__(self):
        super().__init__()
        self.sample_calls = 0

    def SamplePointOfClosestApproach(self, rand, detector_model, interactions, record):
        self.sample_calls += 1
        return smath.Vector3D(0.0, 0.0, 0.0)

    def GenerationProbability(self, detector_model, interactions, record):
        return 1.0


class FlagSecondaryDistribution(distributions.SecondaryInjectionDistribution):
    """Generic secondary distribution that records that it was called."""

    def __init__(self):
        super().__init__()
        self.sample_calls = 0

    def Sample(self, rand, detector_model, interactions, record):
        self.sample_calls += 1

    def GenerationProbability(self, detector_model, interactions, record):
        return 1.0


def _make_injector(dm, energy_dist=None, direction_dist=None, position_dist=None,
                   extra_dists=None, seed=7, n_inject=1000):
    dists = [distributions.PrimaryMass(0)]
    dists.append(energy_dist if energy_dist is not None
                 else distributions.PowerLaw(2.0, 0.5, 5.0))
    dists.append(direction_dist if direction_dist is not None
                 else distributions.IsotropicDirection())
    dists.append(position_dist if position_dist is not None
                 else distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0))
    if extra_dists:
        dists.extend(extra_dists)

    primary_inj = injection.PrimaryInjectionProcess()
    primary_inj.primary_type = NuMu
    primary_inj.interactions = _interaction_collection()
    primary_inj.distributions = dists

    rand = utilities.SIREN_random(seed)
    return injection._Injector(n_inject, dm, primary_inj, rand)


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


# --- unit-level dispatch through the bound base-class methods ---


def test_energy_distribution_dispatch():
    dist = UniformEnergy(1.0, 3.0)
    rand = utilities.SIREN_random(11)
    record = dc.PrimaryDistributionRecord(NuMu)
    # The base class Sample is C++; it calls back into python SampleEnergy.
    dist.Sample(rand, None, None, record)
    assert 1.0 <= record.energy <= 3.0


def test_direction_distribution_dispatch():
    dist = FixedZDirection()
    rand = utilities.SIREN_random(11)
    record = dc.PrimaryDistributionRecord(NuMu)
    dist.Sample(rand, None, None, record)
    assert record.direction == pytest.approx((0.0, 0.0, 1.0))


def test_energy_direction_distribution_dispatch():
    dist = PairEnergyDirection(2.5)
    rand = utilities.SIREN_random(11)
    record = dc.PrimaryDistributionRecord(NuMu)
    dist.Sample(rand, None, None, record)
    assert record.energy == pytest.approx(2.5)
    assert record.direction == pytest.approx((0.0, 0.0, 1.0))


def test_position_distribution_dispatch():
    dist = BoxVertex()
    rand = utilities.SIREN_random(11)
    record = dc.PrimaryDistributionRecord(NuMu)
    dist.Sample(rand, None, None, record)
    x, y, z = record.interaction_vertex
    assert abs(x - BOX_CENTER[0]) <= BOX_HALF_WIDTH
    assert abs(y - BOX_CENTER[1]) <= BOX_HALF_WIDTH
    assert abs(z - BOX_CENTER[2]) <= BOX_HALF_WIDTH

    ir = dc.InteractionRecord()
    ir.interaction_vertex = [x, y, z]
    first, last = dist.InjectionBounds(None, None, ir)
    assert first.GetZ() == pytest.approx(BOX_CENTER[2] - BOX_HALF_WIDTH)
    assert last.GetZ() == pytest.approx(BOX_CENTER[2] + BOX_HALF_WIDTH)


def test_area_distribution_dispatch():
    dist = FlagAreaDistribution()
    rand = utilities.SIREN_random(11)
    record = dc.PrimaryDistributionRecord(NuMu)
    dist.Sample(rand, None, None, record)
    assert dist.sample_calls == 1


def test_secondary_distribution_dispatch():
    parent = dc.InteractionRecord()
    parent.signature.primary_type = NuMu
    parent.signature.secondary_types = [NuMu]
    parent.secondary_momenta = [[1.0, 0.0, 0.0, 0.5]]
    parent.secondary_masses = [0.0]
    parent.secondary_helicities = [-0.5]
    parent.secondary_ids = [dc.ParticleID(1, 1)]

    rand = utilities.SIREN_random(11)

    generic = FlagSecondaryDistribution()
    record = dc.SecondaryDistributionRecord(parent, 0)
    generic.Sample(rand, None, None, record)
    assert generic.sample_calls == 1

    vertex = LengthSettingSecondaryVertex(12.5)
    record = dc.SecondaryDistributionRecord(parent, 0)
    # The base class Sample is C++; it calls back into python SampleVertex.
    vertex.Sample(rand, None, None, record)
    assert vertex.sample_calls == 1
    assert record.length == pytest.approx(12.5)


def test_weightable_distribution_dispatch():
    dist = UnitWeightable(0.75)
    ir = dc.InteractionRecord()
    assert dist.GenerationProbability(None, None, ir) == pytest.approx(0.75)


# --- trampoline defaults ---


def test_name_defaults_to_class_name():
    assert UniformEnergy(1.0, 3.0).Name() == "UniformEnergy"
    assert UnitWeightable().Name() == "UnitWeightable"


def test_name_override_wins():
    assert NamedEnergy(1.0, 3.0).Name() == "not-the-class-name"


def test_equality_defaults_to_identity():
    d1 = UniformEnergy(1.0, 3.0)
    d2 = UniformEnergy(1.0, 3.0)
    assert d1 == d1
    assert not (d1 == d2)


def test_clone_requires_override():
    with pytest.raises(RuntimeError, match="clone"):
        UniformEnergy(1.0, 3.0).clone()

    cloned = CloneableEnergy(1.0, 3.0).clone()
    assert isinstance(cloned, CloneableEnergy)
    assert cloned.emin == pytest.approx(1.0)
    assert cloned.emax == pytest.approx(3.0)


def test_pickle_roundtrip_preserves_dispatch():
    import pickle

    dist = UniformEnergy(1.5, 2.5)
    restored = pickle.loads(pickle.dumps(dist))
    assert isinstance(restored, UniformEnergy)
    assert restored.emin == pytest.approx(1.5)
    assert restored.emax == pytest.approx(2.5)

    rand = utilities.SIREN_random(11)
    record = dc.PrimaryDistributionRecord(NuMu)
    restored.Sample(rand, None, None, record)
    assert 1.5 <= record.energy <= 2.5


# --- end-to-end through the Injector ---


def test_injector_with_python_energy_distribution(detector_model):
    inj = _make_injector(detector_model, energy_dist=UniformEnergy(1.0, 3.0), seed=21)
    events = _generate(inj, 30)
    assert events
    energies = [ev.tree[0].record.primary_momentum[0] for ev in events]
    assert all(1.0 <= e <= 3.0 for e in energies)
    assert max(energies) - min(energies) > 1e-6


def test_injector_with_python_direction_distribution(detector_model):
    inj = _make_injector(detector_model, direction_dist=FixedZDirection(), seed=22)
    events = _generate(inj, 10)
    assert events
    for ev in events:
        e, px, py, pz = ev.tree[0].record.primary_momentum
        p = pymath.sqrt(px * px + py * py + pz * pz)
        assert pz / p == pytest.approx(1.0)


def test_injector_with_python_position_distribution(detector_model):
    # The injector locates the vertex distribution with a dynamic cast, so
    # this also verifies that a python subclass is recognized as a
    # VertexPositionDistribution.
    inj = _make_injector(detector_model, position_dist=BoxVertex(), seed=23)
    events = _generate(inj, 10)
    assert events
    for ev in events:
        x, y, z = ev.tree[0].record.interaction_vertex
        assert abs(x - BOX_CENTER[0]) <= BOX_HALF_WIDTH
        assert abs(y - BOX_CENTER[1]) <= BOX_HALF_WIDTH
        assert abs(z - BOX_CENTER[2]) <= BOX_HALF_WIDTH


def test_injector_with_python_generic_distribution(detector_model):
    inj = _make_injector(detector_model, extra_dists=[InitialTimeShift()], seed=24)
    events = _generate(inj, 10)
    assert events
    for ev in events:
        assert ev.tree[0].record.primary_initial_time == pytest.approx(TIME_CONST)


# --- end-to-end through the Weighter ---


def _physical_process(shared_energy_dist, extra_phys=None):
    primary_phys = injection.PhysicalProcess()
    primary_phys.primary_type = NuMu
    primary_phys.interactions = _interaction_collection()
    dists = [
        distributions.PrimaryMass(0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        shared_energy_dist,
    ]
    if extra_phys:
        dists.extend(extra_phys)
    primary_phys.distributions = dists
    return primary_phys


def test_weighter_with_python_distributions(detector_model):
    energy_dist = UniformEnergy(1.0, 3.0)

    dists = [
        distributions.PrimaryMass(0),
        energy_dist,
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
    ]
    primary_inj = injection.PrimaryInjectionProcess()
    primary_inj.primary_type = NuMu
    primary_inj.interactions = _interaction_collection()
    primary_inj.distributions = dists

    rand = utilities.SIREN_random(25)
    inj = injection._Injector(1000, detector_model, primary_inj, rand)

    events = _generate(inj, 20)
    assert events

    # The python energy distribution instance is shared between the injection
    # and physical processes, so it cancels in the weight by pointer identity.
    # The unshared python weightable contributes its constant factor.
    phys_shared = _physical_process(energy_dist, extra_phys=[UnitWeightable(1.0)])
    weighter = injection._Weighter([inj], detector_model, phys_shared)
    weights = [weighter.EventWeight(ev) for ev in events]
    assert all(pymath.isfinite(w) and w > 0.0 for w in weights)

    # An identical but unshared python energy distribution must give the same
    # weights, now computed by explicitly evaluating both sides in python.
    phys_unshared = _physical_process(UniformEnergy(1.0, 3.0), extra_phys=[UnitWeightable(1.0)])
    weighter_unshared = injection._Weighter([inj], detector_model, phys_unshared)
    weights_unshared = [weighter_unshared.EventWeight(ev) for ev in events]
    assert weights_unshared == pytest.approx(weights)

    # Scaling the physical side by a constant factor scales the weights.
    phys_scaled = _physical_process(energy_dist, extra_phys=[UnitWeightable(2.0)])
    weighter_scaled = injection._Weighter([inj], detector_model, phys_scaled)
    weights_scaled = [weighter_scaled.EventWeight(ev) for ev in events]
    assert weights_scaled == pytest.approx([2.0 * w for w in weights])


def test_injector_with_python_secondary_distribution(detector_model):
    # DummyCrossSection signatures produce the primary type among the
    # secondaries, so a secondary process on NuMu re-enters the chain once.
    energy_dist = UniformEnergy(1.0, 3.0)
    dists = [
        distributions.PrimaryMass(0),
        energy_dist,
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
    ]
    primary_inj = injection.PrimaryInjectionProcess()
    primary_inj.primary_type = NuMu
    primary_inj.interactions = _interaction_collection()
    primary_inj.distributions = dists

    secondary_vertex = LengthSettingSecondaryVertex(3.0)
    secondary_proc = injection.SecondaryInjectionProcess()
    secondary_proc.secondary_type = NuMu
    secondary_proc.interactions = _interaction_collection()
    secondary_proc.distributions = [secondary_vertex]

    rand = utilities.SIREN_random(26)
    inj = injection._Injector(1000, detector_model, primary_inj, [secondary_proc], rand)

    def stop_after_first_generation(tree, datum, index):
        return len(tree.tree) >= 2

    inj.SetStoppingCondition(stop_after_first_generation)

    events = _generate(inj, 5)
    assert events
    assert secondary_vertex.sample_calls > 0
    found_secondary = any(len(ev.tree) > 1 for ev in events)
    assert found_secondary
