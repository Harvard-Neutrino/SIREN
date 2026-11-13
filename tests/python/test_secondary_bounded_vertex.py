"""SecondaryBoundedVertexDistribution with a fiducial volume.

A ray that crosses the fiducial box is bounded by the crossing; a ray that
misses it falls back to the detector's outer bounds. For a decay the vertex
density along the crossing is exponential in the decay length.
"""

import math

import pytest

pytest.importorskip("siren")

from siren import dataclasses as dc
from siren import detector, distributions, geometry, interactions, utilities
from siren import math as smath

N4 = dc.ParticleType.N4
MASS = 0.1  # GeV
MOMENTUM = 1.0  # GeV
ENERGY = math.hypot(MASS, MOMENTUM)
WIDTH = 1e-16  # GeV
WORLD_RADIUS = 100.0  # m
BOX_CENTRE_Z = 50.0  # m
BOX_WIDTH = 2.0  # m


class FixedWidthDecay(interactions.Decay):
    """N4 -> nu gamma with a constant width; only the width is used."""

    def __init__(self):
        super().__init__()
        signature = dc.InteractionSignature()
        signature.primary_type = N4
        signature.target_type = dc.ParticleType.Decay
        signature.secondary_types = [dc.ParticleType.NuLight, dc.ParticleType.Gamma]
        self.signature = signature

    def equal(self, other):
        return self is other

    def TotalDecayWidthAllFinalStates(self, record):
        return WIDTH

    def TotalDecayWidth(self, arg):
        return WIDTH

    def DifferentialDecayWidth(self, record):
        return WIDTH

    def SampleFinalState(self, record, random):
        raise NotImplementedError("not sampled in these tests")

    def GetPossibleSignatures(self):
        return [self.signature]

    def GetPossibleSignaturesFromParent(self, primary):
        return [self.signature] if primary == N4 else []

    def FinalStateProbability(self, record):
        return 1.0

    def DensityVariables(self):
        return []


@pytest.fixture(scope="module")
def setup():
    model = detector.DetectorModel()
    sector = detector.DetectorSector()
    sector.name = "world"
    sector.level = 0
    sector.material_id = 0
    sector.geo = geometry.Sphere(geometry.Placement(smath.Vector3D(0, 0, 0)), WORLD_RADIUS, 0)
    sector.density = detector.ConstantDensityDistribution(0.0)
    model.AddSector(sector)
    box = geometry.Box(
        geometry.Placement(smath.Vector3D(0, 0, BOX_CENTRE_Z)), BOX_WIDTH, BOX_WIDTH, BOX_WIDTH
    )
    vertex = distributions.SecondaryBoundedVertexDistribution(box, 1000.0)
    collection = interactions.InteractionCollection(N4, [FixedWidthDecay()])
    return model, vertex, collection


def record(direction, vertex):
    rec = dc.InteractionRecord()
    signature = rec.signature
    signature.primary_type = N4
    signature.target_type = dc.ParticleType.Decay
    signature.secondary_types = [dc.ParticleType.NuLight, dc.ParticleType.Gamma]
    rec.signature = signature
    rec.primary_mass = MASS
    rec.primary_momentum = [ENERGY] + [MOMENTUM * c for c in direction]
    rec.primary_initial_position = [0.0, 0.0, 0.0]
    rec.interaction_vertex = vertex
    return rec


def xyz(v):
    return (v.GetX(), v.GetY(), v.GetZ())


def test_crossing_bounds_the_vertex(setup):
    model, vertex, collection = setup
    first, last = vertex.InjectionBounds(model, collection, record([0, 0, 1], [0, 0, BOX_CENTRE_Z]))
    assert xyz(first) == pytest.approx((0, 0, BOX_CENTRE_Z - BOX_WIDTH / 2), abs=1e-12)
    assert xyz(last) == pytest.approx((0, 0, BOX_CENTRE_Z + BOX_WIDTH / 2), abs=1e-12)


def test_miss_falls_back_to_outer_bounds(setup):
    model, vertex, collection = setup
    first, last = vertex.InjectionBounds(model, collection, record([1, 0, 0], [10, 0, 0]))
    assert xyz(first) == pytest.approx((0, 0, 0), abs=1e-12)
    assert xyz(last) == pytest.approx((WORLD_RADIUS, 0, 0), abs=1e-9)


def test_decay_density_along_crossing(setup):
    model, vertex, collection = setup
    depth_into_box = 0.3
    rec = record([0, 0, 1], [0, 0, BOX_CENTRE_Z - BOX_WIDTH / 2 + depth_into_box])
    decay_length = MOMENTUM / MASS * utilities.Constants.hbarc / WIDTH
    expected = (
        math.exp(-depth_into_box / decay_length)
        / decay_length
        / -math.expm1(-BOX_WIDTH / decay_length)
    )
    assert vertex.GenerationProbability(model, collection, rec) == pytest.approx(expected, rel=1e-12)
