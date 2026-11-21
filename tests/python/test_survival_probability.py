"""Pre-injection survival probabilities from the native Weighter.

``Weighter.GetSurvivalProbabilities`` (``survival_probabilities`` in Python)
gives the probability that a particle survives from its creation point to the
near edge of its injection region. For a decay of constant width the answer is
exp(-L / lambda), with lambda = beta*gamma*hbarc/width. These tests compare the
native value with that closed form, including survival probabilities far below
double-precision epsilon, where one minus an interaction probability loses all
of its digits, and check that an invalid decay length is rejected.
"""

import math

import pytest

pytest.importorskip("siren")

from siren import dataclasses as dc
from siren import detector, distributions, geometry, injection, interactions, utilities
from siren import math as smath

N4 = dc.ParticleType.N4
MASS = 0.1  # GeV
MOMENTUM = 1.0  # GeV, along +z
ENERGY = math.hypot(MASS, MOMENTUM)
RADIUS = 1.0  # m, fiducial sphere radius
DISTANCE = 50.0  # m, creation point to the near edge of the sphere
BETA_GAMMA = MOMENTUM / MASS


class FixedWidthDecay(interactions.Decay):
    """N4 -> nu gamma with a constant total width; only the width is used."""

    def __init__(self, width):
        super().__init__()
        self.width = width
        signature = dc.InteractionSignature()
        signature.primary_type = N4
        signature.target_type = dc.ParticleType.Decay
        signature.secondary_types = [dc.ParticleType.NuLight, dc.ParticleType.Gamma]
        self.signature = signature

    def equal(self, other):
        return self is other

    def TotalDecayWidthAllFinalStates(self, record):
        return self.width

    def TotalDecayWidth(self, arg):
        return self.width

    def DifferentialDecayWidth(self, record):
        return self.width

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


def width_for_depth(depth):
    """Width whose decay length puts ``depth`` decay lengths in DISTANCE."""
    return depth * BETA_GAMMA * utilities.Constants.hbarc / DISTANCE


def closed_form_survival(width):
    decay_length = BETA_GAMMA * utilities.Constants.hbarc / width
    return math.exp(-DISTANCE / decay_length)


def survival(width):
    """Native survival from the origin to the fiducial sphere's near edge."""
    decay = FixedWidthDecay(width)
    collection = interactions.InteractionCollection(N4, [decay])
    model = detector.DetectorModel()  # decay-only depth needs no material

    centre = smath.Vector3D(0, 0, DISTANCE + RADIUS)
    sphere = geometry.Sphere(geometry.Placement(centre), RADIUS, 0)
    injected = injection.PrimaryInjectionProcess()
    injected.primary_type = N4
    injected.interactions = collection
    injected.distributions = [
        distributions.PrimaryMass(MASS),
        distributions.Monoenergetic(ENERGY),
        distributions.FixedDirection(smath.Vector3D(0, 0, 1)),
        distributions.PrimaryBoundedVertexDistribution(sphere, 10 * DISTANCE),
    ]
    physical = injection.PhysicalProcess()
    physical.primary_type = N4
    physical.interactions = collection
    physical.distributions = [distributions.PrimaryMass(MASS)]

    injector = injection._Injector(1, model, injected, utilities.SIREN_random(1))
    weighter = injection._Weighter([injector], model, physical)

    record = dc.InteractionRecord()
    signature = record.signature
    signature.primary_type = N4
    signature.target_type = dc.ParticleType.Decay
    signature.secondary_types = [dc.ParticleType.NuLight, dc.ParticleType.Gamma]
    record.signature = signature
    record.primary_mass = MASS
    record.primary_momentum = [ENERGY, 0.0, 0.0, MOMENTUM]
    record.primary_initial_position = [0.0, 0.0, 0.0]
    record.interaction_vertex = [0.0, 0.0, DISTANCE + RADIUS]
    tree = dc.InteractionTree()
    tree.add_entry(record, None)

    # The injection region starts where the beam axis enters the sphere, at
    # z = DISTANCE; the closed-form comparison fails for any other entry.
    (value,) = weighter.GetSurvivalProbabilities(tree, 0)
    return value


# Depths from nearly transparent to a survival of 1e-304. Above about 37
# decay lengths, 1 - (1 - exp(-depth)) rounds to exactly zero.
DEPTHS = [1e-8, 1e-3, 0.5, 5.0, 20.0, 30.0, 40.0, 100.0, 300.0, 700.0]


@pytest.mark.parametrize("depth", DEPTHS)
def test_survival_matches_closed_form(depth):
    width = width_for_depth(depth)
    expected = closed_form_survival(width)
    # SIREN forms beta*gamma from the four-momentum, the closed form from
    # p/m, so the depths agree to a few ulp and the survival to depth*ulp.
    assert survival(width) == pytest.approx(expected, rel=1e-14 * max(1.0, depth), abs=0)


def test_stable_particle_survives():
    assert survival(0.0) == 1.0


@pytest.mark.parametrize("width", [-width_for_depth(1.0), float("nan")])
def test_invalid_depth_raises(width):
    with pytest.raises(RuntimeError, match="interaction depth"):
        survival(width)
