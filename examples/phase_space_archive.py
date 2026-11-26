"""Sample and weight a fixed-vertex dipole decay, then reload its archives."""

import math
from tempfile import TemporaryDirectory

import siren
from siren import dataclasses as dc, distributions as dist
from siren import injection as inj, interactions as xs, utilities
from siren.math import Vector3D

N4 = dc.Particle.ParticleType.N4
hnl_mass = 1.0          # GeV
dipole_coupling = 1e-6  # GeV^-1
model = xs.HNLDipoleDecay(
    hnl_mass, dipole_coupling, xs.HNLDipoleDecay.ChiralNature.Majorana)
interactions = xs.InteractionCollection(N4, [model])
detector = siren.detector.DetectorModel()

process = inj.PrimaryInjectionProcess()
process.primary_type = N4
process.interactions = interactions
process.weighting_mode = inj.VertexWeightingMode.Fixed()
process.distributions = [
    dist.PrimaryMass(hnl_mass),
    dist.Monoenergetic(2.0),  # GeV
    dist.PrimaryNeutrinoHelicityDistribution(),
    dist.FixedDirection(Vector3D(0, 0, 1)),
    dist.SphereVolumePositionDistribution(siren.geometry.Sphere(1.0, 0.0)),
]
for signature in model.GetPossibleSignaturesFromParent(N4):
    process.SetPhaseSpace(signature, inj.MultiChannelPhaseSpace([
        inj.PhysicalDecayChannel(model, signature),
        inj.Isotropic2BodyChannel(),
    ], [0.3, 0.7]))

physical = inj.PhysicalProcess()
physical.primary_type = N4
physical.interactions = interactions
physical.weighting_mode = inj.VertexWeightingMode.Fixed()
physical.distributions = process.distributions

n_events = 1000
generator = inj._Injector(n_events, detector, process, utilities.SIREN_random(19))
weighter = inj._Weighter([generator], detector, physical)
events = [generator.GenerateEvent() for _ in range(n_events)]
assert all(len(event.tree) == 1 for event in events)

# A Majorana dipole decay is isotropic, so the physical and mixture densities
# are equal and every event weighs exactly 1/n_events.
weights = [weighter.EventWeight(event) for event in events]
assert math.isclose(sum(weights), 1.0, rel_tol=1e-12)

with TemporaryDirectory() as directory:
    injector_path = f"{directory}/decay.injector"
    weighter_path = f"{directory}/decay"
    generator.SaveInjector(injector_path)
    weighter.SaveWeighter(weighter_path)  # writes decay.siren_weighter

    # The archive stores the RNG state, so the restored injector continues the
    # original stream. The engine passed here is used only for archives that
    # predate stored RNG state.
    restored = inj._Injector(n_events, injector_path, utilities.SIREN_random(0))
    restored_weighter = inj._Weighter([], weighter_path)
    assert restored.InjectionAttempts() == n_events
    assert restored_weighter.EventWeight(events[0]) == weights[0]
