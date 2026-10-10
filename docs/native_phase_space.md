# Native phase-space proposals

A process can sample an interaction's final state from a proposal instead of
the physical model. Register a `MultiChannelPhaseSpace` for an interaction
signature with `SetPhaseSpace`. The injector draws final states from that
mixture, and the weighter divides the physical density by the mixture density,
with both expressed in the same measure.

These classes live in `siren.injection`. The high-level `siren.Injector` and
`siren.Weighter` wrappers keep those names, so the native classes used below are
exposed as `siren.injection._Injector` and `siren.injection._Weighter`.

## Example

This example samples HNL dipole decays at a fixed vertex from a mixture of the
physical decay and an isotropic proposal, then weights and archives them. It
checks sampling and weighting only; there is no beam exposure or detector.

```python
import math
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

generator.SaveInjector("decay.injector")
weighter.SaveWeighter("decay")  # writes decay.siren_weighter

# The archive stores the RNG state, so the restored injector continues the
# original stream. The engine passed here is used only for archives that
# predate stored RNG state.
restored = inj._Injector(n_events, "decay.injector", utilities.SIREN_random(0))
restored_weighter = inj._Weighter([], "decay")
assert restored.InjectionAttempts() == n_events
assert restored_weighter.EventWeight(events[0]) == weights[0]
```

## Weighting modes

`VertexWeightingMode.Fixed()` omits the transport and vertex-position factors
but keeps the probability of selecting the interaction channel. In the example,
the injection and physical processes share their source and position
distributions, so those factors cancel. The default, `Propagated()`, includes
the interaction and position probabilities over the injection bounds of the
vertex distribution. A moving particle that scatters needs material in the
detector model.

## Densities and measures

A channel provides `Sample`, `Density`, `Topology` and `Measure`. `Density` must
be defined at points drawn by any other channel in the same mixture. A mixture
draws channel *i* with probability `weights[i]` (normalized to sum to one), and
its density at a point is the weighted sum over **all** channels, not just the
one that drew it. Nested mixtures and the adapters for existing `Decay` and
`CrossSection` models follow the same rule.

Every density is differential in a declared measure. For example,
`CosThetaRest` is a density per rest-frame `cos(theta)` with a uniform azimuth
integrated out; converting it to `SolidAngleRest` divides by `2*pi`. The reverse
would require integrating an arbitrary density over azimuth, so it is rejected.

Interaction models declare their measure through `Measure()`, or through
`MeasureForSignature()` when channels differ. The default is `Unspecified`: an
undeclared density can be mixed only with channels in the same undeclared
measure, and any conversion raises an error. `HNLDipoleDecay`,
`ElectroweakDecay` and the two-body channels of `HNLDecay` declare
`CosThetaRest`. The other built-in interactions do not declare a measure yet.

The conversion functions also cover two-body scattering measures and the
three-body charts. An incompatible topology, missing kinematic inputs or an
unsupported conversion raises a typed error. Passing a convention to a physical
adapter declares the measure of that model's density; it does not convert a
density that is normalized in a different measure.

## Failed attempts and normalization

Each `GenerateEvent()` call is one attempt. A retryable `InjectionFailure`
returns an empty tree and records the reason, depth and particle in
`GetFailureLedger()`; the failed draw is not replaced. Configuration and measure
errors propagate. Check for empty events and recorded failures before using a
run.

Weights are normalized by the number of **attempts**, including failed ones.
Before generation starts, the configured number of events is used instead.
Generate the whole batch before weighting it, because weights change while the
attempt count grows. An empty tree cannot receive a positive weight.

When several injectors are pooled in one weighter, an injector whose generation
density is zero at an event contributes nothing to that event's weight; an
error is raised only if no injector could have produced the event. A zero
physical density gives weight zero. Negative or nonfinite densities raise
errors, and `EventWeightWithBreakdown` reports them as a NaN total with
per-vertex flags.

Pure decays do not look up material. For a parent at rest, the decay channel is
chosen from partial-width ratios. Material interactions require a finite,
nonzero momentum.

## Archives

Native channels, mixtures, process phase-space maps, RNG state, attempt and
failure counters, and weighters can be saved to archives and pickled. If loading
a file fails, the receiving object keeps its previous configuration. Archives
written before RNG state was stored remain readable; loading one uses the
engine passed by the caller. The failure ledger's example messages and the last
partially generated tree are not saved.

An injector with a stopping-condition callback cannot be saved or pickled,
because the callback cannot be restored; this includes callbacks set through
the `siren.Injector` wrapper. Python-defined interactions and distributions keep
their pickled Python state, and state that cannot be pickled raises an error.
Do not load pickles or native archives from untrusted sources.

## Tests

Run the focused checks against a build whose Python package and native
libraries match the source:

```bash
python -m pytest tests/python/test_native_phase_space.py -q
ctest --test-dir build -R 'PhaseSpace|TwoBodyKinematics|Injector|Weighter' --output-on-failure
```
