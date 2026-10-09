# Native phase-space proposals

A process can now use a sampling distribution different from its physical
interaction. Assign a `MultiChannelPhaseSpace` to an interaction signature with
`SetPhaseSpace`. The injector samples that mixture; the weighter evaluates the
physical and proposal densities in a common measure.

This is the native API in `siren.injection`. The existing `siren.Injector` and
`siren.Weighter` wrappers remain available. This increment does not add a new
high-level simulation interface.

## A complete decay example

This conditional decay example uses the existing HNL dipole interaction. It
checks sampling and weighting, not a beam exposure or detector prediction.
All construction and execution calls use public bindings.

```python
import math
import siren
from siren import dataclasses as dc, distributions as dist
from siren import injection as inj, interactions as xs, utilities
from siren.math import Vector3D

N4 = dc.Particle.ParticleType.N4
model = xs.HNLDipoleDecay(1.0, 1e-6, xs.HNLDipoleDecay.ChiralNature.Majorana)
interactions = xs.InteractionCollection(N4, [model])
detector = siren.detector.DetectorModel()

process = inj.PrimaryInjectionProcess()
process.primary_type = N4
process.interactions = interactions
process.weighting_mode = inj.VertexWeightingMode.Fixed()
process.distributions = [
    dist.PrimaryMass(1.0),
    dist.Monoenergetic(2.0),
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

generator = inj._Injector(1000, detector, process, utilities.SIREN_random(19))
weighter = inj._Weighter([generator], detector, physical)
events = [generator.GenerateEvent() for _ in range(1000)]
assert all(len(event.tree) == 1 for event in events)
weights = [weighter.EventWeight(event) for event in events]
assert math.isclose(sum(weights), 1.0, rel_tol=1e-12)

generator.SaveInjector("decay.injector")
weighter.SaveWeighter("decay")
restored = inj._Injector(1000, "decay.injector", utilities.SIREN_random(0))
restored_weighter = inj._Weighter([], "decay")
assert restored.InjectionAttempts() == 1000
assert restored_weighter.EventWeight(events[0]) == weights[0]
```

`Fixed()` suppresses transport and vertex-position factors; it retains the
probability of selecting the interaction channel. Here the identical source
and spatial distributions cancel. The default `Propagated()` mode retains
interaction and position probabilities over the distribution's injection bounds.
A moving particle that scatters still requires a material geometry.

## Density contract

A channel supplies `Sample`, `Density`, `Topology`, and `Measure`. `Density`
must be evaluable at points drawn by other channels. A mixture samples channel
*i* with normalized probability `weights[i]`; its proposal density is the sum
of **all** weighted component densities at the resulting point. Nested mixtures
and adapters for existing `Decay` and `CrossSection` objects follow the same rule.

Densities carry their integration measure. For example, `CosThetaRest` means a
density per rest-frame `cos(theta)` with a declared uniform omitted azimuth.
Lifting it to `SolidAngleRest` divides by `2*pi`. An arbitrary solid-angle
density cannot be integrated over azimuth by a pointwise conversion, so that
reverse operation is rejected. HNL dipole decays declare `CosThetaRest`; their
existing `FinalStateProbability` values retain their historical meaning.

The native conversion functions also cover compatible two-body scattering
and three-body charts. Incompatible topology, missing kinematic inputs, or
an unsupported conversion produces a typed error. A convention override on a
physical adapter is a declaration by its caller; it does not repair a density
that is normalized in a different measure.

## Failure and normalization contract

`GenerateEvent()` consumes one attempt. A retryable `InjectionFailure` returns
an empty tree and records the reason, depth, and particle in `GetFailureLedger()`.
It does not replace the failed draw. Configuration and measure errors propagate.
Callers should inspect empty events and failures before accepting a run.

Weights use the number of **attempts**, including failed draws. Before generation,
the configured quota supplies the normalization. Generate the batch before
weighting it; weights change while the attempt count is increasing. Empty trees
cannot receive a positive event weight. Zero physical support has weight zero;
negative/nonfinite densities or nonpositive generation support raise errors.
`EventWeightWithBreakdown` reports invalid weights as NaN with diagnostic flags.

Pure decays do not navigate material. At rest, their discrete channel selection
uses partial-width ratios. Material interactions require finite nonzero momentum.
Selected versus inclusive parent-width accounting and creation-point survival
remain separate work; this PR does not introduce those interfaces.

## Persistence and review boundary

Native channels, mixtures, process maps, RNG state, attempt/failure counters,
and weighters support archives and Python pickle. A failed file load preserves
the receiving object's configuration. Old main injector archives remain readable;
they did not store the RNG, so loading those uses the caller's supplied engine.
Failure-ledger exemplars and the last partial tree are transient diagnostics.

Native archives reject live stopping callbacks rather than silently discard them.
That restriction also applies when an existing Python wrapper attempts to pickle
its underlying native injector. Serializable Python interaction/distribution
subclasses retain their trampoline state; non-pickleable state fails explicitly.
Do not load pickle or native archives from untrusted sources.

These are replacement-stack archives, not an interchange format for the unmerged
original stack: mixture tuning state is deliberately absent. Directed proposals,
tuning, beam readers, selected/parent decay widths, and model-authoring conveniences
will be separate increments. Dutta–Kim models belong in BSM-beam.

Run the focused acceptance checks against a matching staged runtime:

```bash
python -m pytest tests/python/test_native_phase_space.py -q
ctest --test-dir build -R 'PhaseSpace|TwoBodyKinematics|Injector|Weighter' --output-on-failure
```
