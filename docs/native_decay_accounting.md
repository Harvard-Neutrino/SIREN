# Native decay accounting

Decay selection, physical branching fractions and propagation lifetime are
separate inputs. These interfaces extend the [native phase-space API](native_phase_space.md).
Model amplitudes, widths and particle masses remain the responsibility of the
physics model; SIREN does not infer missing channels or their rates.

## Widths and channel selection

`PhaseSpaceDecay(signature, masses, partial_width, total_width, physical_channel)`
represents one decay signature. `physical_channel` is its normalized physical
final-state law, with an explicit compatible measure. An importance proposal
belongs on the injection process through `SetPhaseSpace`.

- `partial_width` is this signature's width in GeV. The model's
  `TotalDecayWidthAllFinalStates` returns only its owned partial width.
- `total_width` is the inclusive **parent** width, including omitted signatures.
  `ParentDecayWidth` exposes it. All positive declarations in one collection
  must agree; owned widths cannot exceed it.
- `InteractionCollection.TotalDecayWidthAllFinalStates` returns the parent
  width when explicitly declared. Legacy models can leave `ParentDecayWidth`
  at zero and retain the sum-of-owned-widths contract.
- `SetDecayChannels([signature, ...])` limits generated decays. It retains every
  model for lifetime and physical weighting. `None` restores all signatures;
  empty, duplicate and unavailable selections are rejected before replacement.

For a selected two-body decay with an isotropic physical law, generation is
conditioned on the selection. Its physical weight retains the partial/parent
branching fraction. This also holds for `Fixed()` vertices and parents at rest.
At-rest pure decays use scaled width ratios; material scattering requires
nonzero momentum. Negative or nonfinite selected and unselected widths fail.

The following executable example has a 20% generated branch, another 30%
modeled branch and a 50% unmodeled remainder. It checks accounting, without
representing a specific particle model, beam exposure or detector prediction.

```python
import math
import siren
from siren import dataclasses as dc, distributions as dist
from siren import injection as inj, interactions as xs, utilities
from siren.math import Vector3D

P = dc.Particle.ParticleType
parent_width = 1e-16  # GeV
signature = dc.InteractionSignature()
signature.primary_type = P.N4
signature.target_type = P.Decay
signature.secondary_types = [P.NuMu, P.Gamma]
selected = inj.PhaseSpaceDecay(
    signature, [0.0, 0.0], 0.2 * parent_width, parent_width,
    inj.Isotropic2BodyChannel())
other_signature = dc.InteractionSignature()
other_signature.primary_type = P.N4
other_signature.target_type = P.Decay
other_signature.secondary_types = [P.NuE, P.Gamma]
other = inj.PhaseSpaceDecay(
    other_signature, [0.0, 0.0], 0.3 * parent_width, parent_width,
    inj.Isotropic2BodyChannel())
collection = xs.InteractionCollection(P.N4, [selected, other])
collection.SetDecayChannels([signature])

detector = siren.detector.DetectorModel()
process = inj.PrimaryInjectionProcess(P.N4, collection)
process.weighting_mode = inj.VertexWeightingMode.Fixed()
process.distributions = [
    dist.PrimaryMass(1.0), dist.Monoenergetic(2.0),
    dist.PrimaryNeutrinoHelicityDistribution(),
    dist.FixedDirection(Vector3D(0, 0, 1)),
    dist.SphereVolumePositionDistribution(siren.geometry.Sphere(1.0, 0.0)),
]
physical = inj.PhysicalProcess(P.N4, collection)
physical.weighting_mode = process.weighting_mode
physical.distributions = process.distributions

generator = inj._Injector(100, detector, process, utilities.SIREN_random(19))
weighter = inj._Weighter([generator], detector, physical)
events = [generator.GenerateEvent() for _ in range(100)]
assert all(len(event.tree) == 1 for event in events)
assert all(event.tree[0].record.signature == signature for event in events)
assert math.isclose(sum(weighter.EventWeight(event) for event in events), 0.2,
                    rel_tol=1e-12)
```

## Survival from creation

`VertexWeightingMode.Propagated()` weights interaction and position within the
injection bounds, conditional on reaching their entry. Set
`VertexWeightingMode.PropagatedFromCreation()` on the **physical process** to
include upstream survival from `record.primary_initial_position` to that entry.
The injection process's survival flag does not change its proposal or enable
physical survival. Using the same preset on both processes is convenient, but
the physical process alone controls this factor. This applies to primary and
secondary processes. The native weight breakdown exposes `survival`.

The recorded initial position must still be the intended creation point:

- `PrimaryBoundedVertexDistribution` and `PrimaryPhysicalVertexDistribution`
  preserve a supplied creation point, including external-source coordinates.
- The secondary bounded/physical distributions preserve their parent's vertex.
- `FixedTargetPositionDistribution` generates and retains its own target source
  point. `PointSourcePositionDistribution` retains its configured origin only
  when clipping to the detector world has not moved that origin.
- Sphere/cylinder volume, column-depth, range and decay-range distributions
  replace the initial position with their injection entry. For those samplers,
  this mode gives unit upstream survival; it cannot recover an earlier external
  creation point. Do not use them to charge survival from that lost point.

The direct native API does not reject these sampler/mode combinations. Validate
this contract when constructing a process.

For decay length L, distance d to entry and path length ell inside the bounds,
the integrated physical probability for a branch b is
`b * exp(-d/L) * (1 - exp(-ell/L))`. Only the middle factor is new. A particle
created inside the bounds has unit upstream survival. Material interactions
also contribute to the upstream interaction depth. Do not multiply by another
survival factor in a downstream adapter.

The mode requires geometry bounds and both interaction and position factors.
Inconsistent modes are rejected. Negative/NaN depth, survival below the smallest
normal double, and overflow of the pooled inverse weight raise visible errors.
For example, exp(-700) is supported at ordinary branch fractions; combining it
with a very rare branch can exceed the representable inverse-weight range.
These guards do not cover every product: a product of otherwise representable
positive factors can still underflow to zero and be reported as zero physical
support. Such extreme products require additional numerical diagnostics.

The public `CrossSectionProbability`, `CrossSectionProbabilityWithPhaseSpace`,
`ChannelSelectionProbability` and `FixedVertexChannelSelectionProbability`
helpers describe physical, unconditioned probabilities. They do not apply the
generation-only decay selection; use the process weighter for that density.

## Persistence and compatibility

`InteractionCollection` version 1 stores selection; version 0 loads unrestricted.
`VertexWeightingMode` version 1 stores survival; version 0 loads with it off.
`PhaseSpaceDecay` begins at version 0. `Injector` and `PhysicalProcess` retain
U02's version 1. Native binary archives and pickle preserve channels, modes,
RNG continuation and weights. Live Python model state cannot be silently saved
as a native decay. The public Vertex/Simulation compiler and its configuration
persistence will be added in U04/U05; this slice uses direct native bindings.
