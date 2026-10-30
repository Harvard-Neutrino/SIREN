# Directed phase-space proposals

`Directed` replaces the final-state proposal at a `Vertex`. It does not change
its physical decay width or cross section. Choose a daughter by particle name
(or by secondary index when the same type occurs twice):

```python
proposal = siren.Directed("Gamma", toward=target, fraction=0.98)
vertex = siren.Vertex("N4", model, distributions=distributions,
                      kinematics=proposal, weighting=siren.PropagatedFromCreation())
```

The remaining 0.02 samples the vertex's physical model. Keep a nonzero physical
fallback unless the directed proposal covers the whole physical support.
`fraction=1` explicitly restricts support. Missing physical models and ambiguous
daughter names fail at compilation.

For three-body final states, `channels.toward_3body` makes the factorization
explicit: `strategy="direct"` directs the named daughter and samples the other
two as a pair; `strategy="recursive", spectator="NuLight"` first produces that
spectator and a pair, then directs one member of the pair. `PairMass.uniform()`,
`breit_wigner(mass, width)`, `power_law(index, offset=...)`, and `tabulated(nodes,
cdf)` specify fixed pair-mass-squared proposals. They are proposal densities,
not physical resonance widths or branching factors. The same native channel
supports a stationary target in a three-body scattering final state.

For two-body scattering, `channels.scatter_toward(target, daughter=...)` selects
the outgoing particle. Its Q2, Bjorken-y and recoil-y conventions are explicit;
propagator and tabulated Q2 modes sample the transfer distribution rather than
restricting directions to the target. Unequal-mass processes can have physical
negative Q2. The propagator proposal excludes its pole and therefore generally
requires a physical fallback.

`channels.tile(target, n, daughter=...)` divides a Box into a fixed grid.
`directed_tiling` also provides angular sectors and Boolean union/subtraction
tiles. Volume weights require exact volumes; composite geometries need an
explicit exact-volume callback. Unresolved volume estimates fail rather than
using the bounding box's volume. Containment deduplication never removes a
volume based on sampled points. No adaptive tuning is included here.

Native directed channels and compiled native mixtures support pickle and C++
archives, including their geometry, roles, fixed mapping and mixture weights.
Old archives with obsolete branch-density semantics are rejected explicitly.
A Python physical callback still cannot be silently persisted as native state;
the existing authoring persistence checks apply.

Momenta use the shared `FourMomentum` and `RecordKinematics` interfaces.
Decay frames retain declared masses and accept primary energy rounding only
within the existing relative 2e-5 shell tolerance. Input record components are
preserved. A sampler's exhausted rejection loop is `SamplingFailure`, distinct
from a kinematically forbidden event; generation must keep those failed attempts
in its accounting.
