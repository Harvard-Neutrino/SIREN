# Vertex configuration

`Vertex` declares one particle's physical models, proposal distributions,
phase-space channels and secondary expansion. `siren.injection.Injector`
compiles these declarations into the native processes.
`siren.injection.Weighter(generator)` uses the same declarations for the
physical target. Existing Injector keyword arguments use the same validation.

Run the complete, data-free example with the matching installed SIREN package:

```sh
python resources/examples/authoring_vertex.py
```

It defines an explicit two-body sampler, reads a fixed creation point and
momentum, samples a decay in a sphere, and checks the summed weight against
`exp(-1) * (1 - exp(-0.2))`. This is a transport check for a toy model;
it is not a beam-rate prediction. Model authoring and closure requirements
are described in [authoring.md](authoring.md).

## Physical laws and proposals

- `interactions` holds the models used for generation. `physical_interactions`
  optionally supplies different models for weighting; otherwise they are shared.
- `distributions` and `position` describe the proposal. `physical` contains
  physical flux factors. Supply each explicitly; the compiler never widens a
  physical energy interval to match the proposal.
- `kinematics` accepts `channels.physical()`, `channels.isotropic(daughter)` or
  a `channels.Channel` factory. Weighted sums form a normalized mixture:
  `0.3 * channels.isotropic(0) + 0.7 * channels.physical()`.
- `decay_channels` restricts generated decays. Keep competing models in
  `interactions` so the inclusive parent width still controls propagation.
- `event_factor` on Weighter multiplies the whole-tree physical density. It
  does not change the proposal or repair missing support.

`PowerLaw` bounds apply both to its proposal density and to its physical flux.
Use separate distributions when the physical and proposal ranges differ.
A zero physical density yields zero weight; missing proposal support and
invalid densities remain errors. Equal physical and proposal factors retain
the native support checks even when they cancel numerically.

`Propagated()`, `Fixed()` and `ExternalBounds()` choose the physical weighting
contract. `PropagatedFromCreation()` additionally includes survival from the
creation point to the injection bounds. The Python Weighter forwards this
choice to its physical process. Native callers must set the physical mode
explicitly; an injection-only flag does not enable physical survival.

Creation survival requires a sampler that retains the source point. Primary
bounded/physical samplers retain an external source, secondary bounded/physical
samplers retain the parent vertex, and FixedTargetPosition defines its own
source. The compiler rejects volume, column-depth, range, decay-range,
point-source and unverified custom position samplers with this mode.
Point-source clipping can move the origin, so use a bounded sampler with an
explicit source for this case.

## Expansion and generation

Declare secondary Vertices in `secondaries`. The parent uses
`expand=(siren.expand.child("N4"),)` to continue a named daughter, or
`siren.expand.depth_below(n)` for a depth limit. Unlisted daughters terminate.
Duplicate secondary types and unreachable or missing declarations fail during
configuration. An optional `continue_if(tree, parent, index)` can veto expansion.
An empty expansion list ends the chain at that vertex. When the primary and
a secondary share a particle type, each retains its own expansion policy.

`generate(N)` returns N successful trees, subject to `max_attempts`. The native
attempt count includes failed attempts and accumulates across calls until
`reset()`. Compute weights after all batches have been generated, using that
combined denominator. The configured event count is independent of retry limits
and remains the default for later calls and loaded injectors.
Use `on_failure="raise"` to expose non-geometric failures immediately and
`on_shortfall="raise"` to require the requested count. `report()` exposes the
attempt and failure counts. Validation probes use a separate random stream.
Probe diagnostics are warnings; an inapplicable template does not certify closure.

Use either Vertex process declarations or legacy process keywords. Mixing them,
unknown Weighter override keys, and pooled processes with different weighting
modes are errors. The first injector defines the shared physical target for
`Weighter(first, second, ...)`; explicit physical overrides replace that target.
Construct a new Injector to change generation distributions after it is built.

## Persistence

Native process modes, mixtures and declarative expansion rules survive saving
and pickling. Injector archive version 2 adds expansion rules; version 1 already
contains the random-engine state. Version 0 archives use the caller's fresh
random engine. SecondaryExpansion version 1 distinguishes primary and secondary
policies; its version-0 archives retain their original policy by particle type.

`Injector.save()` stores generation processes. A Vertex's separate physical
declarations require `Weighter.save()`. Python interaction/distribution
subclasses and `event_factor` cannot be saved as native physics; unsupported
state raises explicitly. A picklable stopping function or `continue_if` predicate
can travel in a Python Injector pickle, separately from the native state.
Native archive callbacks
remain unsupported. Loading restores the random stream rather than restarting
from its seed. `Injector.load(path)` constructs an injector; `injector.load(path)`
replaces an existing instance. The constructor's `max_attempts` setting limits
each `generate()` call and survives Python pickling. Native archives restore
the default retry policy. Weighter pickling resolves its native generation and
physical processes first; incomplete configurations and live injector callbacks
therefore fail explicitly, even if the weighter has not yet been used.

[Compatibility changes](configuration-compatibility.md) records intentional
differences from the legacy interface and the initial U04 draft.
