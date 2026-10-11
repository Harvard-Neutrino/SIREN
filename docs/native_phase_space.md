# Phase-space proposals

Register a proposal for each interaction signature you want to sample through
`SetPhaseSpace`. The injector draws final states from that proposal, and the
weighter compares its density with the physical density in a common measure.

For a two-body decay, register a mixture on an existing injection process:

```python
from siren import injection as inj

process.SetPhaseSpace(signature, inj.MultiChannelPhaseSpace(
    [inj.PhysicalDecayChannel(model, signature), inj.Isotropic2BodyChannel()],
    [0.3, 0.7],
))
```

The [runnable archive example](../examples/phase_space_archive.py) sets up
a Majorana HNL dipole decay, generates and weights a batch, and reloads the
injector and weighter. Run it with `python examples/phase_space_archive.py`.
It uses the C++ bindings `siren.injection._Injector` and
`siren.injection._Weighter` directly, with a fixed vertex and no detector or
beam exposure.

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
its density at a point is the weighted sum over **all** channels. Nested
mixtures and the adapters for existing `Decay` and `CrossSection` models follow
the same rule.

Supply the channels when constructing a mixture. Omitting weights gives uniform
probabilities. Call `mixture.SetWeights(weights)` to replace and normalize the
weights, or pass an empty list to restore uniform probabilities. Python's
`channels` and `weights` properties return list snapshots. Use
`injection.ValidateChannelDensities(mixture, random, detector, record)` to probe
how peer channels evaluate densities on one another's samples.

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
physical density from a producing injector gives weight zero. Physical factors
are not evaluated for injectors without proposal support; their diagnostic
fields are NaN with a `physical factors not evaluated` flag, while a supported
event can still have a finite total. Channel conventions and cancelled names
remain available for skipped injectors, but per-channel densities are marked
unevaluated. Negative or nonfinite evaluated densities raise
errors, and `EventWeightWithBreakdown` reports them as a NaN total with
per-vertex flags. Matching physical and proposal distributions cancel their
finite positive density factors, but their support and validity are still
checked for every event. Cancellation also requires matching interaction
collections: equal position distributions can have different densities when
their interaction rates differ.

Pure decays do not look up material. For a parent at rest, the decay channel is
chosen from partial-width ratios. Built-in fixed, cone and isotropic direction
distributions assign unit probability to the single momentum of a massive
on-shell parent at rest, where the sampled direction is unobservable. Custom
direction distributions should use the same unit probability mass for valid
stationary records, and zero for invalid stationary records. Material
interactions require a finite, nonzero momentum.

## Archives

Built-in channels, mixtures, process phase-space maps, RNG state, attempt and
failure counters, and weighters can be saved to archives and pickled. If loading
a file fails, the receiving object keeps its previous configuration. Archives
written before RNG state was stored remain readable; loading one uses the
engine passed by the caller. The failure ledger's example messages and the last
partially generated tree are not saved.

RNG continuation requires the same C++ standard library. New version-2 RNG
archives identify their encoding; incompatible or malformed legacy state is
refused explicitly. See [random-engine archives](random_archives.md).

An injector with a stopping-condition callback cannot be saved or pickled,
because the callback cannot be restored; this includes callbacks set through
the `siren.Injector` wrapper. Python-defined interactions and distributions keep
their pickled Python state, and state that cannot be pickled raises an error.
Do not load pickles or SIREN archives from untrusted sources.

## Migrating phase-space callers

Construct `MultiChannelPhaseSpace(channels, weights)` with a nonempty channel
list. Omit weights or pass an empty list for a uniform mixture. Construction
and `SetWeights(weights)` normalize coefficients and check their validity.
Replace empty construction followed by assignment and `Normalize()` with
direct construction. Use C++ `Channels()` and `Weights()` or Python
`channels` and `weights` for inspection. Use `SetWeights()` to update
coefficients and construct a new mixture to replace channels. The third
constructor argument, `allow_incompatible`, and C++ `allow_incompatible_`
field have been removed. Every mixture needs a common topology and measure.

`ValidateChannelsDetailed()` still reports convention diagnostics. The
`ValidateChannelDensities` method is now a free function taking the mixture as
its first argument, in `siren::injection` and `siren.injection`. It samples
each channel and probes the other channels at those points.

In C++, include `SIREN/injection/FinalStateProbability.h` in place of
`SIREN/injection/WeightingUtils.h`. Use `ProcessFinalStateConvention` in
place of `SelectedFinalStateConvention`. Use
`ProcessFinalStateProbability` for the process's weighting factor in the
requested convention, in place of the `SelectedFinalStateProbability`
overloads and the convention-taking overloads of `CrossSectionProbability`
and `CrossSectionProbabilityWithPhaseSpace`. Configure the process with its
weighting mode, interaction collection and any signature-specific proposal.
`ProcessGenerationProbability` supplies the injector's generation factor in
the process's own convention. The collection-level diagnostics
`CrossSectionProbability`, `CrossSectionProbabilityWithPhaseSpace`,
`ChannelSelectionProbability` and `FixedVertexChannelSelectionProbability`
remain available without a convention argument, including in Python.

Mixture archives retain their version-0 layout and exact normalized weights.
The legacy `AllowIncompatible` field is read for compatibility and new
archives write it as false. A legacy mixture with that field set to true
loads if its topology, measure and weights pass validation. A mixture with
fatal convention diagnostics fails to load.

## Tests

Run the focused checks against a build whose Python package and compiled
libraries match the source:

```bash
python -m pytest tests/python/test_native_phase_space.py -q
ctest --test-dir build -R 'PhaseSpace|TwoBodyKinematics|Injector|Weighter' --output-on-failure
```

Momentum values and raw record conversions are documented in [four-momentum and records](four_momentum.md).

For selected decay channels, inclusive parent widths, and opt-in upstream
survival, see [native decay accounting](native_decay_accounting.md).
