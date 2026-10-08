# Error reference

SIREN error messages that originate from the C++ engine end with a tag like
`[siren-docs: errors#configuration]`. Find the token after `errors#` in the
message, then read the matching section below.

## errors#configuration

**What raises it.** Injector/Weighter process configuration mismatches,
raised as `ConfigurationError`:

- the `Weighter`'s primary physical process does not match the `Injector`'s
  primary process ("does not match the injector primary process");
- a secondary physical process does not match its injector's secondary
  process, or has no matching injector secondary process at all
  ("has no matching injector secondary process");
- the injector's secondary processes and the weighter's secondary processes
  are not in one-to-one correspondence ("no one-to-one mapping between
  injection ... and physical ... secondary processes");
- a `MultiChannelPhaseSpace` has a `weights` list whose size does not match
  its `channels` list, has a weight that is not finite and non-negative, has
  weights that do not sum to 1, or has no channels at all
  ("call Normalize() or use the validating constructor" /
  "MultiChannelPhaseSpace has no channels").

**The fix.** Ensure the physical processes passed to the `Weighter` mirror
the `Injector`'s processes exactly (same primary type, same set of secondary
types, matching collections). For a hand-built channel mixture, call
`MultiChannelPhaseSpace.Normalize()` (or use the validating constructor)
before using it -- do not assign `weights` and `channels` separately without
normalizing. See `docs/quickstart.md`'s closure check and this repository's
config/closure gauges for a way to catch this before a full run.

## errors#measure-compat

**What raises it.** A `MeasureCompatibilityError` when phase-space densities
cannot be compared pointwise in one convention. Common cases are:

- channel topologies disagree ("Topology mismatch: channel 0 ... is
  `<Topology>` but channel `i` ... is `<Topology>`");
- a channel's `Measure()` cannot be converted to the mixture's common measure
  within the shared topology ("Measure incompatibility: channel `i` ... uses
  `<Measure>` which is not convertible to `<Measure>` within `<Topology>`
  topology");
- code asks to integrate an explicit azimuth at one sampled point ("cannot
  marginalize an explicit azimuth pointwise; evaluate the mixture in an
  explicit-azimuth measure");
- `PhysicalProcess::SetPhaseSpace` is given a same-topology, same-family
  mixture measure that the model density cannot reach pointwise ("Phase space
  registered for a signature is not convertible from that interaction
  model's final-state convention"); or
- multiple models for one signature declare conventions with no common
  pointwise density measure.

The scattering measures distinguish an azimuth-integrated marginal from its
explicit-azimuth completion: for example `MandelstamQ2` is `dQ^2`, while
`MandelstamQ2Phi` is `dQ^2 dphi`. Conversion is directional. A marginal that
declares the omitted azimuth uniform can be lifted to its explicit completion
by multiplying the density by `1/(2*pi)`. The reverse operation is an
integral, not a pointwise Jacobian, and therefore throws.

Mixtures elect a common measure that every channel can reach in that
direction. Consequently, one explicit-azimuth channel forces an explicit
common measure and integrated channels lift into it. Other convertible
coordinate changes auto-convert silently.

**The fix.** Give every channel in a mixture the same `Topology()` and a
`Measure()` that can convert into the elected common measure. Check
`PhaseSpaceDensityConvertible(topology, from_measure, to_measure)` in that
direction; do not assume convertibility is symmetric. For a physical adapter
whose model metadata cannot express the density precisely, use the
`(model, signature, convention)` constructor to declare it explicitly. See
`docs/units_conventions.md` for the rest-frame-vs-lab-frame distinction
between `SolidAngleRest` and `SolidAngleLab`, another common source of an
inconvertible pairing.

## errors#weight-calc

**What raises it.** A `WeightCalculationError` when `Weighter::EventWeight`
cannot form a usable weight for an event: the accumulated
`generation_probability` is `<= 0` or non-finite, or the accumulated
`physical_probability` is non-finite ("unusable probabilities for injector
... generation_probability=..., physical_probability=...").

An event with `physical_probability == 0` (out-of-physical-support kinematics)
is NOT an error: it weights to exactly `0.0` by design, since that drives the
weight reciprocal to `+inf` correctly. Only a non-finite or non-positive
generation probability, or a non-finite physical probability, raises.

**The fix.** If you see this for an event that should be in-support, check
the offending injector's distributions and the event's kinematics: a
distribution returning a non-finite or negative density anywhere in its
support will surface here the first time that region is sampled. Run
`siren.check_closure` (`docs/quickstart.md`) on the model responsible for the
vertex to find which distribution's density disagrees with its sampler.
