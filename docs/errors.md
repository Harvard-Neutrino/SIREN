# Native generation and weighting errors

Error messages end with a tag such as `[siren-docs: errors#configuration]`
that names a section below. Examples are in
[Native phase-space proposals](native_phase_space.md).

<a id="configuration"></a>
## Configuration errors

`ConfigurationError` or `AddProcessFailure` means a process cannot be used as
configured: for example, it has no vertex distribution, its primary and
secondary processes do not match, its mixture weights are invalid, or it asks
for a material interaction at zero momentum. Correct the configuration before
generating more events.

<a id="measure-compat"></a>
## Measure compatibility errors

`MeasureCompatibilityError` means two densities have no supported common
measure. Declare the measure each model's density is actually differential in
and use a compatible proposal. Changing the declared label does not convert a
density.

<a id="weight-calc"></a>
## Weight calculation errors

`WeightCalculationError` is raised when no injector could have produced an
event, when a probability is negative or nonfinite, when the event is an empty
failed tree, or when the weight overflows. A zero physical density is valid and
gives weight zero. `EventWeightWithBreakdown` reports an unusable total as NaN.

<a id="injection-failure"></a>
## Injection failures

`InjectionFailure` marks one failed sampling attempt, such as a kinematically
closed final state or a ray that misses a required volume. Inspect the
injector's failure ledger, and keep the attempt in the normalization.
