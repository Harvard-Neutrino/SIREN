# Native generation and weighting errors

See [the native phase-space contract](native_phase_space.md) for complete examples.

## Configuration

`ConfigurationError` or `AddProcessFailure` means a process cannot be used as
configured: for example, no vertex distribution, mismatched primary/secondary
processes, invalid mixture weights, or material scattering at zero momentum.
Correct the configuration before generating more events.

## Measure-compat

`MeasureCompatibilityError` means two densities do not have a supported common
pointwise measure. Declare the actual measures and use a supported proposal;
changing a label alone does not convert a density.

## Weight-calc

`WeightCalculationError` rejects missing proposal support, invalid probabilities,
empty failed events, or weight overflow. Zero physical support is valid and
receives zero weight. The breakdown API records unusable totals as NaN.

## Injection failure

`InjectionFailure` identifies a failed sampling attempt, such as a kinematically
closed final state or no intersection with a required volume. Inspect the
injector's failure ledger and retain the attempt in normalization.
