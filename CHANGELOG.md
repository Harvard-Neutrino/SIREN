# Changelog

## Unreleased

### Breaking (loud-by-design)

- Misconfigured primary/secondary processes now raise AddProcessFailure instead of terminating the interpreter with exit(0).
- Weighter initialization mismatches raise ConfigurationError instead of a debug-only assert or printed message; missing secondary types in weighting raise instead of returning empty/zero results.
- MultiChannelPhaseSpace: channel/weight length mismatches and un-normalized weights raise ConfigurationError at construction or use, instead of silently assigning leftover probability to the last channel; non-convertible measure combinations raise MeasureCompatibilityError instead of a one-shot stderr warning (allow_incompatible=True opts out).
- Zero or non-finite generation density during weighting raises WeightCalculationError instead of returning an infinite weight. Zero physical density still yields weight 0.

### Physics-affecting fixes

- Event weights are normalized by the realized injected-event count instead of the configured target count (falls back to the configured count when weighting precedes generation).
- ConvertDensity now rejects non-Decay2Body SolidAngleLab conversions as unconvertible (MeasureCompatibilityError) instead of silently applying the parent-rest-frame two-body boost Jacobian in the wrong frame; Decay2Body SolidAngleLab conversions are unchanged.

### Fixed

- Zero-argument Injector.ResetInjectedEvents() overload; repairs the Python reset_injected_events() wrapper.
- Injector.load() secondary-process iteration bug (missing .items()).

### Added

- Typed exceptions exported via siren.utilities and registered with RuntimeError base.
- ConvertDensity pybind binding.
- Fixed-seed golden-physics regression harness (tests/python/test_golden_regression.py).
