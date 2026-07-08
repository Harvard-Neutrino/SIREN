# Changelog

## Unreleased

### Breaking (loud-by-design)

- Model authors must implement `sample(record, random)`; declaring a measure no longer silently selects isotropic sampling. Isotropic two-body decays can explicitly call `sample_isotropic`.
- Misconfigured primary/secondary processes now raise AddProcessFailure instead of terminating the interpreter with exit(0).
- Weighter initialization mismatches raise ConfigurationError instead of a debug-only assert or printed message; missing secondary types in weighting raise instead of returning empty/zero results.
- MultiChannelPhaseSpace: channel/weight length mismatches and un-normalized weights raise ConfigurationError at construction or use, instead of silently assigning leftover probability to the last channel; non-convertible measure combinations raise MeasureCompatibilityError instead of a one-shot stderr warning (allow_incompatible=True opts out).
- Nonpositive/nonfinite generation probabilities and negative/nonfinite physical probabilities raise WeightCalculationError in event and per-process weighting. Vertex factors are checked before multiplication; inverse-weight and final-weight overflow also raise. Zero physical density still yields weight 0.

### Physics-affecting fixes

- Event weights are normalized by the realized injected-event count instead of the configured target count (falls back to the configured count when weighting precedes generation).
- ConvertDensity now rejects non-Decay2Body SolidAngleLab conversions as unconvertible (MeasureCompatibilityError) instead of silently applying the parent-rest-frame two-body boost Jacobian in the wrong frame; Decay2Body SolidAngleLab conversions are unchanged.

### Fixed

- Closure refuses certification when densities depend on sampler-written parameters unavailable to the reference, sizes missing reference storage, uses independent sample variances and a joint covariance test in angular comparisons, and omits unmeasured worst-bin diagnostics.
- DarkNews native and legacy defaults compare by identity, matching authoring models. Sampler audits accept bound instance methods and lambdas while rejecting the default; isotropic mass errors name the model and `SecondaryMasses` hook.
- Authoring models use identity equality by default; the override audit checks `equal`, and DarkNews trampolines forward explicit Python equality overrides.
- Closure uses an independent joint angular reference, reports incomplete checks explicitly, rejects invalid densities, and rechecks mutable models. Coordinate diagnostics use the declared lab, parent-rest, or collision-CM frame.
- Serialization guards recognize Python subclasses even when they define no methods.
- Python weighting rejects invalid final weights, including negative values that underflow during float conversion, subclass results, and custom generation-batch results. Valid zero weights remain usable in diagnostics.
- Weighter save guards inspect fully initialized injectors, including compiled expansion callbacks and Python sampling models.
- EventWeightWithBreakdown follows the scalar weight guards, reporting invalid vertex probabilities and arithmetic overflow with flags and a NaN total. Valid zero physical support still gives a zero total.
- Weighter archives carry a magic+version header tied to the class version, load into a temporary so a failed parse cannot half-mutate the live weighter, name the file in load errors, and still read headerless version-0 archives.
- Injector archives load into a temporary and move-assign, so a failed parse cannot half-mutate the live injector, and the load errors name the file and which parse (headered or headerless) failed, matching the weighter. The archive header version is now tied to the class version on the save side of both the injector and the weighter, rather than a hand-tracked constant. Enum fields validate their range on load: VertexWeightingMode's bound source and the HNL decay/dipole channel enums (HNLDecay and HNLDipoleDecay ChiralNature, HNLDipoleFromTable HelicityChannel) throw a named runtime error on an out-of-range value instead of silently accepting a corrupt archive. DarkNewsDecay's load_and_construct is now static, the form cereal requires, so a concrete subclass would deserialize through it. The trampoline cereal load unpickles the Python state bytes once instead of twice, so a model's __setstate__ side effects no longer run twice on reload.

### Added

- Weighter accepts an optional event_factor(tree) for whole-event physical corrections, with scalar/batch/explain agreement, in-memory copying, and explicit serialization guards.
- Typed exceptions exported via siren.utilities and registered with RuntimeError base.
- ConvertDensity pybind binding.
- Fixed-seed golden-physics regression harness (tests/python/test_golden_regression.py).
