# Changelog

## Unreleased

### Breaking (loud-by-design)

- Misconfigured primary/secondary processes now raise AddProcessFailure instead of terminating the interpreter with exit(0).
- Weighter initialization mismatches raise ConfigurationError instead of a debug-only assert or printed message; missing secondary types in weighting raise instead of returning empty/zero results.
- MultiChannelPhaseSpace: channel/weight length mismatches and un-normalized weights raise ConfigurationError at construction or use, instead of silently assigning leftover probability to the last channel; non-convertible measure combinations raise MeasureCompatibilityError instead of a one-shot stderr warning (allow_incompatible=True opts out).
- Zero or non-finite generation density during weighting raises WeightCalculationError instead of returning an infinite weight. Zero physical density still yields weight 0.

### Physics-affecting fixes

- Bounded vertex distributions (PrimaryBoundedVertexDistribution, SecondaryBoundedVertexDistribution) fail the attempt (InjectionFailure, NoPathThroughVolume) when the sampled ray does not cross the configured fiducial volume, instead of silently placing the vertex along the full unrestricted path through the detector model. The silent fallback generated out-of-volume (e.g. dirt) vertices while the density and bounds reported the in-volume normalization branch; both now agree that a ray missing the volume is outside the distribution's support.
- Event weights are normalized by the number of generation ATTEMPTS instead of the number of successfully injected events (falls back to the configured target count when weighting precedes generation). A failed attempt -- a sampled ray missing a fiducial volume, a kinematically forbidden draw -- is a legitimate zero-weight sample from the generation density, and normalizing by successes inflated every surviving weight by 1/efficiency (a factor ~74 for an isotropic pi+ -> mu+ nu_mu decay chain constrained to the SBND active volume). Failure-free chains are unchanged (attempts equals successes), and the pre-redesign convention (normalize by the requested attempt budget) is restored in spirit. The golden regression archive was re-blessed for this fix: the pinned chain has 598 attempts for its 500 accepted events, so its weights move by exactly 500/598 with per-event shape untouched.
- ConvertDensity now rejects non-Decay2Body SolidAngleLab conversions as unconvertible (MeasureCompatibilityError) instead of silently applying the parent-rest-frame two-body boost Jacobian in the wrong frame; Decay2Body SolidAngleLab conversions are unchanged.
- Distribution cancellation between an injection process and its physical process matches by value (WeightableDistribution::operator==) rather than shared-pointer identity, as the in-tree documentation everywhere already claimed. Weights change only for distributions that are value-equal but distinct instances, which previously evaluated on both sides of the ratio: that ratio was already 1 for them, except for the 0/0 edge pairs that previously threw WeightCalculationError and now cancel cleanly. Distributions that are literally shared (the common case) compare equal exactly as before, so their weights do not move. The documented Weighter.load()-with-live-injectors path, which binds deserialized physical distributions against live injection distributions, now cancels them as intended (deserialized and live instances are never pointer-equal, so identity matching could never fire).

### Fixed

- Weighter archives carry a magic+version header tied to the class version, load into a temporary so a failed parse cannot half-mutate the live weighter, name the file in load errors, and still read headerless version-0 archives.
- Injector archives load into a temporary and move-assign, so a failed parse cannot half-mutate the live injector, and the load errors name the file and which parse (headered or headerless) failed, matching the weighter. The archive header version is now tied to the class version on the save side of both the injector and the weighter, rather than a hand-tracked constant. Enum fields validate their range on load: VertexWeightingMode's bound source and the HNL decay/dipole channel enums (HNLDecay and HNLDipoleDecay ChiralNature, HNLDipoleFromTable HelicityChannel) throw a named runtime error on an out-of-range value instead of silently accepting a corrupt archive. DarkNewsDecay's load_and_construct is now static, the form cereal requires, so a concrete subclass would deserialize through it. The trampoline cereal load unpickles the Python state bytes once instead of twice, so a model's __setstate__ side effects no longer run twice on reload.

### Added

- Typed exceptions exported via siren.utilities and registered with RuntimeError base.
- ConvertDensity pybind binding.
- Fixed-seed golden-physics regression harness (tests/python/test_golden_regression.py).
