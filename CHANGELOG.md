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
- MesonProduction three-body rest-frame sampler: the acceptance is weighted by the E_phi kinematic band width, so the accepted (E_nu, E_phi) density is proportional to |M|^2 on the Dalitz plane as FinalStateProbability and total_width assume. A uniform-in-band proposal with unweighted acceptance carries an extra 1/band(E_nu) factor that oversamples both Dalitz edges, biasing MesonThreeBodySIRENDecay and BiasedMesonThreeBodyDecay event samples (E_nu-marginal chi2/dof vs the physical marginal: 130 unweighted, 1.2 weighted; the Recursive2Body FinalStateProbability was verified against the engine Jacobian to 1e-14 and is unchanged).
- Directed 2-body Overlap regime: the sampler now draws from the smaller of the kinematic and bounding cones, which provably contains the cone-intersection lens, so sampled daughter directions cover the whole lens the density reports as 1/omega_eff. The retired tight-cap proposal could enclose only part of the lens (a reproduced geometry left 57% of the lens unreachable), truncating the directed channel's support below its stated density and breaking Sample/Density closure for DetectorDirected2BodyChannel and DetectorDirected3BodyChannel.
- Fixed-mode vertices charge the channel-selection probability (selected rate over total rate) in both the physical and generation final-state factors whenever the vertex's interaction collection offers more than one channel, matching the rate-weighted channel selection the injector always performs. Single-channel collections are untouched (the factor short-circuits to exactly 1.0), and weighting raises WeightCalculationError when competing channels have no computable total rate.
- The PionKaon flux loader converts tabulated counts to nu/m^2/GeV/POT with the BNB factor order (divide by the 50 MeV bin width, multiply by MeV per GeV and by cm^2 per m^2). Dividing by the combined product understated the tabulated flux by exactly 1e14; current in-tree consumers sample only the normalized shape, which is unaffected.
- Vector-portal (Dutta-Kim) coherent upscattering chi N -> chi' N includes the coherent nuclear enhancement Z^2 in the cross section; the factor was computed but never applied, leaving single-nucleon strength despite coherent nuclear kinematics (argon: 324x). Kinetic mixing enters as epsilon^2 exactly once and the Helm form factor is unchanged, so the Q2 spectrum shape and the off-shell narrow-width density (where Z^2 cancels in differential over total) are unaffected. Applies to VectorPortalUpscatteringXS and VectorPortalOffShellXS via their shared kinematics core.

### Fixed

- Zero-argument Injector.ResetInjectedEvents() overload; repairs the Python reset_injected_events() wrapper.
- Injector.load() secondary-process iteration bug (missing .items()).
- Injector archives carry FailedEvents, and process archives carry the vertex weighting mode (each with a class-version bump; older archives load with the previous defaults). PhysicalProcess copy and move preserve the weighting mode.

### Added

- Typed exceptions exported via siren.utilities and registered with RuntimeError base.
- ConvertDensity pybind binding.
- Fixed-seed golden-physics regression harness (tests/python/test_golden_regression.py).
