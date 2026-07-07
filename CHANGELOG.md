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
- BiasedMesonThreeBodyDecay biases the (nu V1) pair direction in the meson rest frame into the detector cone, keeping the physical Recursive2Body factorization, so FinalStateProbability is the physical density times the constant 4pi/Omega_cone and the importance weight is exactly Omega_cone/4pi. The replaced lab-direction construction forced the neutrino collinear with V1, violating four-momentum conservation by about 0.1 GeV and returning a zero generation density for most of its own boosted samples (biasing every pion-decay event weight in scan_cone_width.py and VectorPortal_SBND_dk2nu.py). Corrected records conserve four-momentum to machine precision and the class passes check_closure alongside the rest of the model family in fixed-cone configurations; the adaptive cone correlates direction with energy, outside the single-coordinate gauge's jurisdiction, and is certified by the per-event identity FinalStateProbability = physical density x 4pi/Omega(event).
- The chi-flux builders (compute_chi_flux, compute_chi_flux_from_dk2nu) spread each V1 -> chi chi' decay over its exact relativistic box spectrum: two-body kinematics fix E_chi_rf = (m_V1^2 + m_chi^2 - m_chi_prime^2)/(2 m_V1), and the isotropic decay makes the chi lab energy uniform on [gamma(E_rf - beta p_rf), gamma(E_rf + beta p_rf)], deposited with exact partial-bin overlap (beta -> 0 collapses to the single bin at gamma E_rf). The tabulated builder converts each per-GeV source density to a node weight over the node's own energy interval before spreading, so the output stays per GeV and its integral matches the input folding; the dk2nu builder deposits per-meson counts unchanged. Collapsing every chi to a single energy with a massless-partner rest-frame convention produced per-node spikes with the wrong center (for m_V1 = 17 MeV and m_chi = m_chi_prime = 8 MeV at E_V = 5 GeV the box is 1.66 to 3.34 GeV against a single point); the per-V1 multiplicity is unchanged.

### Fixed

- Zero-argument Injector.ResetInjectedEvents() overload; repairs the Python reset_injected_events() wrapper.
- Injector.load() secondary-process iteration bug (missing .items()).
- Injector archives carry FailedEvents, and process archives carry the vertex weighting mode (each with a class-version bump; older archives load with the previous defaults). PhysicalProcess copy and move preserve the weighting mode.
- check_closure reports moment_z keyed by the coordinate it measured instead of stamping one z-score under every declared DensityVariable name.

### Added

- Typed exceptions exported via siren.utilities and registered with RuntimeError base.
- ConvertDensity pybind binding.
- Fixed-seed golden-physics regression harness (tests/python/test_golden_regression.py).
- Results.save(pot=...) records protons-on-target as the pot attribute of the saved HDF5 Events group; pot with save_hdf5=False raises ValueError.
