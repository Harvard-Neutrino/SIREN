# SIREN Project Log

## 2026-05-26: Dutta-Kim chain infrastructure + end-to-end validation

### Completed

- **Fixed CrossSectionDistributionRecord pybind copy policy**: removed `shared_ptr` holder from pybind registration. CSDR is always stack-allocated and passed by reference; the holder caused a copy-policy failure when the pyDecay trampoline crossed the C++/Python boundary.

- **Extracted PhaseSpaceConvention to dataclasses**: moved the enum from `injection/PhaseSpaceChannel.h` to `dataclasses/PhaseSpaceConvention.h` so both `interactions` and `injection` can use it without circular dependencies.

- **Added Convention() to Decay/CrossSection base classes**: virtual method with DensityVariables-based fallback + warning. All Dutta-Kim models (VectorPortal.py, MesonProduction.py) override explicitly. PhysicalDecayChannel/PhysicalCrossSectionChannel now delegate to the model's Convention().

- **Implemented automatic Jacobian conversion in MultiChannelPhaseSpace::Density**: dispatches BjorkenXY<->MandelstamST, Recursive2Body<->Dalitz<->HelicityAngles via PhaseSpaceJacobian.h transforms. Not exercised by Dutta-Kim chain (conventions match at every vertex), but ready for mixed-convention use.

- **Fixed DetectorDirected2BodyChannel Sample/Density consistency**: Sample() was falling back to isotropic when directions were kinematically forbidden, but Density() returned 0 for those events. Fix: resample directions from target geometry until kinematically valid.

- **End-to-end chain test with SBND**: generates 5 events through pion->V1->chi->chi'->V1_signal with detector-directed biasing at each vertex, then verifies finite positive weights. Uses SBND GDML geometry with Ar40 targets.

- **EventWeightBreakdown diagnostic**: `Weighter.breakdown(event)` returns per-vertex VertexWeight with interaction_probability, position_probability, physical_probability, generation_probability. Immediately identifies which vertex causes inf/NaN weights.

- **Exposed injection bounds and weighting utils to Python**: PrimaryInjectionBounds, SecondaryInjectionBounds, CrossSectionProbabilityWithPhaseSpace, ChannelSelectionProbability now accessible from pybind.

### Key design decisions

- Fiducial box for biasing MUST be placed at the detector center (use Placement), not at the coordinate origin.
- Every multi-channel SHOULD include a physical/isotropic fallback channel (~1% weight) for efficiency, though it is not required for correctness after the Sample/Density fix.
- Stopping condition must prevent chi->chi'->chi recursion: only process chi from V1_prod decay, not from chi' decay.
- Models report Convention() explicitly (Custom for off-shell/non-factorizable processes). The base class fallback infers from DensityVariables but warns.

## 2026-05-26: Topology/Measure refactor + test suite

### Completed

- **Split PhaseSpaceConvention into Topology + Measure**: two orthogonal enums replace the single convention enum. Topology (Decay2Body, Decay3Body, DecayNBody, Scatter2to2, Scatter2to3, Unspecified) describes the structural shape. Measure (SolidAngleRest, SolidAngleLab, Recursive2Body, DalitzPair, HelicityAngles, MandelstamQ2, BjorkenXY, Unspecified) describes the density parameterization. Topology mismatches are hard errors; measure mismatches within the same convertibility group are auto-Jacobian converted.

- **Added SolidAngleRest<->MandelstamQ2 Jacobian**: `|dQ2/d(cos_theta)| = 2*p_CM^2`. Enables combining detector-directed and physics-model scattering channels in the same MultiChannelPhaseSpace for Scatter2to2 topology.

- **Comprehensive phase space validation test suite**: 28 C++ tests + 31 Python tests covering Jacobian invertibility, normalization integrals, cross-measure agreement, closure tests, sampling-density consistency, phase-space coverage, boundary conditions, convertibility group tables, and mock physics model integration.

### Key design decisions

- Topology determines compatibility (hard error on mismatch). Measure determines convertibility (auto-Jacobian within same group).
- SolidAngleRest is in a DIFFERENT convertibility group from Recursive2Body/Dalitz within Decay3Body topology. They use fundamentally different phase-space factorizations and cannot be pointwise converted.
- For Scatter2to2, all four measures (SolidAngleRest, SolidAngleLab, MandelstamQ2, BjorkenXY) are in the same group and fully interconvertible.
- Legacy PhaseSpaceConvention is kept with TopologyFromConvention/MeasureFromConvention converters. Convention() on PhaseSpaceChannel has a default implementation derived from Measure().
- All Decay/CrossSection.SampleFinalState implementations MUST produce lab-frame momenta (boost from rest frame to lab using parent 4-momentum). This is the established convention in VectorPortal.py, MesonProduction.py, and the C++ decay models.

### Open work

- RotatedFinalStateChannel (generic biasing by rotating the full final state in the rest frame) is designed but not implemented.
- RestFrameSolidAngle<->LabFrameSolidAngle conversion for Decay2Body is implemented but requires knowing which daughter the density refers to (daughter_index).
- MesonThreeBodySIRENDecay reports Custom (lab-frame variables) instead of Recursive2Body. Refactoring to use recursive 2-body variables would enable mixing with DetectorDirected3BodyChannel.
- Scalar Primakoff model (paper Model ii) not yet implemented.
- Cross-measure conversion in MultiChannelPhaseSpace for Scatter2to2 (SolidAngleRest + MandelstamQ2) is implemented but not yet exercised by any real physics model combination.

## 2026-05-31: Option B (volume-aware directed channels) + directed-fallback optimizer fixes

Branch `interface-redesign`. Builds on the architecture-alignment work, the R5 regime tests (`b38dd2c6`), and the nestable sub-mixture channel (`4618d106`). Full research summary + design: `docs/volume_aware_directed_channel_plan.md`. Commits `a75d1864`, `b1a6b402`, `b35f305f`, `f49e41a3`, `87f515f4` (interleaved with the user's vertex-distribution commits `b60d38a3`, `b6b2ef84`).

### Completed

- **Disjoint-tiling directed channels ("Option B")** (`a75d1864`). Overlapping/nested directed channels are redundant (the mixture already counts overlap once, so it is an efficiency/optimizer-conditioning problem, not bias). Strategy: present the reachable region as a DISJOINT partition of directed sub-channels (a non-degenerate basis). Reuses the existing `BooleanGeometry(UNION/SUBTRACTION)` fed to the EXISTING directed channels via the virtual `Geometry` interface, so C1 (Sample==Density) holds with zero channel changes and 3-body/scatter get it for free; bound `Geometry::GetWorldBoundingBox` to Python for the builder. New `python/directed_tiling.py` (grid cells, subtraction shells, union, angular sectors, `dedup_contained`, `cluster_by_cone_overlap`), all wrapped in `NestedMixtureChannel`. New `DetectorDirectedAngularSectorChannel` (the one new channel; tiles the cone into (u,phi) bins). New `volume` ctor arg on `DetectorDirected2BodyChannel` (skip the MC volume estimate + viability guard for composite tiles). Opt-in `--target-set {all_13,fiducial_only,tiled,union}` on the example.

- **Canonical Kleiss-Pittau default** (`b1a6b402`). `optimize_chain_weights` defaults to `update_rule="alpha_sqrt_W"`. The memoryless `sqrt_W` cannot turn off a channel whose density already covers the support (near g=f every covering channel has W~1, so alpha~sqrt(W) parks it); the multiplicative canonical rule decays it. Toy regression in `test_optimizer_kp.py`.

- **DirectingActive variance attribution** (`b35f305f`). A directed channel returns the same isotropic 1/4pi density in its fallback regimes; `DirectingActive(record)` reports per-event whether it genuinely directs vs sits in that shared fallback. The example's `attribute_directing_variance` (`--attribute`) splits each vertex's variance into W_directing / W_floor / W_phys; `--pion-energy` widens the V1 cone to make directing active at the primary vertex.

- **discount_fallback optimizer option** (`f49e41a3`, default True). Credit a directed channel's KP variance only on directing-active events, so a pure-fallback director gets W->0 and is driven to min_weight regardless of the update rule (consolidating isotropic coverage onto the physical channel). Exact because the fallback density is exactly 1/(4pi) and the regime depends only on the EVENT (parent kinematics), not the sampled daughter direction.

- **Group directed channels** (`87f515f4`). `DirectingActive` is now a virtual on `PhaseSpaceChannel` (default true), overridden by the directed channels and aggregated by `NestedMixtureChannel` ("any member directs"). `group_directed_channels(mc)` (optimize.py) wraps a flat mixture's directed channels into one NestedMixtureChannel, preserving g(x) exactly; the optimizer then has a single "direct vs physical" knob that floors at one min_weight (N-fold smaller residual) and converges faster. Opt-in `--group-directed`.

- **Verification.** C1 closure holds for spatial/angular/union/subtraction tilings (`UnitTest_DirectedAngularSector` + `test_directed_regimes.py`). 2-body ESS (separated regions, full sphere): union +44% / tiled +32% vs a single gap-spanning fiducial, both beat a flat nested set. Monoenergetic chain A/B: all target sets give an identical rate (unbiased); the collimated chain is directing-inert (R2 confirmed) so the ESS win lives in the 2-body tests; `--optimize` drives the inert directed channels to min_weight (physical -> 0.98), eff ~90%. Python suite 497 passed; C++ at baseline.

### Key design decisions

- Disjoint tiling, not union, is the strategy; union/`BooleanGeometry` are only construction tools (dedup nested inputs, disjointify overlaps).
- Reuse `BooleanGeometry` rather than write new union math in the injection layer: the CSG ray-walk IS the interval-merge, so C1 is automatic via virtual dispatch.
- `DirectingActive` default-true on the base means "always a genuine proposal"; only directed channels override it. This lets a generic group/optimizer treat physical/isotropic channels as always-credited without special-casing.
- The three optimizer levers stack: `alpha_sqrt_W` (correct fixed point) + `discount_fallback` (zero the fallback credit) + grouping (one knob, N-fold smaller residual).

### Open work

- `optimize_chain_weights` does NOT recurse into a group's inner weights (only `optimize_multichannel_weights` does, via `_submixture`). Fine for the collimated chain (in-group weights stay uniform isotropic, which is correct there); a non-collimated multi-target chain that must rank in-group targets would need the `nested`/`_submixture` handling mirrored into `optimize_chain_weights`.
- The dk2nu (variable-energy, mixed-regime) chain has not been measured with `attribute_directing_variance` after optimization. That measurement is the trigger for per-event channel masking ("option 3"): masking fallback directors per-event in the mixture selection (C1-clean because the regime is event-determined), a step toward MadNIS-style point-dependent alpha_i(x). Only worth building if the post-optimization W_floor at the primary vertex is a meaningful share of the vertex variance.
- A genuinely non-collimated full-chain ESS A/B (lighter chi / lower boost) to show the separated-target gain end-to-end is not yet done; the gain is currently demonstrated only at the 2-body level.

## 2026-06-02: Channel-weight optimizer moved into C++ + failure-mode policy

Branch `interface-redesign`. The multichannel Kleiss-Pittau weight optimizer was refactored from a Python add-on (`python/optimize.py` reaching into the C++ objects, re-deriving everything) into a first-class C++ capability on `MultiChannelPhaseSpace`, then extended with inter-target tuning and a selectable failure-handling policy. Commits `d2a5bc52`, `601999c5`, `60ccd050`, `49b87d0c`, `17197812`, `baa8366f`. Behavior-preserving by design (the KP estimator is unbiased for any alpha); every phase gated on closure (`E_g[f/g]->1`) + C++/Python equivalence.

### Completed

- **Phase A -- per-channel density breakdown** (`d2a5bc52`). `MultiChannelPhaseSpace::DensityBreakdown` returns each channel's alpha-weighted, common-measure contribution (sum == `Density()`); `Density()`/`DensityBreakdown()` share one `ComputeContributions` loop. `optimize.py` reads the breakdown instead of re-evaluating every channel from Python -- removes N Python->C++ calls per point and fixes a latent measure mismatch (the old code compared an UNCONVERTED `g_i` against the common-measure `g`). C++ + Python sum-invariant tests.

- **Phase B -- KP accumulate/update in C++** (`601999c5`). `Accumulate`/`UpdateWeights`/`ResetAccumulators` move the per-channel variance accumulation + the Kleiss-Pittau update onto the mixture (mutable, non-serialized accumulators). The per-channel statistic stays the BARE form `W_i = mean(w^2 g_i/g)` -- NOT alpha-weighted; alpha-weighting would shift both rule fixed points and break behavior preservation + `test_optimizer_kp`. C++ `UpdateWeights` reproduces python `_kp_update` exactly (kept as the reference oracle); deterministic equivalence test + inner/outer count lockstep test. `optimize_multichannel_weights` collapses to sample -> Accumulate -> UpdateWeights.

- **Phase C -- accumulate during generation** (`60ccd050`). Injector `GetPhaseSpaces`/`AccumulateEventToMixtures`/`AccumulateSelectionToMixtures` route each tree record to its mixture in C++ via the existing signature map (depth-0 -> primary process, else `secondary_process_map[primary_type]`), replacing the Python `_sig_match` + `_collect_phase_spaces` tree-walking. The `Sample`-selected channel index is genuinely unused (the KP statistic sums over all channels), so no provenance threading / serialization change (option ii, not i). Failure penalty carried as a parallel per-channel selection accumulator folded into `UpdateWeights`. `optimize_chain_weights` becomes a thin generate -> weight -> feed -> update loop; dead `_credit_directing`/`_submixture` removed (`_kp_update` kept).

- **Inter-target (nested-group) tuning** (`49b87d0c`). `optimize_chain_weights` gains `recurse_nested` (default True): it now descends into grouped directed channels (`NestedMixtureChannel`) and tunes the inner per-target weights with the total-event-weight statistic, not just the outer "direct vs physical" knob. **Closes the 2026-05-31 open item.** The C++ recursion already existed (used by `optimize_multichannel_weights`); this wires it through the chain path. Failure penalty stays outer-only.

- **Selectable failure_mode** (`17197812`, default "throughput"; example4 `--failure-mode` + labeled converged-weight dump in `baa8366f`). The chain failure handling in `UpdateWeights` is now a policy knob: "throughput" (`W_i *= 1-f_i`, down-weight lossy channels), "ignore" (no adjustment), "coverage" (`W_i /= 1-f_i`, the original). Default changed coverage -> throughput.

### Key findings (failure-mode sweep, DuttaKim SBND, tiled tile-n 5, group-directed)

- Regime set via monoenergetic pion energy (floor frac = inert-fallback share of directed variance): 0.6 GeV non-collimated (floor~0), 1.5 GeV intermediate (~0.55), 3.0 GeV collimated (~0.84).
- **Non-collimated:** the modes diverge maximally (primary directed-group weight 0.07 / 0.36 / 0.89 for throughput / ignore / coverage; ESS 82 / 81 / 63%). The original `coverage` penalty is the WORST -- it up-weights lossy channels, fighting the success-weighted statistic which already discounts failures.
- **Collimated:** all three identical (directed group inert; `discount_fallback` floors it; no `fail_select` to act on).
- **Production dk2nu flux:** all three within noise (ESS ~50%); floor frac ~0 = directing genuinely useful, kept at ~0.47. So the default change is safe on production, AND the post-opt floor frac being ~0 means the "option 3" per-event masking is NOT triggered -- **resolving another 2026-05-31 open item**.
- **Convergence:** ignore and throughput share the same fixed point (lossy channel -> off); throughput reaches it ~2-3 iterations sooner. The "ignore slower with more samples" effect is real and is the Jensen/geometric-mean bias of the multiplicative update under sampling noise (`sqrt(W)` biased down ~ `1/batch`); throughput's deterministic `(1-f)` deflation is batch-immune. At `opt-iterations 10` the chain is NOT fully converged (~20 iters needed) -- production runs likely under-optimize.

### Key design decisions

- Per-channel KP statistic stays BARE (`W_i = mean(w^2 g_i/g)`); the only numeric change vs the old code is `g_i` is now measure-CONVERTED (no-op where measures agree, a latent-bug fix otherwise).
- Phase C routes by the C++ signature map (the index is never used by the estimator); option i (threading the `Sample` index into tree datums) was rejected -- it would touch the serialized `InteractionTreeDatum`.
- `failure_mode` default = "throughput" at every level (C++ `UpdateWeights`, `optimize_chain_weights`, example4); the one coverage-specific test requests that mode explicitly.
- Accumulators are transient runtime state, NOT serialized (`phase_space_map_` is absent from `Process::serialize`).

### Open work / followups

- Production convergence: characterize the `opt-iterations`/`damping` trade-off (the chain is mid-decay at 10 iters); real runs likely need more iterations or higher damping to reach the optimum.
- A clean non-collimated full-chain tiled-vs-single ESS A/B is still only partial (target differentiation shown via `recurse_nested` inner spread ~3000, not a clean ESS A/B).
- Phase D (AdaptiveMapping) and Phase E (narrow the conversion surface) from `architecture_alignment_plan.md` remain; Phase D is the natural next layer (intra-channel shapes self-tune, same "sampler owns its accumulation + refinement" pattern as the now-C++ weight optimizer; the `Mapping1D` Accumulate/Refine hooks are already in place).
- Verification: Python suite 519 passed; C++ at baseline (the documented pre-existing data-file failures + flaky direction histogram, none touched).
