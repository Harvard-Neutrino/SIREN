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

## 2026-06-02: MiniBooNE detector + BNB beamline wired into the SBN loader

Branch `interface-redesign`. MiniBooNE became a first-class SBN-loader detector: `load_detector("SBN", detector="MiniBooNE")` now returns the BNB beamline (`BooNE_50m.gdml`) + Fermilab site geology + a full MiniBooNE detector/vault enclosure, composited in BNB-frame coordinates -- the same machinery ICARUS/SBND use. Closes the "G4BNB stops at the 50 m absorber; nothing models the detector at 541 m" gap. Enclosure dimensions verified against the primary MiniBooNE papers via a deep-research pass.

### Completed

- **Surveyed position pinned.** Tank center at BNB `(0, 1.89614, 541.34) m`, from G4BNB `bsim::Location(0., 189.614, 54134.0, "MiniBooNE")` in `g4bnb .../NuBeamOutput.cc:136` (comment "PRD says 541 meters..Refined by Zarko"): on-axis in x, 1.896 m vertical, 541.34 m downstream. Matches the value already in `sbn_geometry.py` (which had lacked the source comment). (Source nailed down 2026-06-05 -- see that entry: both the offset and the baseline are published in flux PRD 79 (2009) 072002 Sec. II, and 189.614 cm is the survey-precision form of the published 1.9 m. The original parenthetical here -- "decay-pipe axis ~1.9 m above the tank" -- was dropped: it stated a direction the PRD does not give and that disagrees with the stored +y.)

- **SBN loader extended** (`resources/detectors/SBN/SBN-v1/`): `sbn_geometry.py` gains a `MiniBooNE_local` frame + translation edge to BNB (cited to NuBeamOutput.cc:136), and `DETECTORS["MiniBooNE"]` now mirrors the LArSoft detectors (native frame, origin at tank center). `detector.py` gains a `MiniBooNE` spec (`unwrap=False`, locally generated, no remote URL). `sbn_loader.py` gains `ensure_miniboone_gdml()` + `_build_miniboone_gdml()` which emit the enclosure GDML.

- **Standalone model removed + consumers ported.** Deleted top-level `resources/detectors/MiniBooNE/` (the idealized oil-sphere `densities.dat`/`materials.dat`). example2/example4 switch to `load_detector("SBN", detector="MiniBooNE")` and build the 5 m fiducial inline (`Sphere(5.0, 0.0)`). Legacy example2 ported via a new `SIREN_Controller(detector_model=...)` constructor arg (lets the deprecated controller accept a pre-built GDML composite; `GetFiducialVolume` guarded for the no-`densities.dat` case).

- **Detector spec from the primary papers** (deep-research, adversarially verified): tank inner radius 6.096 m (40 ft); opaque optical barrier at 574.6 cm splitting inner-signal / 35 cm veto oil; Marcol 7 oil rho 0.845 g/cm3, n_D 1.4684, CH2 (H:C ~2); 1280 inner + 240 veto 8-inch PMTs (not modelled geometrically); 13.7 m (45 ft) cylindrical vault; >= 3 m dirt overburden; 541 m baseline / 1.9 m offset. Sources: NIM A 599 (2009) 28 = arXiv:0806.4201 (detector); PRD 79 (2009) 072002 = arXiv:0806.1449 (flux).

- **Cylindrical vault enclosure GDML.** Radial onion: inner oil (r<5.746) -> veto oil (5.746-6.096) -> ~1 cm carbon-steel shell (6.096-6.106) -> air vault cavity (to 6.858) -> 0.47 m concrete wall -> dirt. Vertical stack (sphere on supports, ~1.0 m floor clearance): lower slab / vault air / upper slab / oil chimney / electronics room (air + walls + roof slab) / 3 m dirt berm, carved into a dirt world block. Cylinders are GDML `<tube>`s rotated 90 deg about x to stand vertical; placed flat with the detector oil/steel/cap emitted last so they win their overlaps with the vault air.

### Key design decisions

- **Paper > diagram.** Where the papers give a number it wins: vault inner diameter = 13.716 m (45 ft) and overburden = 3 m (>= 3 m), both NIM A 599 Sec. 1.4, replacing larger enclosure-diagram estimates (~15.26 m, ~2.5 m). Un-published dimensions (heights, wall/slab thicknesses, tank clearance, room, cap) are scaled off the enclosure cross-section (anchor: 40 ft sphere = 913 px) and labelled ESTIMATES in-code.
- Steel-shell thickness ~1 cm is back-computed from the ~37 t shell mass over a 12.2 m sphere (no published thickness). Concrete rho 2.30, dirt rho 1.90 assumed.
- Enclosure placed `unwrap=False` (the dirt world block is a real sector; the surrounding SBN till/atmosphere takes over outside it). Local grade lands at BNB y +8.38 m vs the SBN site grade +7.62 m -- the enclosure overrides the site geology locally.
- The detector frame stays axis-aligned with BNB (identity rotation) so example beam directions (`FixedDirection([0,0,1])`) are unchanged; verticality is handled by per-`<tube>` rotation inside the sub-GDML.

### Verification

- `load_detector("SBN", detector="MiniBooNE")` -> `DetectorOrigin (0, 1.89614, 541.34) m`, identity rotation. Density probes confirm every radial layer (oil 0.845 / steel 7.86 / air / concrete 2.30 / dirt 1.90) and the full vertical stack (floor clearance, slabs, oil chimney, room air, roof, 3.000 m berm). SBND/ICARUS unaffected.
- Tests: `test_sbn_geometry.py` (+`test_miniboone_position`), `test_sbn_loader_offline.py`, `test_controller.py` (+`test_prebuilt_detector_model_is_used`) -- 104 pass.

### Open work / followups

- Detector support legs (6 columns) and the side overflow tank omitted (negligible mass). Berm modelled as a dirt block, not the sloped 31.69 deg frustum.
- Grade reconciliation: diagram-derived local grade (+8.38 m BNB) vs SBN site grade (+7.62 m) differ by ~0.8 m; the enclosure overrides locally, the SBN site model is unchanged.
- The full beamline->detector composite is large (~1000+ sectors); not yet profiled for injection throughput.
- Deep-research report (cited evidence per quantity) saved to the workflow output; not committed.

## 2026-06-04: Interactive geometry-viewer controls + SIREN-backed picker

Branch `interface-redesign`, `python/visualization.py` (uncommitted module). Added interactive controls to `siren.visualization.view()` and replaced pyg4ometry's crash-prone right-click picker with a SIREN geometry query.

### Completed

- **Robust picker (the headline fix).** pyg4ometry's `MouseInteractorNamePhysicalVolume.rightButtonPressEvent` segfaults on the ~1000 m site-geology blocks / world: it loops `vtkImplicitPolyDataDistance` over every appended mesh, none matches on large/degenerate meshes, so `di` stays -1 and `appPolyData.GetInput(-1)` reads out of bounds (C-level crash). Neither `buildPipelinesAppend` (segfault) nor `buildPipelinesSeparate` (AttributeError) gives a working picker -- not patchable from SIREN. Replaced the whole interactor with a custom `vtkInteractorStyleTrackballCamera` subclass whose right-click does `vtkCellPicker.Pick -> GetPickPosition` then queries SIREN (`GeoPositionToDetPosition(GeometryPosition(...)) -> GetContainingSector / GetMassDensity`), exactly like `at()`. Result shown in a `vtkTextActor` and printed. `SetInteractorStyle` fully unseats the buggy style, so it can never fire.
- **Keyboard controls** (one `KeyPressEvent` observer on the same style): `h` help overlay, `b` world bbox, `l` material legend, `c` x/y/z section outlines, `g` gas/low-density visibility, `n`/`m` global opacity, `v` restore all, `o` orientation cube, `k` clip widget; shift+right-click hides the clicked material. Built-in `w`/`s`/`r`/`e`/`q` left to VTK (handled via the separate CharEvent, so our KeyPress observer does not clobber them).
- **In-window overlays:** material colour legend (`vtkLegendBoxActor`, densest-first), auto-scaled axis triad (was a fixed 20 mm = invisible at this mm-scale scene), world bbox outline, pick/help/hint text.
- **New `view()` kwargs:** `legend`, `bounding_box`, `picker`. **New `to_web(model, path)`** wrapping pyg4ometry `exportGLTFScene` (.gltf/.glb) / `exportThreeJSScene` (.html); raises a clear "pip install pygltflib [jinja2]" error if deps are missing.
- Refactored shared export+load into `_load_registry` / `_build_viewer`; the SIREN query is isolated in `_world_to_sector` (unit-testable).
- **Surface boundary fix (from live testing), direction-aware.** A cell-picker hit lands exactly *on* a sector boundary, where the point-only `GetContainingSector` (BVH + `IsInside`) is ambiguous -- so a right-click sometimes named the volume *outside* the clicked surface. The clean SIREN-native fix (no geometric nudge): pass the view ray (camera->point) and use the direction-aware overload `GetContainingSector(GetIntersections(p0, dir), p0)`, which `SectorLoop`s the ray's segments and selects the sector being *entered* (`start_point == 0`, `DetectorModel.cxx:1372`) -- the front (clicked) volume. Public pybind overloads used: `GetIntersections(DetPos, DetDir)`, `GetContainingSector(IntersectionList, DetPos)`, `GeoDirectionToDetDirection`. The geometry-frame view ray is rotated into the detector frame first (identity for MiniBooNE, general otherwise). `_world_to_sector` falls back to the point-only query when no direction is given (so `at()` is unchanged).
- **Density must come from the entered sector, not `GetMassDensity(intersections, p0)`.** That overload (`DetectorModel.cxx:845`) tie-breaks the *same* boundary the opposite way from `GetContainingSector` (`<=0 and >=0` matches the segment whose `end_point==0`, i.e. the sector the ray is *exiting* at p0 -- often the outer vacuum/air), so paired with the entered-sector name it gave the wrong density (e.g. veto-oil volume but steel 7.86, or vault-air volume but steel; "often vacuum" when the exit side is world vacuum). Fix: read `sector.density.Evaluate(Vector3D(geo_xyz))` straight off the sector `GetContainingSector` returned (the non-template `GetMassDensity(inter,p0)` is itself just `sector.density->Evaluate(p0)`), guaranteeing volume/density agree. Verified across every MiniBooNE shell: inward/outward rays now give the entered volume AND its material density (oil 0.845 / steel 7.86 / vault air 0.001225 / wall 2.30); interior points still match `at()`.
- **WORLD box was half-size (pre-existing `to_gdml` bug).** The world `<box>` was sized `wh = max(abs(lo),abs(hi)) + margin` per axis -- correct as a *half*-extent, but GDML `<box>` x/y/z are FULL lengths (as the module docstring and `_solid_xml`'s Box branch already use), so the world spanned only +/-wh/2 and reached ~half as far as needed; the outer sectors poked out. Fixed to `wh = 2*(max(abs(lo),abs(hi)) + margin)`. Verified: for MiniBooNE the geometry reaches 900/200/1300 m from the origin and the world half-extent is now 902/202/1302 m (contains all sectors + the 2 m `world_margin`), vs the buggy 451/101/651 m.
- **bbox-fallback solids were half-size too (same `to_gdml` bug, found chasing an ICARUS "nested volumes shifted" report).** First ruled out the GDML *parser*: SIREN's flattened sector placements match pyg4ometry's accumulated global transforms EXACTLY (delta 0) for the standalone ICARUS GDML (incl. 2-level `volFloorComplete`->`volBuildingFloor_*` and `volDetEnclosure`->`volCryostat`), for synthetic rotated-mother/rotated-daughter/grandchild nests, and the default composite preserves every inter-volume distance. So `DetectorModel::BuildVolume`'s composition (`child_pos = parent.LocalToGlobalPosition(child.pos)`, `child_rot = Q_parent * child.rot.conjugated()`) is correct. The shift was in the viewer's GDML *export*: `to_gdml`'s bbox fallback (any non Box/Sphere/Cyl/Cone/Polycone solid -- here `BooleanGeometry`, used by the cavity/exp-hall walls/floor/warm-vessel/shield/therm-ins) wrote `_bbox`'s HALF-extents as GDML `<box>` FULL dims, so those volumes rendered at half size, shrunk toward their (correct) centroid while the neighbouring exact `Box` volumes (enclosure/overburden/cryostat/roof) stayed full size -- reading as "shifted around". Fixed to `2*half`. Verified on composite ICARUS: volDirtUnderground 100x57.6x100 m, volFloorComplete 17.7x3.55x37.5 m, etc., dims and centres now correct.

### Key design decisions

- **Units: VTK/pyg4ometry are in mm, SIREN geometry in m.** The exported GDML carries `lunit="m"`, which the reader scales x1000 when meshing (confirmed: an 800 m box -> +/-400000 mm bounds). So the picker divides `GetPickPosition()` by `_GDML_MM_PER_M = 1000.0` before querying SIREN. There is no rotation between the GDML/VTK frame and SIREN geometry coords (`to_gdml` writes geo positions directly).
- Always install the custom style (even when `picker=False` or a bare GDML path is passed) so pyg4ometry's broken picker is unconditionally replaced; with no model, right-click just reports the 3-D point.
- Per-material actor map keys off `str(VisualisationOptions)` (the exact key `buildPipelinesAppend` uses for its merged body actors); gas cut reuses `_material_vis_options`'s rho<=0.05 threshold.

### Verification

- Headless (no GL): `_world_to_sector` matches `at()` exactly on tank center / origin / a geology block / a far universe point (no crash); full `view()` build sequence (load -> colour map -> addLogicalVolume -> addClipper -> buildPipelinesAppend -> scene bounds -> scaled axes -> clipper widget -> `_install_controls`) runs clean; control state correct (18-entry legend, 29 body / 9 gas actors, our style installed); `coloured=False`, picker-disabled, and `to_web` dep-error paths all behave. Existing viz functions + 546-test pytest collection unaffected.
- **NOT yet run live:** the interactive `v.view()` window needs a display and cannot be exercised headless (segfaults / no GL context) -- run on the desktop: `m = siren.load_detector("SBN", detector="MiniBooNE"); siren.visualization.view(m)`, then right-click an outer block to confirm no segfault.

### Open work / followups

- Live desktop test of the window + right-click picker on the large outer volumes (the specific case that segfaulted before).
- `pip install pygltflib` (+ jinja2) needed before `to_web` works; the glTF path uses random PBR colours (pyg4ometry limitation) -- a coloured export would need the display-bound `exportVtkGLTFScene`.
- pyg4ometry's clipper pipeline prints harmless `vtkContourLoopExtraction: Input contains no points` for non-intersecting meshes (only with `clipper=True`); not ours, not suppressed.

## 2026-06-05: Sourced the MiniBooNE BNB vertical offset (189.614 cm)

Branch `interface-redesign`. Deep-research pass to find a primary reference for the BNB-beam-axis-to-MiniBooNE-tank-center vertical offset hardcoded as `189.614` cm (the `g4bnb` `bsim::Location` tuple carried into the SBN loader). The in-code comment cited only "PRD says 541 m.. refined by Zarko" and gave no source for the vertical value; this entry pins it down and corrects the 2026-06-02 note above.

### Findings

- **The offset IS published.** MiniBooNE flux paper PRD 79 (2009) 072002 [arXiv:0806.1449], Sec. II states verbatim: *"The axis of the beam, defined by the center of the decay pipe, is displaced vertically from the center of the MiniBooNE detector by 1.9 meters"* -- plus the 541 m target-to-detector-center baseline (also in the detector NIM A 599 (2009) 28 = arXiv:0806.4201). So the very "PRD" the code already names contains the offset, not just the baseline. Earlier note wrongly implied PRD gives only 541 m.
- **189.614 cm has no finer published source.** The exact digits are the survey-precision form of the published 1.9 m, carried from the MiniBooNE `BooNEG4Beam` MC into `SBNSoftware/G4BNB` `src/NuBeamOutput.cc:136` (comment "Refined by Zarko, e-mail April 3rd"). `git log` confirms G4BNB's initial commit is by Zarko Pavlovic himself (MiniBooNE flux physicist) -- the value is his refinement, not a separately published measurement.
- **Underlying measurement.** Fermilab Alignment & Metrology survey of the beamline + sphere-fit tank center in the Fermilab site coordinate system; methodology in Oshinowo, FERMILAB-CONF-02-425 (IWAA 2002) and -03-500 (horn/target, 2003). The numeric coordinates live in the internal 8 GeV Beam TDR (the PRD's ref [6]) and Detector TDR -- not public. Coordinate convention (origin at tank center, z = beam direction, y = vertical) is pinned in the detector NIM paper.
- **Uniquely MiniBooNE.** All other BNB detectors are at y=0 in `NuBeamOutput.cc` (MicroBooNE/SciBooNE/T600 on-axis; SBND offset only in x). The 1.9 m vertical is MiniBooNE's specific surveyed offset, not a global beam tilt.

### Changes

- `sbn_geometry.py` (MiniBooNE_local -> BNB transform): rewrote the comment to quote the PRD sentence, explain 189.614/54134 cm as the survey-precision form, attribute to Pavlovic + the Fermilab survey, and add a sign caveat. The transform provenance string now leads with the PRD citation.
- `sbn_loader.py`: corrected "refined ... by the G4BNB survey" (G4BNB is a simulation, not a survey) -> "from the MiniBooNE/G4BNB beam MC"; noted both numbers are published in the PRD.
- Corrected the 2026-06-02 bullet (above), which asserted "decay-pipe axis ~1.9 m above the tank" -- a direction the PRD does not state and that disagrees with the stored +y.

### Open work / followups

- **Sign unresolved.** Stored as +y (BNB +y = up per the BNB frame), which places the tank center 1.9 m *above* the beam axis; the PRD gives only the magnitude. Confirm the up/down sense against the dk2nu/`bsim` +y convention if it matters physically. The coordinate is a faithful copy of the authoritative G4BNB value, so the placement matches G4BNB regardless of how the prose reads.

## 2026-06-05: Sourced the SBND BNB horizontal offset (73.78 cm)

Branch `interface-redesign`. Companion to the MiniBooNE pass above: deep-research pass to pin a well-documented provenance for the SBND position hardcoded as `[0.7378, 0.0, 110.0]` m in the `SBND_LArSoft -> BNB` transform. The in-code comment said only "G4BNB SBND Location (73.78, 0, 11000) cm" with no source for the off-axis x = 73.78 cm. Conclusion: the value is correct and current, and traces to an internal design document with a clear correction history.

### Findings

- **Primary source confirmed.** `SBNSoftware/G4BNB` `src/NuBeamOutput.cc:145`: `bsim::Location(73.78, 0., 11000., "SBND")` -- detector origin in the BNB frame: x = +0.7378 m (off-axis), y = 0 (on-axis vertically, unlike MiniBooNE), z = 110 m downstream. The stored tuple is a faithful copy. Same value in the G4BNB `README.md` detector table (`sbnd | 73.78 | 0 | 11000`).
- **Authoritative provenance = SBN-DocDB 20891.** The 73.78 cm is the as-designed horizontal offset. SBND Release Notes 09.21.00: *"previous flux files (config F) have the beam spot at a 45.7 cm offset ... According to the current design this offset should instead be -73.78 cm, as mentioned in #79 and described in docdb 20891. The new files (called config H) were produced by Zarko and have the corrected -73.78 cm offset."* The general flux-window methodology is in DocDB 5672; the offset itself is DocDB 20891 (both restricted, `sbn-docdb.fnal.gov`).
- **Correction history (why provenance matters here).** The offset was wrong twice before: a uBooNE-frame ~1.3 m, then 45.7 cm in flux configs F/G, corrected to -73.78 cm in config H (Z. Pavlovic, sbndcode PR #95, ~Apr 2021). SBND flux-files wiki: *"the centre of the window has been changed from (X, Y) = (45.7, 0) cm to (X, Y) = (-73.78, 0) cm, as in the current design."*
- **G4BNB's own default lagged.** G4BNB kept SBND on-axis (0, 0, 11000) in its dk2nu metadata until commit `bdd5f29e` (2025-03-06, J. Paton): "Changed the default position of sbnd ... to be correctly off axis." The off-axis value lived in the flux-production / sbncode layer for ~4 years before being baked back into G4BNB.
- **Independent public confirmations (all agree on 73.78 cm / 110 m).** sbncode `FluxReader/job/run_fluxreader_sbnd.fcl` (`XShift: -73.78 # cm, detector shift along X w.r.t. beamline`); sbncode `EventGenerator/MeVPrtl/config/bnb_kaon_sbnd.fcl` (`[..., 73.78, 0, 11000] # detector in beam frame`, commented `[-73.78, 0, -11000] # beam in detector frame`); G4BNB README; SBND flux-files wiki + release notes.
- **Sign reconciled (no bug).** SIREN stores +0.7378 = SBND origin in the BNB frame (detector-in-beam-frame), matching G4BNB `bsim::Location` and the bnb_kaon_sbnd.fcl "detector in beam frame" tuple. The wiki/FluxReader -73.78 cm is the SAME offset as beam-center-in-detector-frame -- the inverse, fully consistent.

### Changes

- `sbn_geometry.py` (SBND_LArSoft -> BNB transform): replaced the one-line comment with a sourced block citing NuBeamOutput.cc:145, DocDB 20891, the 45.7 -> -73.78 config-H correction (Pavlovic, PR #95), the sign reconciliation, and two caveats. The provenance string now leads with the file:line and DocDB. Added a parallel SBND line to the module-docstring References (DocDB 20891) beside the existing ICARUS DocDB 22998 line.

### Open work / followups

- **Design value, not a published survey.** Unlike MiniBooNE's 1.9 m (published in PRD 79 (2009) 072002), 73.78 cm has no peer-reviewed citation; its root is the internal SBN-DocDB 20891. The as-built survey behind the design is not public.
- **Offset is to the SBND coordinate origin** (LArSoft world origin = cathode plane), not the active-LAr center -- the latter is carried separately in `DETECTORS["SBND"].center_native = [0, -0.59, 2.92]`. Off-axis-ness of the active volume combines both.
- **Cardinal-direction label unverified.** SIREN annotates BNB +x as "horizontal west"; the numeric tuple is faithful to G4BNB by construction, but whether "west" is the correct physical sense of G4BNB +x was not independently pinned (same class of caveat as MiniBooNE's +y up/down).
