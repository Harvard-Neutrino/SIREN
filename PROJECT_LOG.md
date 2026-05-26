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

### Open work

- Automatic Jacobian conversion for RestFrameSolidAngle<->LabFrameSolidAngle not yet implemented (requires knowing which daughter the density refers to).
- MesonThreeBodySIRENDecay reports Custom (lab-frame variables) instead of Recursive2Body. Refactoring to use recursive 2-body variables would enable mixing with DetectorDirected3BodyChannel in a multi-channel without relying on Custom fallback.
- Scalar Primakoff model (paper Model ii) not yet implemented.
