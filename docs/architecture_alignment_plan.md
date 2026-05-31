# Implementation Plan: Architecture Alignment (MadGraph / Achilles)

Status: living implementation plan. Companion to
`phase_space_biasing_design.md`. This document turns the architectural
evaluation of SIREN against MadGraph5/MadEvent and Achilles into a
concrete, phased, testable work plan.

Provenance: the evaluation found that SIREN's core importance-sampling
machinery (the multichannel mixture `g = sum_i alpha_i g_i`, the weight
`w = f/g`, the Kleiss-Pittau alpha optimizer, and closure as the operational
analog of unweighting efficiency) is the same well-validated machinery both
peers use. Where SIREN diverges, almost every divergence traces to one root
cause: SIREN's importance target is detector geometry, which neither peer
has. The detector-directed subsystem is a genuine contribution and is out of
scope for "alignment" -- we keep and invest in it. This plan addresses the
divergences that are NOT forced by the detector-geometry goal, i.e. the
places where the two mature generators show a safer or more capable design.

Seven work items, from the evaluation's prioritized list:

| # | Item | Priority | Effort | Risk | Depends on |
|---|---|---|---|---|---|
| 1 | C1 by construction (shared Mapping object) | Highest | L | Med | Phase B |
| 2 | G6: ConvertDensity throws on unsupported pairs | High | S | Low-Med | -- |
| 3 | Adaptive intra-channel sampling (VEGAS-style) | Med-High | L | Med | Phase B, 5 |
| 4 | R2: collapse redundant geometric channels | High | M | Low-Med | -- |
| 5 | Verify Kleiss-Pittau alpha-update rule | Med | S | Low | -- |
| 6 | Narrow the measure-conversion surface | Med | M | Low | 2, 4 |
| 7 | Keep / promote per-vertex weight breakdown | Low-Med | S | Low | -- |

Effort: S = hours-to-2 days, M = 2-5 days, L = 1-2 weeks.

---

## 0. Guiding constraints (read first)

These shape every item below.

**C0.1 -- Density must be evaluable at arbitrary points.** `PhaseSpaceChannel`
(`projects/injection/public/SIREN/injection/PhaseSpaceChannel.h:65-73`)
requires `Density(detector, record)` to return `g_i(x)` for ANY valid point,
not just points this channel sampled. This is non-negotiable: the mixture
density `g = sum_i alpha_i g_i(x)` evaluates every channel at the single
sampled point (`PhaseSpaceChannel.cxx:340-360`). Therefore "C1 by
construction" (#1) CANNOT mean merging Sample and Density into one call.
It means: both must route through ONE shared parameterization object whose
forward map (draw) and density (evaluate-at-x) cannot independently drift.
This is exactly Achilles' `Mapper` (GeneratePoint / GenerateWeight as the
forward and inverse of one declared bijection) and MadEvent's
single-routine-per-channel generation.

**C0.2 -- Unbiasedness is preserved regardless of tuning.** Kleiss-Pittau
guarantees the integral estimate is unbiased for ANY valid alpha and ANY
proposal density that covers the support. Items #3, #4, #5 change variance
and throughput only; they cannot change the answer. This makes them safe to
A/B test by comparing total weight per POT (must agree within statistics)
plus `Var(log W)` and effective sample size (the thing we want to improve).

**C0.3 -- Closure is the regression oracle.** Every change is gated on the
per-vertex closure tests `E_g[f/g] -> 1` (`tests/python/test_dutta_kim_chain.py`)
and the normalization/closure suite (`test_phase_space_validation.py`)
staying green. Closure is SIREN's operational equivalent of unweighting
efficiency; it is what caught the original C1 violation.

**C0.4 -- Build/install/test discipline (CLAUDE.md).** After any C++ change:
`source env.sh`; `cmake --build .` and `cmake --install .` from `build/`
(never skip install, never test the build-dir copy); `ctest
--output-on-failure`; then `python -m pytest tests/python/ -q` from the
source root. ASCII only in source. Do not reference Claude in commits.

**C0.5 -- Pre-existing failures are not regressions.** `UnitTest_DISFromSpline`,
`UnitTest_DipoleFromTable`, `UnitTest_HNLDISFromSpline` fail from missing
data files on this machine; `Cone.SampleDistributionPhi` is a flaky
fixed-seed histogram bound. Baseline before starting: 39/42 C++ pass, Python
suite green.

---

## Phasing and sequencing

```
Phase A  (independent quick wins, de-risk later work)
   #2  G6 throw            #5  KP verify
   #4  R2 collapse         #7  breakdown test
        |
Phase B  (foundation -- prerequisite for #1 and #3)
   Mapping1D abstract base + Python bindings
        |
        +-----------------+
        |                 |
Phase C (#1)         Phase D (#3)
   C1 via shared       AdaptiveMapping + optimizer-loop refinement
   Mapping in models
        |
Phase E  (#6)  narrow measure-conversion surface (needs A's audit + B)
```

Rationale for ordering:
- Phase A items are mutually independent and independent of the Mapping
  refactor. Do them first: they are cheap, and #2 + #4 SHRINK the
  conversion surface and channel counts, making Phases C-E cheaper and
  safer to test.
- Phase B is the single shared prerequisite for the two large items.
- #1 (Phase C) and #3 (Phase D) both consume the Mapping1D base; they can
  proceed in parallel once B lands, but #3 benefits from #5 (the optimizer
  loop is where adaptive refinement hooks in).
- #6 (Phase E) is final cleanup: it depends on #2's audit (which conversions
  actually fire) and #4's collapse (which removes most cross-measure mixes).

---

## Item 2 -- G6: ConvertDensity throws on unsupported conversions

Phase A. The fastest correctness win; also an audit that informs #6.

**Why.** `ConvertDensity` (`PhaseSpaceChannel.cxx:59-306`) silently returns
its input unchanged on the final fall-through (`:305`) and at ~10 internal
precondition guards (e.g. `:77` missing masses, `:90` `:99` degenerate
momenta, `:140` `:149` missing `bjorken_y`). A caller cannot distinguish
"no conversion needed" from "conversion not implemented" from "inputs
missing." Since this density forms the denominator of `w = f/g`, a silent
identity is the same latent-corruption class as a C1 violation. MadGraph and
Achilles cannot have this bug: they never convert between named measures
(common base parameterization).

**Current state.** Three semantically distinct `return density;` cases are
conflated:
1. Legitimate identity: `from == to` (`:66`) and the trivial
   Recursive2Body<->HelicityAngles same-measure case (`:269`). MUST keep.
2. Precondition failure: handled conversion, but the record lacks required
   inputs (the internal guards). Today silent; should be loud.
3. Unsupported pair: no branch matches, falls to `:305`. Today silent;
   should be loud.

**Target.** Distinguish the three. (1) returns unchanged. (2) and (3) throw a
`std::runtime_error` with a descriptive message naming `from`, `to`,
topology, and (for case 2) the missing field.

**Steps.**
1. Add a unit test `UnitTest_PhaseSpaceConvertDensity` (or extend
   `PhaseSpaceChannels_TEST.cxx`) that enumerates every (topology, from, to)
   triple actually exercised by the two chains and the Python suite, and
   asserts each returns finite positive density. This is the audit; capture
   its output as the "supported set."
2. Refactor `ConvertDensity`: tag legitimate identities explicitly; replace
   the precondition guards and the final fall-through with throws. Keep a
   single early `if (from == to) return density;` at the top.
3. Add tests asserting an UNSUPPORTED pair (e.g. a Scatter2to3 SolidAngleLab
   target, currently unhandled) throws, and that a handled pair with a
   deliberately blanked input field throws with the field name.
4. Build/install; run ctest; run pytest. The closure suite must stay green
   (it only exercises the supported set; if a throw fires there, we found a
   real latent gap -- fix the conversion or the channel measure, do not
   re-silence).

**Risk / rollback.** Low-Med. The only way this breaks a green path is if a
chain silently relied on the identity fallback to paper over a genuine
mismatch -- which is precisely the bug we want surfaced. Rollback is a
one-line revert to the non-throwing fall-through, but prefer to fix forward.

**Done when.** ConvertDensity throws on (2)/(3); the supported-set audit test
documents exactly which conversions the chains use; full suite green.

---

## Item 4 -- R2: collapse redundant geometric channels

Phase A. High pass-rate-per-cost win; independent of all other items.

**Why.** The optimizer drives all 13 geometric targets to near-uniform
because it cannot tell them apart -- the textbook Kleiss-Pittau degeneracy
signature (when channels are near-identical, all `W_i` are equal and any
correct KP rule yields uniform alpha). MadGraph and Achilles never see this
because each channel is genuinely distinct by construction (one Feynman
diagram each). R2 realigns SIREN's channel basis with the one assumption KP
needs to be useful.

**Current state.** `resources/examples/example4/DuttaKim_SBND_full_chain.py`
defines `build_geometric_targets(detector_model, fiducial)` (~line 170-216)
returning 13 geometries (1 fiducial + 4 spheres + 8 cylinder segments).
`geo_list = list(targets.values())` (~:348) is looped (`for target in
geo_list`) at the three primary/off-shell channel builders (~:355, :410,
:478), producing 1 physical + 13 directed = 14 channels per biased vertex.
The on-shell intermediate vertices (upscatter, chi' decay) were already
collapsed to physical-only in `fa32b8d0`. Not yet done for the primary
(pion) vertex or the off-shell scatter.

**Target.** `physical + 1 directed(fiducial)` at the primary vertex and the
off-shell scatter; single physical/isotropic channels where directing is
inert. Roughly 7-14x fewer density evaluations per event.

**Steps.**
1. Add a `target_set` parameter (or module-level constant) selecting
   `{all_13, fiducial_only}` so the change is one localized switch and the
   old behavior remains reproducible for the A/B.
2. Point the three channel builders at `[targets["fiducial"]]` when
   collapsed.
3. Re-measure: run `Testing/variance_probe.py` per-vertex `Var(log W)`
   attribution and the chain's total weight per POT, effective sampling
   fraction, and density-eval count, for `all_13` vs `fiducial_only`, both
   chains, AFTER re-running `optimize_chain_weights`.
4. Acceptance gate (C0.2): total weight per POT agrees within statistics;
   `Var(log W)` and effective sampling unchanged-or-better; density-eval
   count drops as expected. Record the numbers in
   `onshell_chain_efficiency.md` / `geometric_sampling_plan.md`.

**Risk / rollback.** Low-Med. Collapsing shifts tuned weights and end-to-end
weight values, so closure-test tolerances must still hold (they should --
closure is unbiased). Rollback = flip `target_set` back to `all_13`.

**Done when.** Both chains run with `fiducial_only`; A/B shows negligible
accuracy change with the expected throughput gain; numbers recorded.

---

## Item 5 -- Verify the Kleiss-Pittau alpha-update rule

Phase A. A precise, checkable variance-efficiency item (NOT correctness).

**Why.** The canonical KP update (Kleiss-Pittau 1994, Eq. 7) and Achilles
(`Multichannel.cc`, `channel_weights[i] * pow(train[i], beta)`) are
MULTIPLICATIVE: `alpha_i_new ~ alpha_i * sqrt(W_i)`, with
`W_i = <(g_i/g) w^2>`. Its fixed point is `W_i = const` -- the proven
variance minimum (every channel contributes equal variance). SIREN computes
the correct bare `W_i` (`optimize.py:318`, `var_contrib[i] += w2 * gi / g`,
no alpha factor) but then sets `alpha_i_new ~ sqrt(W_i)` (`:376`, and the
single-mixture variant at `:100`), DROPPING the leading `alpha_i`. SIREN's
fixed point is therefore `alpha_i ~ sqrt(W_i)`, i.e. `W_i ~ alpha_i^2` --
NOT the variance minimum; high-alpha channels are left carrying more
variance. The feedback sign is correct (high-variance channels still get
up-weighted) and the integral is unbiased either way, so this is a
"leaves variance on the table" issue, not a bug. The `damping` blend does
not change the fixed point (only the path).

**Current state.** `python/optimize.py`: `optimize_multichannel_weights`
(:20-120, update at :100) and `optimize_chain_weights` (:123+, update at
:376, with the per-channel failure-penalty at :329-374). Both use
`new_alpha = sqrt(var_contrib)`.

**Target.** Selectable update mode; default to the canonical multiplicative
rule if the A/B confirms it helps.

**Steps.**
1. Add a `update_rule` parameter: `"sqrt_W"` (current) vs
   `"alpha_sqrt_W"` (canonical, `new_alpha[i] = old[i] *
   sqrt(var_contrib[i] / n_matched)`). Preserve the failure-penalty term in
   both.
2. Add a toy unit test (`tests/python/test_optimizer_kp.py`): a two-channel
   mixture with analytically known `W_i`. Assert the canonical rule
   converges to `W_i = const` and the current rule converges to
   `alpha ~ sqrt(W_i)`; assert both leave the integral unbiased.
3. A/B on both chains: compare final `Var(log W)`, effective sampling, and
   convergence speed (iterations to a stable D = max|W_i - W_j|). Watch the
   `min_weight` floor interaction (the canonical rule can shrink a channel
   faster).
4. Adopt the better default; keep both modes.

**Risk / rollback.** Low. Unbiased either way; switching modes is a parameter
flip. Only risk is re-tuning the example's optimization hyperparameters.

**Done when.** Toy test documents the fixed-point difference; A/B numbers
recorded; default chosen with evidence.

---

## Item 7 -- Keep / promote the per-vertex weight breakdown

Phase A. Low effort; protects an asset both peers lack.

**Why.** `log W = sum_v log(phys_v / gen_v)` with `Weighter.breakdown()`
gives per-vertex variance attribution -- strictly better diagnostics than
MadGraph's or Achilles' global unweighting efficiency / D-convergence.
Neither peer attributes variance to a vertex of a multi-stage chain. This is
how the upscatter was found to be 88% of the variance. Do not regress it;
make it first-class.

**Steps.**
1. Add a regression test asserting `sum_v breakdown_v == log(total weight)`
   (per-vertex factors compose to the total) on a small chain.
2. Promote `Testing/variance_probe.py` to a formal test (this is gap P5):
   a fast, low-statistics per-vertex `Var(log W)` attribution that asserts
   the dominant-vertex ranking is stable, so future refactors that shift
   variance are caught.

**Risk / rollback.** Low. Pure test addition.

**Done when.** Breakdown-sum invariant tested; P5 probe promoted.

---

## Phase B -- Mapping1D abstract base + Python bindings (foundation)

Prerequisite for #1 and #3. No behavior change by itself.

**Why.** `InvariantMassMapping.h` already implements the C0.1 pattern: each
mapping (`BreitWignerMapping` :25, `PowerLawMapping` :71, `TabulatedMapping`
:119, `PropagatorMapping` :227, `UniformMapping` :264) is a struct bundling
`Forward(r)` (draw), `Inverse(x)`, and `Density(x)` (evaluate) on ONE
object. The C++ directed scattering channel already uses this correctly
(`DetectorDirectedScatteringChannel.cxx`: Forward at :457/:460, Density at
:523/:526 -- same mapping type both sides). But: (a) there is NO common base
class (duck-typed structs), and (b) NONE are exposed to Python. So the
Python physics models cannot route through a shared mapping and instead
hand-code sampling and density separately -- the C1 drift surface.

**Target.** A polymorphic `Mapping1D` base the models can hold and the
adaptive mapping can extend, bound to Python.

**Steps.**
1. Define `class Mapping1D` in `InvariantMassMapping.h` with pure virtuals
   `double Forward(double r) const`, `double Inverse(double x) const`,
   `double Density(double x) const`, plus virtual no-op hooks
   `void Accumulate(double x, double weight)` and `void Refine()` (default
   empty; only AdaptiveMapping overrides them -- this keeps #3 a true
   drop-in without a second interface).
2. Make the five existing structs derive from `Mapping1D` and mark methods
   `override`. They are header-only and stateless-per-call, so the change is
   mechanical; the inline instantiations in the directed scattering channel
   are unaffected (they call concrete types).
3. Add pybind bindings in `projects/injection/private/pybindings/` exposing
   `Mapping1D` (base) and each concrete mapping with its current constructor
   (`TabulatedMapping(s_nodes, cdf_nodes, s_min, s_max)`,
   `PropagatorMapping(m2, x_min, x_max)`,
   `PowerLawMapping(nu, m2, s_min, s_max)`,
   `BreitWignerMapping(M, Gamma, s_min, s_max)`,
   `UniformMapping(s_min, s_max)`), and `.Forward/.Inverse/.Density`.
4. Add `UnitTest`-level and Python tests: forward/inverse round-trip,
   `integral of Density == 1`, `Density(Forward(r))` consistency, for each
   concrete type (closes a documented coverage gap for the non-adaptive
   mappings too).
5. Build/install/ctest/pytest.

**Risk / rollback.** Low. Additive: a base class plus bindings; the perf cost
of one virtual call is negligible next to the Python calls and geometry
already in `Density`. If virtual dispatch is undesirable in a measured hot
path, the fallback is to keep the structs duck-typed and bind only the
concretes (the models then hold a concrete mapping); this loses runtime
swappability for #3 but is otherwise equivalent.

**Done when.** Mappings are a polymorphic family, bound to Python, with
round-trip/normalization tests.

---

## Item 1 -- C1 by construction (Phase C)

The highest-leverage change. Make Sample == Density structural for the
density-bearing variables, not a prose contract validated after the fact.

**Why.** Today a Python model implements `SampleFinalState` and
`FinalStateProbability` as two separately-authored methods (the duck type in
`pyCrossSection.h:33-58` / `pyDecay.h:33-57`; thin C++ adapters in
`PhysicalChannelAdapters.cxx`: decay Sample->SampleFinalState `:66`,
Density->FinalStateProbability `:74`; cross section `:153`/`:161`). The
density-shape agreement between them is maintained BY HAND. This is exactly
where the original 88%-of-variance bug entered (uniform Q^2 draw vs
propagator-peaked reported density). The drift surface, verified in
`resources/processes/DarkNewsTables/VectorPortal.py`:
- `VectorPortalUpscatteringXS`: `_sample_Q2` (rejection from `dsigma/dQ2`)
  vs `FinalStateProbability` reconstructing Q2 from momenta and recomputing
  `diff_xsec_Q2`.
- `VectorPortalOffShellXS`: `_sample_s_pair` (Breit-Wigner inverse-CDF) and
  `_sample_Q2` vs `FinalStateProbability` independently recomputing the BW
  and the dQ2/dcos Jacobian; the chi_prime width is recomputed in both.
- `DarkPhotonDecay`: angular rejection envelope `1 + beta_e^2 cos^2 theta`
  vs density recomputing the same with a separately-derived `beta_e`.

Both peers make this class of bug impossible: Achilles' `Mapper` and
MadEvent's generated routine produce the point and its density from one
expression.

**Target.** Each density-bearing variable in a model is owned by ONE
`Mapping1D` object (from Phase B), stored on the model at construction.
`SampleFinalState` draws via `mapping.Forward(r)`; `FinalStateProbability`
evaluates via `mapping.Density(x)` where `x` is extracted from the record.
The density shape has a single source of truth and cannot drift.

**Design notes.**
- For a PHYSICAL channel the mapping IS the physics shape. Build it once at
  construction from the exact differential rate: e.g. tabulate `dsigma/dQ2`
  into a `TabulatedMapping` (exact), or use `PropagatorMapping` if the shape
  is exactly `1/(Q^2+m^2)^2`. Sampling then needs no rejection loop (faster)
  and the reported density is the same object.
- Multi-variable models (off-shell: s_pair, cos_theta_sub, Q2) hold one
  mapping per factorized variable; the joint density is the product, the
  joint draw is sequential. Order the factorization to match `Measure()`.
- RESIDUAL C1 surface: the mapping guarantees the density-SHAPE consistency
  (the historical bug). It does NOT by itself guarantee that the
  momentum-builder (variable -> momenta in Sample) and the extractor
  (momenta -> variable in FinalStateProbability) are mutual inverses. Add a
  cheap runtime self-check (below) to cover that smaller surface.

**Steps.**
1. (Phase B done.) Refactor `VectorPortalUpscatteringXS` first (single
   variable, the original bug site): construct `self._q2_mapping` once;
   `_sample_Q2` becomes `self._q2_mapping.Forward(random.Uniform(0,1))`;
   the Q2 factor of `FinalStateProbability` becomes
   `self._q2_mapping.Density(Q2_from_record)`. Remove the rejection loop and
   the duplicated `diff_xsec_Q2` density path.
2. Add a generic C1 self-check utility (Python + a C++ analog for
   `ValidateChannelDensities`): after `Sample`, extract the density
   variables from the record, re-evaluate `FinalStateProbability`, and
   assert the drawn point's density equals the reported density to tight
   tolerance (round-trip). Wire it into the existing per-vertex closure
   tests.
3. Migrate `VectorPortalOffShellXS` (s_pair via `BreitWignerMapping` or
   `TabulatedMapping`; Q2 via `PropagatorMapping`), then `DarkPhotonDecay`
   (angular shape via a small tabulated/analytic mapping). The remaining
   isotropic decays already satisfy C1 trivially (constant density) -- leave
   them, but route through `UniformMapping` for uniformity if cheap.
4. Optional interface hardening: extend the model duck type with an
   introspection method `Mappings()` returning the per-variable mappings, so
   `PhysicalCrossSectionChannel`/`PhysicalDecayChannel` (and tooling) can
   assert a model that declares a `Measure` actually backs each density
   variable with a mapping. Keeps C1 checkable at construction, not just at
   closure time.
5. Build/install/ctest/pytest after EACH model migration (one model per
   commit), gating on closure.

**Tests.** Strengthen the per-vertex closure tests
(`test_dutta_kim_chain.py`) with the round-trip self-check from step 2.
Add a targeted regression: deliberately perturb a mapping parameter and
assert closure breaks (proves the test has teeth). Keep
`test_dutta_kim_models.py` topology/measure assertions.

**Risk / rollback.** Med. Touches the physics models; a wrong mapping build
shows immediately as a broken closure (the oracle). Mitigate by migrating
one model per commit, upscatter first (best-understood, single variable).
Rollback per model is independent.

**Done when.** Each density-bearing variable in the VectorPortal models is
backed by one shared `Mapping1D`; the round-trip self-check is in the
closure tests; the perturbation test confirms teeth; full suite green.

---

## Item 3 -- Adaptive intra-channel sampling (Phase D)

Catch up to the field's second adaptive layer; retire the manual-mapping tax.

**Why.** MadGraph and Achilles both run TWO adaptive layers: inter-channel
KP alpha AND intra-channel VEGAS grids; the frontier (i-flow / MadNIS, same
Achilles authors) replaces VEGAS with normalizing flows. SIREN adapts ONLY
the inter-channel alpha. Its per-channel densities are fixed analytic maps,
and the recent `TabulatedMapping` / `PowerLawMapping` / `PropagatorMapping`
work is, viewed from the field, hand-rolled per-variable reinvention of what
1-D VEGAS does automatically. Every new model owes a hand-tuned mapping per
peaked variable, and a forgotten one silently re-creates a variance problem.

**Target.** An `AdaptiveMapping : Mapping1D` -- a VEGAS-style 1-D adaptive
grid that refines from sampled weights -- usable anywhere a `Mapping1D` is,
and refined inside the existing optimization loop alongside the KP alpha
update.

**Design.**
- State: `bin_edges` (n+1 monotone over [x_min, x_max]) and
  `cumulative_weight` (the running second-moment-per-bin estimator, exactly
  the VEGAS damping target). `Forward`/`Inverse`/`Density` are piecewise on
  the grid (same math as `TabulatedMapping`, which already does piecewise
  CDF inversion -- reuse it).
- Adaptive hooks (the Phase B base already declares them as no-ops):
  `Accumulate(x, weight)` bins the contribution `f^2/g`-style; `Refine()`
  rebins to equalize per-bin contribution (standard VEGAS rebinning, with
  the usual smoothing/damping factor).
- Integration: `optimize_chain_weights` already samples batches and iterates
  (`python/optimize.py:123+`). In each iteration, after computing total
  weights, call `Accumulate` on the adaptive mappings for the sampled
  points, then `Refine` between iterations -- the SAME warm-up loop that
  does the KP alpha update. This is why #3 depends on #5 being settled
  (one loop, two adaptations, composed like the peers).
- Freeze the grid for production (read the grid back via a getter; persist
  alongside the tuned alphas).

**Steps.**
1. Implement `AdaptiveMapping` reusing `TabulatedMapping`'s piecewise CDF
   inversion for the static part; add binning + rebinning for the adaptive
   part. Unit-test: a known peaked target (e.g. `1/(x+a)^2`) -- assert the
   refined grid concentrates bins at the peak and `Var(1/g)` drops across
   refinements toward the analytic optimum.
2. Add an `Adaptive` option to the channel mass/Q2 mode enums
   (`DetectorDirected3BodyChannel`, `DetectorDirectedScatteringChannel`) and
   store adaptive state as channel members; bind to Python.
3. Extend `optimize_chain_weights` to drive `Accumulate`/`Refine` on any
   adaptive mappings found in the chain's channels, in the same iteration
   loop as the alpha update.
4. Demonstrate the tax reduction: replace one hand-tuned mapping in the
   example (e.g. the pion `s_X` PowerLaw) with `AdaptiveMapping` and show it
   reaches comparable-or-better `Var(log W)` with no manual shape tuning.
5. Build/install/ctest/pytest; A/B the demonstrated vertex (C0.2 gate).

**Risk / rollback.** Med. New capability with a new control-flow hook in the
optimizer. Contained: adaptive mappings are opt-in per channel; if disabled
the chain is byte-for-byte the Phase C result. A poorly converged grid only
costs variance (unbiased), and closure still gates correctness.

**Done when.** `AdaptiveMapping` exists, refines in the optimizer loop, and
the demonstrated vertex matches/beats its hand-tuned mapping without manual
shape tuning.

---

## Item 6 -- Narrow the measure-conversion surface (Phase E)

Final cleanup. NOT a removal of the measure framework -- that abstraction is
downstream-justified by the detector-directed feature (lab-frame geometric
channels mixed with rest-frame physical channels genuinely need boost
reconciliation, which neither peer faces). The goal is to shrink the
EXERCISED conversion set and make the unsupported set loud, moving toward
the peers' common-base-parameterization ideal where practical.

**Why.** Both peers avoid an N x N conversion matrix entirely (Achilles: one
[0,1]^n hypercube base; MadEvent: per-node Jacobians composed from the
diagram). SIREN reifies named measures and converts pairwise
(`PhaseSpaceChannel.cxx` ConvertDensity + `CommonMeasure` :380). After #4
collapses redundant channels, most mixtures become single-measure and
ConvertDensity is rarely called (it only fires when `ch_meas !=
common_meas`, `:353`). The main remaining cross-measure case is the pion
vertex's Recursive2Body cross-factorization (spec=2 directed vs spec=0
physical), which round-trips through Dalitz (`:240-264`).

**Steps.**
1. Use #2's audit (the supported-set test) as the authoritative list of
   conversions the chains actually need. Anything outside it now throws (#2).
2. Align the pion-decay factorization (gap P6): give the physical channel
   the same factorization indices as the directed channel so no
   cross-factorization conversion is needed -- eliminating the
   Recursive2Body<->Dalitz round trip at that vertex.
3. Document the SUPPORTED conversion set as a short explicit table in
   `phase_space_biasing_design.md` (and as the throw's allowlist), so the
   surface is visible and bounded rather than implied by code.
4. Re-run the full suite; confirm the exercised conversion count dropped and
   no throw fires on a legitimate path.

**Risk / rollback.** Low. Cleanup gated by the same closure oracle.

**Done when.** The pion round-trip is gone, the supported set is a documented
allowlist, and the exercised-conversion count is measured and minimal.

---

## Cross-cutting validation (every phase)

1. `cmake --build . && cmake --install .` from `build/` (C0.4).
2. `ctest --output-on-failure` -- expect 39/42 (C0.5 baseline); no NEW
   failures.
3. `python -m pytest tests/python/ -q` -- fully green.
4. For variance/throughput items (#3, #4, #5): A/B total weight per POT
   (must agree within statistics -- C0.2), `Var(log W)`, effective sampling
   fraction, density-eval count; record in the relevant companion doc.
5. For correctness items (#1, #2, #6): the per-vertex closure tests and the
   normalization/closure suite are the oracle; a perturbation test proves
   the relevant test has teeth.

## Definition of done (overall)

- C1 is structural for the VectorPortal density variables (#1): one shared
  mapping per variable, round-trip self-checked, perturbation-tested.
- ConvertDensity cannot silently mis-report a density (#2); the supported
  conversion set is a documented allowlist (#6).
- An adaptive mapping exists and refines in the optimizer loop (#3),
  retiring at least one hand-tuned mapping with no accuracy loss.
- Redundant geometric channels are collapsed (#4); the KP update rule is
  verified against the canonical form with recorded A/B evidence (#5).
- The per-vertex weight breakdown is regression-guarded and promoted (#7).
- `phase_space_biasing_design.md` gap list (G6, R2, P5, P6, and the KP and
  adaptive notes) is updated to reflect what landed.

---

## Implementation status (changelog)

Correctness-first scope (Phases A + B + C). All items below verified against
the full Python suite (475 passed, 12 network-skipped) and the C++
phase-space tests; the only C++ failures are the documented pre-existing ones
(DISFromSpline / DipoleFromTable / HNLDISFromSpline data files, and the flaky
fixed-seed `Cone.SampleDistribution*` histogram bound), none touched by this
work.

- **A1 (G6) -- DONE.** `ConvertDensity` (`PhaseSpaceChannel.cxx`) now throws
  (`ThrowUnconvertible`) on unsupported measure pairs and on
  supported-but-missing-input cases (bjorken_y, pair invariant mass, <2
  decay masses), instead of silently returning the input. Legitimate
  degenerate-physics no-ops kept (parent-at-rest Rest==Lab, etc.). New C++
  tests `ConvertDensity.ThrowsOnUnsupportedMeasurePair` /
  `SameMeasureMixtureDoesNotConvertOrThrow`; existing `AutoConversion` tests
  and the Python closure suite still pass (no false throw on any exercised
  conversion).
- **A3 (KP rule) -- DONE.** `python/optimize.py` factors the update into
  `_kp_update` with a selectable `update_rule` ("sqrt_W" default vs canonical
  "alpha_sqrt_W"), used by both optimizers. New `tests/python/test_optimizer_kp.py`
  pins the fixed-point difference (canonical -> W_i=const; memoryless ->
  alpha ~ c_i^(1/3)) and unbiasedness of both. Default unchanged pending the
  real-chain A/B (folded into A2).
- **A4 (P5) -- DONE.** `test_breakdown_sum_invariant` guards that
  `Weighter.breakdown()` reconstructs the C++ EventWeight up to one global
  constant (found empirically to be exactly ln(N_gen), the 1/N MC
  normalization) -- i.e. it captures every event-varying factor, which is
  what makes the per-vertex variance attribution valid.
- **B (Mapping1D) -- DONE.** `InvariantMassMapping.h` gains an abstract
  `Mapping1D` base (Forward/Inverse/Density + no-op Accumulate/Refine hooks
  reserved for the deferred adaptive map); the five concrete structs derive
  from it (zero call-site changes). Bound to Python in
  `pybindings/injection.cxx` (`class_<Concrete, shared_ptr<Concrete>,
  Mapping1D>`; models consume, no trampoline). New `tests/python/test_mappings.py`
  covers polymorphism, round-trip, normalization, and the no-op hooks.
- **C (C1 by construction) -- headline DONE.** The verified off-shell
  `s_pair` Sample!=Density bug is fixed: `VectorPortalOffShellXS` now routes
  both `_sample_s_pair` (Forward) and the `FinalStateProbability` BW factor
  (Density) through one shared `BreitWignerMapping` via `_s_pair_mapping`, so
  the reported density is the [s_min,s_max]-renormalized form the sampler
  actually draws from (the old code reported a full-line 1/pi BW that
  integrated to (hi-lo)/pi < 1). New regression test
  `test_offshell_s_pair_sample_matches_density`. **Scope finding:** the
  upscatter Q2 and DarkPhotonDecay angular variables are already
  C1-compliant (sampler and density share one physics function with exact
  integral normalization), so they were intentionally NOT migrated to
  `TabulatedMapping` (which would add tabulation error for no correctness
  gain). The residual off-shell C1 surface (the hand-coded `dQ2/dcos`
  Jacobian and the momentum builder/extractor inverse) is exercised by the
  existing per-vertex closure tests; a dedicated round-trip self-check is an
  optional follow-up.
- **A2 (R2 collapse) -- DONE.** `build_geometric_targets` gains a
  `fiducial_only` switch (physical + one fiducial-directed channel per biased
  vertex instead of 1 + 13).  On-shell A/B (monoenergetic 2 GeV pion, ~700
  events, light optimization): `all_13 sqrt_W` -> 14 primary channels,
  sum_w=5.54e-19, Var(logW)=0.151, eff 87.7%; `fiducial_only sqrt_W` -> 2
  channels, sum_w=5.45e-19, Var(logW)=0.144, eff 87.7%.  The rate agrees
  within ~1.6% (MC error ~1.4% at this N) and Var(log W)/eff are unchanged,
  so collapsing 13 -> 1 is accuracy-neutral at ~7x fewer primary-vertex
  density evaluations -- confirming the optimizer cannot distinguish the
  near-degenerate targets (the Kleiss-Pittau degeneracy premise).  The switch
  is opt-in (default still all_13); flipping production runs to fiducial_only
  is a one-line change at the `build_geometric_targets` call site.
- **A3 real-chain A/B -- DONE.** On the same chain the canonical
  `alpha_sqrt_W` (Var(logW)=0.159, eff 87.2%) and the memoryless `sqrt_W`
  (0.144, 87.7%) are equivalent within noise -- expected, since the collapsed
  mixture's channels are near-degenerate and both rules drive to ~uniform.
  Default kept as `sqrt_W`; the canonical rule stays available for chains
  with genuinely distinct, variance-imbalanced channels.

**Remaining (deferred follow-ups, Phases D and E):**
- **Phase D (AdaptiveMapping)** and **Phase E (narrow the conversion
  surface, pion P6)** -- the Mapping1D base + its no-op Accumulate/Refine
  hooks are already in place for D.
```
