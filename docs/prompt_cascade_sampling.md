# Supported targeting for prompt on-shell cascades

This development branch supplies three complementary changes. The two-body
directed-decay branch probability is now proportional to the inverse lab/rest
solid-angle Jacobian. Scattering retains its uniform branch selection. Existing
weighted samples must be regenerated; the change deliberately has no legacy
branch mode. Do not evaluate an old generation density using the new code.

`RestFrameEnvelope2BodyChannel(target, daughter_index)` samples a conservative
rest-frame polar interval union and azimuth envelope of the target's bounding
cone. It includes both inverse-boost branches and both boost-axis poles. An
empty or numerically unresolved envelope uses normalized isotropic sampling,
with `DirectingActive` false. Exact detector misses are ordinary zero scores.
The public facade is `channels.rest_frame_envelope(daughter, target)`.

`OnShellCascadeChannel(target, pair_mass, kappa=0,
orientation_weights=(.1,.45,.45), first_daughter_probability=.5, volume=-1)` handles
prompt `P -> spectator + R -> spectator + d1 + d2`. Its first implementation
requires equal pair daughter masses and open two-body phase space at both steps.
The final-state order is spectator, d1, d2. The orientation weights select
isotropic, envelope, and inverse-Jacobian volume proposals. A strictly positive
isotropic component is required. The internal cosine proposal is
`(1 + kappa*c*c)/(2*(1 + kappa/3))`, with `kappa >= -1`; its meaning is a
proposal, not an amplitude assumed by the framework. Model physics must supply
its own normalized angular law, widths and branching fractions.

The channel reconstructs the entire configuration and its intermediate pair
four-momentum. Both targeting densities enter every weight:

```
q(x) = p_internal(x) * 4*pi * (rho*Q1(x) + (1-rho)*Q2(x)).
```

`InternalDensity` returns the normalized internal proposal with isotropic global
orientation. It is not an independent physical oracle. There is no `ds_pair`
factor: `Measure.OnShellCascade(pair_mass, spectator=0, pair_first=1,
pair_second=2)` is a four-dimensional angular measure at a fixed mass. The mass
is part of measure identity. Different masses, factorization orderings and
continuous three-body measures are not pointwise convertible. The measure's
archive version is 1 for the added constraint field; version-0 measures still
load. This version change is unrelated to the directed branch repair.

The public facade is `channels.on_shell_cascade(...)`. Native channels and
mixtures support archive/pickle round trips. Generic Python decay trampolines
are registered in the shared serialization header so physical/proposal mixtures
can preserve model state. Physical Python models must implement their own
pickle reconstruction when their native base does not provide it.

`interaction_parameters` retains `cascade_pair_energy`, `cascade_pair_px`,
`cascade_pair_py`, `cascade_pair_pz`, `cascade_pair_mass`,
`cascade_internal_kappa`, and `cascade_target_daughter`. The pair's parent is
the primary at this vertex; its daughters are secondaries 1 and 2. These fields
describe provenance, not an independently weighted propagating particle.
The polynomial internal correlation is preserved; a consumer requiring a
general spin-density matrix or a displaced mediator must use another suitable
physical representation. This implementation does not supply either.

## Source allocation

`siren.source_importance.SourceImportanceTable` freezes nonnegative physical
source weights and pilot mean scores. A baseline fraction in `(0,1]` preserves
every populated row, including chance-zero pilot cells. `sample(rng, size)`
returns row indices and the exact normalized physical/proposal ratio. Normalize
scores by every attempted source draw. Preserve source exposure separately.

For production, use `table.to_distribution(keys, rows, metadata=definition)`.
This builds the existing native `PrimaryExternalDistribution` with physical
`weight` values and independent row-sampling probabilities. Declare that same
object on the injection and physical sides of a `Vertex`. The native weighter
then owns the source correction and absolute normalization; the driver must
not multiply a second `p/q` factor. Row counts, normalized physical weights,
finite values and the frozen metadata fingerprint are checked. The caller
must preserve source ordering and identify its content hash in metadata.
Equal-weight rows cannot reveal reordering from their weights alone.

The stored definition includes source identity, physics, geometry and scoring.
`load(path, metadata)` rejects a changed definition. Include masses, relevant
couplings, decay law, lifetime, row mapping, analyzed daughters, material,
response, timing and cuts as applicable. A sampler-only change need not
invalidate the physical pilot expectation. Freeze the pilot before production
and use an independent random stream.

## Native production and persistence

`injection.PhaseSpaceDecay(signature, masses, partial_width, total_width,
physical_channel)` is a native `Decay` for an explicitly normalized physical
phase-space law. It owns its masses, widths and physical channel. Supply the
biased proposal separately through `Vertex(kinematics=...)`; substituting that
biased law as the physical channel changes the model. The adapter lives in
injection because injection already depends on interactions and owns
PhaseSpaceChannel. Placing this dependency in interactions would create a cycle.

`TotalDecayWidthAllFinalStates` is the sum of widths owned by this model,
which for this adapter is its partial width. `ParentDecayWidth` declares the
whole-parent width, including omitted branches. Other `Decay` implementations
opt in by overriding that method; its default zero preserves their additive
contract. Python models reach it through both trampolines, `siren.DecayModel`
and `decay_model_base(base=siren.interactions.DarkNewsDecay)`, and the authoring
base rejects near-miss spellings of the name. A collection counts a declared parent width once, rejects inconsistent
declarations or owned widths exceeding it, and uses it for the lifetime.
The physical weighter includes omitted branches in both Fixed and Propagated
probabilities. Forced generation remains conditional on modeled channels.
For decay-only processes this supplies sum(modeled partial widths)/parent width;
with material competition it supplies modeled rate/total rate. Do not add an
external branching normalization for the same decay. Width declarations do not
change the normalized final-state law.

The adapter supports native save/load and pickle. A standard
`Simulation(...).run(on_failure="raise", on_shortfall="raise")` returns weights
normalized over all injection attempts, including geometric misses. `Results`
and `Weighter.explain` use the same native weight path.

Pure `expand.child(...)` and `expand.depth_below(...)` declarations compile to
native `SecondaryExpansion` rules. Injector archive version 3 persists them;
versions 0–2 remain readable. A terminal vertex can declare
`expand=(expand.depth_below(0),)`. Arbitrary `continue_if` predicates retain the
Python callback path and remain ineligible for native archives. Native event
trees and weighters reload without importing a Python physics callback.

For an inclusive expected number of daughter interactions, generate separate
single-daughter chains and sum their estimates. Expanding two independently
interacting daughters in one tree multiplies their probabilities and represents
a coincidence. An at-least-one-event probability needs a separately specified
observable including the overlap; it is not generally the sum of rates.

## Validation boundary

Tests cover independent forward boosts, branch balance, angular coverage,
normalization, on-shell constraints, physical moments, both-daughter densities,
constrained-measure rejection, and persistence. The workspace evidence is
`artifacts/idm-directed-sampling-production-20260923`. Its source table and
cylinder scores are controlled fixtures; physical CCM argon rates and production
timing/selection remain unqualified.

## Review repairs and numerical domain (2026-09-24)

The three directed decay channel archive versions are now 1. Their version-0
payloads fail explicitly; regenerate events because their generation density
changed. OnShellCascade, RestFrameEnvelope2Body, PhaseSpaceDecay and
SecondaryExpansion also use version 1 and reject their earlier development
payloads. The latter used ambiguous wildcard semantics and the decay adapter
previously required an external branching factor. New archives preserve mixture
weights bitwise without repeated normalization. Measure version-0 reading is
unchanged: it does not authorize reading old directed samplers.

The new cascade/envelope channels use the declared parent mass and three-momentum
for a stable known-mass boost. Input energy must agree within 2e-5 relative,
allowing ordinary tabulated float32 inputs. The daughters conserve the on-shell
four-momentum `(hypot(m, |p|), p)`. The event record keeps the parent's input
energy, because the parent belongs to the upstream vertex and that vertex's
densities (for example a source table's `emin`) read it. For a tabulated row the
daughters' total energy can therefore differ from the recorded parent energy by
up to the 2e-5 tolerance. Proposal and physical sampling behave the same way,
and Density uses the same convention. Larger discrepancies and closed phase
space are `InjectionFailure(KinematicallyForbidden)`; they enter the attempt
ledger. This is an explicit tabulation policy, not acceptance of arbitrary
off-shell physics.

A physical model paired with these channels must rebuild the parent rest frame
from the same on-shell parent, not from the recorded energy. The difference is
not small at high boost: a float32 pi0 row at gamma 1000 has an invariant mass
sqrt(E^2 - p^2) about 6% away from m, and at gamma 3000 the rounded row can be
massless, because an energy error dE shifts E^2 - p^2 by 2 E dE. `PhaseSpaceDecay` and the BSM-beam `PromptPionCascade`
follow this convention. A model that reads `primary_momentum[0]` directly gets
a wrong rest-frame angle for such rows. For parents that are exactly on shell
the two conventions agree.

The pair-mass check in Density allows for the rounding of the boost. A lab
daughter built from a rest-frame vector carries an absolute error of order
eps * gamma * (E* + p*). For a pair emitted backward from a fast parent this is
far larger than eps times the pair's own lab energy. The tolerance therefore
has a term eps * gamma * E_A * (E + |p|) of the lab pair, in addition to the
rounding of the invariant itself. Daughter energies rebuilt in the parent frame
use the matching error, amplified by gamma on inversion. With this, a pi0 of
40-135 GeV decaying to a 1-3 MeV A' that is aimed at a target 3-6 degrees off
the beam axis keeps a positive density for its own samples.

A nonanalytic geometry requires `volume=` when its volume proposal has nonzero
weight. The constructor resolves that volume before generation; the native
channel, Python facade and BSM builder all accept it. An envelope-only mixture
does not require a volume. "Validates" means: for a box, cylinder, sphere,
ellipsoid, elliptical tube, cone, torus, Trd, parallelepiped, polycone or generic
polycone with a simple R-Z profile (with their cuts and azimuthal segments), or a
closed, outward-oriented triangle mesh of one connected piece (divergence theorem,
used only when the ray estimate below resolves the mesh and agrees with it, since a
closed surface need not bound one solid), a supplied volume must agree with the
exact volume to 1e-4, and the exact value is then used; for any other solid (Boolean composites,
multi-piece meshes and the rest) it must not exceed the bounding box and must
agree with a chord-integration estimate of the solid to the larger of 2% and five
standard errors. The estimate casts rays along 48 fixed directions spread over a
hemisphere, each on a jittered 48 x 48 grid over the smallest rectangle around
the bounding box's projection, and adds the lengths of the chords inside the
solid; the volume is the mean over directions and its standard error their
scatter, so a direction running along a thin wall (which it misses, or meets in
a rare long chord) widens the error instead of biasing the result. The supplied
volume cannot be checked, and is refused, when more than 4 of the directions
have fewer than 20 rays meeting the solid (a speck or an empty solid); when a
crossing distance is not finite, or the crossings of more than a few rays do not
alternate entering and exiting (an open, doubled, nested or self-intersecting
surface, whose inside test and chords disagree); when the solid's inside test,
which the channel samples points with, disagrees with its crossings; when more
than 1% of the chords are shorter than 5e-9 m, where `BooleanGeometry` merges
crossings and loses material; when parts of the solid could hold more than 1% of
the volume that the rays missed (a small core inside a thin shell, which every
ray can miss while the shell looks resolved): each primitive with an exact volume
inside the solid's box is estimated on its own, and its unexplained shortfall
counts; any other part (a connected mesh piece, an intersection or subtraction
node) met by fewer than 10 rays counts with its bounding box; or when the standard
error exceeds 5%. The rays cannot see material that no ray meets inside a part
they do meet, such as a small boss on a thin plate selected by intersecting or
subtracting broad operands, so for Boolean solids the check is a safeguard, not a
guarantee. Any volume-mode target is also refused when the layer the point
sampler's inside test treats as outside (the last 1e-9 m before each exit along
+z) holds more than 0.5% of its volume, as for walls or plates a few nm thick.
The grid is laid out in the solid's own frame, so a rotated slender solid is
covered as tightly as an unrotated one. Such a
target needs an analytic shape or an envelope-only proposal. For well-resolved
solids the tolerance is about 2%; for thin curved shells it can be 5% or more.
The check casts about 110,000 rays plus inside tests on a sub-grid (tens of
milliseconds for simple solids; about 1.4 s for a union of 216 spheres, whose
every intersection query visits all of them). `SetVolume` on the directed
channels applies the same check, and loading an archive refuses a stored volume
these rules reject instead of reading it as a different density.
A nonanalytic volume wrong by less than that tolerance cannot be detected, so it
should still come from an independent calculation: the volume sets the proposal density while the
points come from the solid, and a wrong value biases every weight. Source metadata use one canonical JSON
representation (including NumPy scalars and converted keys) and exact save paths.

Independent quadratic-root envelope areas constrain normalization to 2e-9
relative. A bounded-weight cascade check has an absolute 1.2% tolerance and
explicitly rejects a 6% density scaling. These supplement the earlier broad
stochastic stress tests, whose tolerances alone did not establish percent-level
accuracy. Domain-specific high-boost and persistence regressions are separate.

### Configuration and caller state

The current cascade channel supports spectator index 0 and pair indices 1, 2,
with **exactly equal** pair masses. It does not round unequal masses into the
same model. `PhaseSpaceDecay` rejects an incompatible mass specification at
construction. Direct sampling of a mismatched event remains a classified
`KinematicallyForbidden` failure. The Python `on_shell_cascade` facade checks
matching models' constrained measure, pair mass and daughter masses when it
compiles a vertex. Its optional `spectator=` and `pair=` particle names make the
intended ordering explicit and reject unsupported permutations. Equal masses
alone cannot establish which pair the model physically intends.

Construct constrained measures through `Measure.OnShellCascade(...)`.
`pair_mass` is read-only in Python, and measures now pickle through their
validated native archive. This permits a Python model to retain its measure.
`PhaseSpaceDecay.SampleFinalState` preserves parameters added to the mutable
distribution record before sampling; the channel can still update its own keys.

A built injector's `stopping_condition = None` resets to the native default,
which stops all secondaries. Native rule and callback setters replace each other.
Clearing rules cannot reactivate a prior callback. Version-3 injector archives
own the expansion policy, including the default when no rules are stored;
loading them replaces an existing callback. Older archives retain their legacy
caller-policy behavior. The Python stopping-condition property follows direct
engine changes; public Python callbacks retain their existing pickle path.

## Follow-up review fixes (2026-09-24)

**One decay frame everywhere.** `Isotropic2BodyChannel`, `DetectorDirected2BodyChannel`,
`DetectorDirectedAngularSectorChannel` and the rest/lab conversion in `ConvertDensity`
now build the parent frame from the declared mass and three-momentum, as the cascade
and envelope channels do, and reject a parent more than 2e-5 off its mass shell
(`InjectionFailure(KinematicallyForbidden)` in Sample, density 0).
`DetectorDirected3BodyChannel` does the same for decays and, with the beam's declared
mass (zero allowed), for 2->3 scattering, and its inner boosts use the known masses of
the complementary system and of the parent. For a parent whose
energy already equals `hypot(m, |p|)` nothing changes. A float32 row at gamma 3000 used
to give NaN momenta in the isotropic sampler, and the directed density and the rest/lab
Jacobian differed by up to a factor 2 between the two frames. The three-body helpers in
`siren.three_body` still read the recorded energy.

**Canonical tables (opt-in).** `siren.dist.on_shell_rows(keys, rows)` rebuilds each row's
energy as `hypot(m, |p|)`, rejects rows more than 2e-5 off their shell (or whose rebuilt
energy is not finite) and tables that repeat a kinematic column, and keeps the input
energy in an `E_table` column, which the distribution stores in each record's interaction
parameters. Another name may be chosen with `keep_input_energy`, but not one the
distribution interprets (`weight`, `t0`, a kinematic or position column, or its own
`PrimaryExternalDistribution_*` bookkeeping). Use it for rounded tables before constructing `PrimaryExternalDistribution`
when physical models or other code read `primary_momentum[0]` directly: then every
consumer, including the `emin` cut, sees the same parent.

**Mass check limits.** At high boost the pair-mass check can only resolve masses to the
precision a lab-frame double carries: at gamma 3000 a 1 MeV pair mass shifted by 0.1%
can pass. It rejects other constraints, not near-identical masses; do not use it as a
mass discriminator.

**Other contracts.** A parent at rest with cross sections in its collection raises
`ConfigurationError` (material is located along the direction of motion) instead of
crashing; decays at rest are unchanged. `Injector.SaveInjector(path)` refuses a live
stopping callback, which the archive cannot hold; `SaveInjector(path, True)` writes it
anyway, and the legacy controller does so with a `RuntimeWarning`. Pickling an
`Injector` whose callback was set on the raw engine raises `NotSerializableError`
with that explanation. `PhaseSpaceMeasure` fields are read-only in Python; build
measures with the factories. `SourceImportanceTable.save`/`load` also accept binary
file objects. Point-source, column-depth and range position distributions now give
their decay-length query the parent's momentum; before, every decay of a decaying
primary was placed at the start of the path.

**Lab boost precision.** The older two-body channels (and the isotropic regime of the
shared directed step) boost with the declared parent mass instead of rebuilding it from
E^2 - |p|^2, so their frame is exactly the one their densities use. A double-precision
lab four-vector still fixes rest-frame quantities only to about eps*gamma^2 (3e-10 at
gamma 1000, 3e-9 at gamma 3000), whatever the boost; code that rebuilds rest-frame
energies from lab daughters must allow for that, as the cascade tolerances do.
