# Four-momentum values and records

`SIREN/math/FourMomentum.h` supplies `siren::math::FourMomentum`, a small
value around rk's `P4`. rk implements arithmetic, particle boosts and pair
invariant mass. SIREN supplies explicit construction, checked particle/frame
operations and record conversion. There is no new persisted representation.
Components are `(E, px, py, pz)` in GeV, with metric `(+,-,-,-)`.

```cpp
#include "SIREN/dataclasses/RecordKinematics.h"
using siren::math::FourMomentum;
using namespace siren::dataclasses;

auto recorded = PrimaryMomentum(record);  // Copy; preserves the input energy.
auto parent = FourMomentum::FromMass(record.primary_mass, recorded.Momentum());
auto daughter = FourMomentum::FromMass(daughter_mass, rest_momentum);
SetSecondaryMomentum(record, 0, daughter.BoostToLab(parent));
```

`FromComponents(array)` preserves all four supplied doubles, including signed
energy and spacelike momentum transfers. It does not validate or project a
mass shell. Addition and subtraction also preserve their component results.
`MassSquared()` is the signed component invariant from rk, including negative
values. It still suffers cancellation in `E*E - p*p`; replacing a spelling of
that expression does not improve precision.

`FromMass(mass, three_momentum, negative_energy=false)` constructs an on-shell
particle and retains the nonnegative mass. `Mass()` returns that retained mass;
for a raw component value it uses rk's reconstructed invariant and rejects
spacelike/nonfinite inputs. `HasKnownMass()` distinguishes a retained mass from
a raw component value. Converting to an array and back discards the retained
mass: callers with independently known mass must use `FromMass` explicitly.

`BoostToLab(frame)` and `BoostToRest(frame)` return values. The frame must have
positive mass and positive energy. Particle energy may have either sign,
including massless particles; a massless particle cannot define a rest frame.
Spacelike particles/frames are rejected with `std::domain_error`, even for a
zero boost, instead of reaching rk assertions or silently clamping the mass.
Nonfinite physical inputs or results also raise `std::domain_error`. rk's
intermediate squared magnitudes limit the representable range; finite doubles
whose squares overflow are not supported as physical momenta.

`InvariantMass(a, b)` uses `rk::invMass` with each particle's retained mass
where available. Its inputs must have real particle masses, and the pair must
have a finite nonnegative invariant. It avoids reconstructing *individual*
known masses but rk's cross-dot product can still cancel for highly collinear,
high-energy pairs. It is not a general high-precision invariant solver.

Forward boosts retain mass but round components. Inverse boosts amplify input
component error by the frame's gamma; at large gamma, do not expect bitwise
component round trips or recover tiny masses from arrays. Use an independent
reference and bounds scaled to the boost when validating such calculations.

`RecordKinematics.h` provides `PrimaryMomentum`, `SecondaryMomentum`, explicit
`SetPrimaryMomentum`/`SetSecondaryMomentum`, and `InteractionVertex`. Reads
return independent values; writes change only momentum arrays, not declared
masses or other metadata. Secondary indexing is checked with `std::out_of_range`.
`HasSecondaryStorage`/`RequireSecondaryStorage` check the type, mass and momentum
array lengths for a channel. The latter retains the contextual runtime error.
Archive fields, versions, existing array access and record ownership are unchanged.

The isotropic two-body channel validates its decay domain before construction.
Its existing `2e-5` relative parent-energy tolerance allows rounded beam rows;
the declared mass and three-momentum define the decay frame. Daughters conserve
that on-shell frame while the supplied parent record retains its original
energy for upstream densities. Threshold equality remains allowed, and invalid
parents or daughter masses retain the channel's `InjectionFailure` contract.
The first daughter now carries its known mass through the boost instead of
reconstructing it by subtraction; the second is stored by momentum conservation.
