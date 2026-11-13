# Four-momentum values and records

`SIREN/math/FourMomentum.h` provides `siren::math::FourMomentum`, a value
type backed by rk. Components are `(E, px, py, pz)` in GeV with metric
`(+,-,-,-)`.

```cpp
#include "SIREN/dataclasses/RecordKinematics.h"
using siren::math::FourMomentum;
using namespace siren::dataclasses;

auto recorded = PrimaryMomentum(record);
auto parent = FourMomentum::FromMass(record.primary_mass, recorded.Momentum());
auto daughter = FourMomentum::FromMass(daughter_mass, rest_momentum);
SetSecondaryMomentum(record, 0, daughter.BoostToLab(parent));
```

- `FromComponents(array)` preserves all four supplied doubles, including
  signed energy and spacelike transfers. Addition and subtraction return raw
  components. `MassSquared()` is signed and subject to cancellation in
  `E*E - p*p`. `Mass()` rejects spacelike and nonfinite raw inputs, including
  a collinear sum whose invariant rounds negative.
- `FromMass(mass, three_momentum, negative_energy=false)` constructs an
  on-shell particle. `Mass()` returns the retained mass and `HasKnownMass()`
  reports its presence. An array round trip retains only the components. Use
  `FromMass` when you have an independently known mass.
- `BoostToLab(frame)` and `BoostToRest(frame)` require retained masses for
  both inputs. The frame needs positive mass and energy. The particle may be
  massless or have either energy sign. Unsupported physical inputs or nonfinite
  results raise `std::domain_error`, even for a zero boost. Finite inputs whose
  intermediate squares overflow fall outside rk's supported range.
- `InvariantMass(a, b)` uses `rk::invMass` and uses retained masses where
  available. Inputs need real particle masses and a finite nonnegative pair
  invariant. The cross-dot product can still lose precision for highly collinear
  high-energy pairs. Inverse boosts amplify component rounding at large gamma.

`RecordKinematics.h` provides `PrimaryMomentum`, `SecondaryMomentum`,
`SetPrimaryMomentum`, `SetSecondaryMomentum` and `InteractionVertex`. Reads
return independent values and setters update momentum arrays. Secondary indexing
raises `std::out_of_range` on an invalid index. `HasSecondaryStorage` and
`RequireSecondaryStorage` check the lengths of the type, mass and momentum
arrays. Records store the existing component arrays and mass metadata.

The isotropic two-body channel accepts threshold equality and a relative
parent-energy discrepancy up to `2e-5` for rounded beam rows. The declared mass
and three-momentum define the on-shell decay frame. The input record retains its
energy for upstream densities. The first daughter is boosted with its known
mass, and the second is computed by momentum conservation. Invalid parents or
daughter masses raise `InjectionFailure`. An intermediate overflow can raise
`std::domain_error`.
