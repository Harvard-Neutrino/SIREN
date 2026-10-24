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

The common boost retains nonzero boosts below the former adapter's
`|p|² < 1e-30` cutoff. Parent energy follows rk's square-root calculation,
which can differ by a few ulp from the former `hypot` calculation. The RNG
draw sequence is unchanged.

## Python authoring

The same value is `siren.math.FourMomentum`. Existing `siren.utilities.rk`
exports and array fields remain available. Python checks raise `ValueError`
for invalid particle/frame domains and `IndexError` for secondary bounds.

```python
import pickle
from siren import dataclasses, math

record = dataclasses.InteractionRecord()
record.primary_momentum = [1.0, 2.0, 3.0, 4.0]
raw = record.get_primary_momentum()
assert raw.mass_squared() == -28.0
record.set_primary_momentum(raw + raw)
assert raw.components() == [1.0, 2.0, 3.0, 4.0]

parent = math.FourMomentum.from_mass(2.0, math.Vector3D(0, 0, 2e9))
daughter = math.FourMomentum.from_mass(0.001, math.Vector3D(0, 0, 0))
lab = daughter.boost_to_lab(parent)
assert lab.mass() == 0.001
assert pickle.loads(pickle.dumps(lab)).mass() == 0.001
```

Methods are `from_components`, `from_mass`, `components`, `energy`, `momentum`,
`mass`, `mass_squared`, `has_known_mass`, `boost_to_lab`, `boost_to_rest` and
`invariant_mass(a, b)`, with `+`/`-` for component arithmetic. Record methods are
`get_primary_momentum`, `set_primary_momentum`, `get_secondary_momentum(index)`
and `set_secondary_momentum(index, momentum)`. Reads return values; setting a
momentum does not update the record's mass metadata. Pickling the value preserves
both raw components and retained known mass, while record pickles keep the
existing array representation and do not acquire an implicit mass-shell claim.
