# Three-body decay utilities

`from siren import three_body` provides the shared numerical machinery for
`M -> l + nu + phi`, where `nu` is massless. `l` is the spectator; `nu` and `phi`
form the pair in `Measure.Recursive2Body`. Masses and energies are in GeV, and
four-vectors are arrays ordered `(E, px, py, pz)`.

The functions accept the existing physics-object contract, without requiring a
base class: `m_M`, `m_l`, `m_phi`, `E_nu_max`, `_matel_sq(E_nu, E_phi)`, and
`_E_phi_limits(E_nu)`. The last method returns the allowed `phi` energy interval,
or `(None, None)` outside support. The matrix element is the spin-summed/averaged
squared amplitude in the model's convention. It must be nonnegative and
independent of orientation. The supported masses satisfy
`m_M > m_l + m_phi > 0`, with both daughter masses nonnegative.

| Function | Purpose |
| --- | --- |
| `dalitz_band(m_M, m_l, m_phi, E_nu)` | Allowed rest-frame `E_phi` interval; usable as `_E_phi_limits`. |
| `dalitz_width(physics)` | Integrate the specified channel's squared amplitude over the energy plane. |
| `find_max_weight(physics)` | Estimate a rejection envelope with the existing grid scan and 20% padding. |
| `sample_energies(physics, max_weight, random)` | Sample `(E_nu, E_phi)` from that amplitude. |
| `build_rest_momenta(physics, E_nu, E_phi, random)` | Return `(P_nu, P_l, P_phi)` at isotropic orientation. |
| `boost_four_vector(P_parent, P_rest, M_parent)` | Boost a daughter into the parent's lab frame. |
| `boost_to_rest(P_parent, P_lab)` | Apply the inverse pure boost without rotating the spatial basis. |
| `final_state_probability(physics, total_width, pdgid_nu, pdgid_phi, record)` | Normalized density in the specified Recursive2Body measure. |

For example, an existing physics object can be sampled directly:

```python
from siren import three_body

width = three_body.dalitz_width(physics)
bound = three_body.find_max_weight(physics)
energies = three_body.sample_energies(physics, bound, random)
p_nu, p_l, p_phi = three_body.build_rest_momenta(physics, *energies, random)
```

The width integral uses
`dGamma = |M|² dE_nu dE_phi / (64 pi³ m_M)`. The probability divides by the
supplied channel width and returns density with respect to
`ds_pair dOmega_pair dOmega_sub`. It expects an on-shell record with each of
the two distinct pair-daughter PDGs appearing once; their order in the record
may vary. A partial channel width is the appropriate denominator for conditional
sampling of that channel. Propagation lifetimes and channel selection remain
the model's responsibility.

The width integral requests relative accuracy of `1e-6` with no absolute
error floor, so small couplings do not relax the quadrature's error target.
Numerical integration error carries through to the normalized density.

The energy proposal samples `E_nu` uniformly, then `E_phi` uniformly inside its
band. Its acceptance weight is `|M|² * band_width`. The grid estimate is not a
proof that this weight is bounded; a model with narrow structure should supply
its own envelope. Invalid bounds raise `ValueError`. Invalid acceptance weights,
observed envelope overruns, and exhaustion of 10,000 proposals raise `RuntimeError`.
No substitute event is returned on rejection exhaustion.

These utilities do not provide a general three-massive-body sampler or a spin
density. BeamDecays' muon model adds its polarization-dependent orientation and
density factor explicitly. The generic closure gauge still reports three-body
coverage as incomplete; these functions do not change that certification.
