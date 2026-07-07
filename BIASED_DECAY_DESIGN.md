# BiasedDecay Adapter: Design and Worked Example

## Problem Statement

SIREN generates secondary particles (BSM decay products, scattered
particles) with kinematics sampled from the physical distribution.
For a distant detector, the fraction of secondaries that point toward
the detector is geometrically suppressed by the detector solid angle
(~10^-4 for CCM). This makes unbiased sampling extremely inefficient.

**Goal**: Force secondary particles toward the detector and weight
correctly, without modifying the physics model classes.

## Architecture: Decorator Pattern on Decay

```
                       +--------------------+
                       |    DirectionBias   |
                       | (detector geometry)|
                       +---------+----------+
                                 |
+---------------+    +-----------+-----------+
| ThreeBodyDecay|    |     BiasedDecay       |
| (physics)     |--->| implements Decay iface|
|               |    | wraps a physical Decay|
+---------------+    +-----------------------+
                                 |
               same Decay interface:
               SampleFinalState()
               FinalStateProbability()
               TotalDecayWidth()  <- delegates
               TotalDecayLength() <- delegates
```

- **Injection process** uses `BiasedDecay` (biased sampling + biased
  generation probability)
- **Physical process** uses the original `ThreeBodyDecay` (physical
  probability)
- **Weighter** computes `physical / generation` as usual, no changes
  to the weighting pipeline

## Worked Example: pi+ -> nu + mu + phi

### Setup

Parent: pi+ with lab-frame 4-momentum (E_pi, **p**_pi).
Daughters: neutrino (index 0), muon (index 1), scalar mediator phi
(index 2). We want phi to point toward the detector.

Rest-frame phase space (5 DOF):
1. E_phi* -- energy of phi in the pion rest frame (determines M_12^2
   of the nu+mu system)
2. cos(theta_phi*), phi_phi* -- direction of phi in the pion rest
   frame
3. cos(theta_mu^12), phi_mu^12 -- orientation of the muon in the
   nu+mu rest frame

### Unbiased Sampling (what ThreeBodyDecay currently does)

```python
# 1. Sample E_phi* from dGamma/dE_phi* via inverse CDF
E_phi_rest = sample_E_phi_rest(primary_type, lepton_type, E_meson)

# 2. Sample phi direction isotropically in rest frame
cos_theta_phi_rest = uniform(-1, 1)
phi_phi = uniform(0, 2*pi)

# 3. Construct phi rest-frame 3-momentum
p_phi_rest_mag = sqrt(E_phi_rest^2 - m_phi^2)
p_phi_rest = p_phi_rest_mag * (sin*cos, sin*sin, cos)

# 4. Construct nu+mu system, sample internal angles
M_12 = sqrt(m_pi^2 + m_phi^2 - 2*m_pi*E_phi_rest)
cos_theta_star = uniform(-1, 1)     # mu angle in nu+mu rest frame
phi_star = uniform(0, 2*pi)
# ... build lepton/neutrino momenta in M_12 rest frame ...

# 5. Boost nu+mu to pion rest frame, then everything to lab
```

### Biased Sampling (what BiasedDecay does)

The key insight: we reparameterize the phi kinematics from rest-frame
variables (E_phi*, cos_theta_phi*, phi_phi*) to lab-frame variables
(E_phi_lab, theta_phi_lab, phi_phi_lab). Given lab-frame variables,
the rest-frame variables are recovered by inverse boost.

```python
# 1. Sample lab-frame DIRECTION of phi toward detector
#    g(Omega_lab) = uniform over detector solid angle
cos_theta_phi_lab = sample_toward_detector()
phi_phi_lab = sample_azimuth_toward_detector()

# 2. Sample lab-frame ENERGY of phi uniformly in kinematic range
#    E_phi_lab in [E_min(theta_lab), E_max(theta_lab)]
E_phi_lab = uniform(E_min, E_max)

# 3. Construct lab-frame phi 4-momentum
p_phi_lab_mag = sqrt(E_phi_lab^2 - m_phi^2)
p_phi_lab = p_phi_lab_mag * direction_from(theta_phi_lab, phi_phi_lab)

# 4. Inverse boost to pion rest frame
p_phi_rest_4 = lorentz_boost([E_phi_lab, *p_phi_lab], -beta_pi)
E_phi_rest = p_phi_rest_4[0]

# 5. Check kinematic bounds on E_phi_rest
#    If out of bounds, reject and resample

# 6. Construct nu+mu system from E_phi_rest (same as unbiased)
M_12 = sqrt(m_pi^2 + m_phi^2 - 2*m_pi*E_phi_rest)
cos_theta_star = uniform(-1, 1)
phi_star = uniform(0, 2*pi)
# ... build lepton/neutrino momenta, boost to lab ...
```

Steps 4-6 are identical to the unbiased code. The BiasedDecay only
replaces how phi's direction and energy are chosen (steps 1-3),
then delegates the nu+mu subsystem sampling to the wrapped decay
(or reimplements it identically, since it is model-independent
2-body kinematics within the M_12 system).

### The Generation Probability

The physical `DifferentialDecayWidth` is (from the student's code):

```
dGamma / (dE_phi_rest * dOmega_rest) = dGamma/dE_phi_rest / (4*pi)
```

where `dGamma/dE_phi_rest` encodes the matrix element integrated over
the nu+mu subsystem angles. The `/ (4*pi)` factor comes from the
isotropic rest-frame angular sampling.

The biased `FinalStateProbability` must return the generation density
**in the same measure** as the physical one, so the ratio cancels
correctly.

We sampled 3 variables:
- Omega_phi_lab from g(Omega_lab) = 1 / Delta_Omega_det
- E_phi_lab from h(E_lab) = 1 / (E_max - E_min)
- The nu+mu subsystem angles (cos_theta_star, phi_star) from 1/(4*pi)

The generation density in the original rest-frame measure
(E_phi_rest, Omega_rest, Omega_star) is:

```
p_gen = g(Omega_lab) * h(E_lab) * |J|^{-1} * (1 / 4*pi)
```

where |J| is the Jacobian determinant of the transformation from
rest-frame variables to lab-frame variables:

```
|J| = |d(E_phi_lab, Omega_phi_lab) / d(E_phi_rest, Omega_phi_rest)|
```

This Jacobian factors as:

```
|J| = |dE_lab/dE_rest| * |dOmega_lab/dOmega_rest|
```

because the energy and direction transforms are coupled through the
boost. The individual factors (from the student's existing code,
lines 1875-1876):

```
dE_rest/dE_lab = 1 / (gamma * (1 + beta/beta_phi * cos_theta_rest))
dOmega_rest/dOmega_lab = gamma^2 * (1 + beta * cos_theta_rest)^2
```

So:

```
|J| = [gamma * (1 + beta/beta_phi * cos_theta_rest)]
    * [1 / (gamma^2 * (1 + beta * cos_theta_rest)^2)]
    = (1 + beta/beta_phi * cos_theta_rest)
    / (gamma * (1 + beta * cos_theta_rest)^2)
```

And the biased `FinalStateProbability` returns:

```
p_gen = [1 / Delta_Omega_det]
      * [1 / (E_max - E_min)]
      * |J|^{-1}
      * [1 / (4*pi)]
```

The nu+mu subsystem factor (1/4*pi) cancels between physical and
generation since both sample it identically. So effectively:

```
weight = [dGamma/dE_phi_rest / (4*pi)]
       / [1/(Delta_Omega_det * (E_max - E_min)) * |J|^{-1} * 1/(4*pi)]

       = dGamma/dE_phi_rest * Delta_Omega_det * (E_max - E_min) * |J|
```

### What BiasedDecay needs from ThreeBodyDecay

Only `dGamma/dE_phi_rest` -- the differential decay width as a
function of the rest-frame energy. This is what the student's
`differential_decay_width_rest` method already computes.

BiasedDecay does NOT need:
- The matrix element parameterization
- The Dalitz plot variables
- The coupling constants
- The form factors

It only needs the 1D energy spectrum in the rest frame, which is
model-specific but already exposed as a method on the decay class.

### What BiasedDecay computes itself (model-independent)

- Detector solid angle and direction sampling (from DirectionBias)
- Lab-frame energy bounds for a given lab direction
- The inverse boost: lab 4-momentum -> rest-frame 4-momentum
- The Jacobian |d(E_lab, Omega_lab) / d(E_rest, Omega_rest)|
- The nu+mu subsystem kinematics (pure 2-body, no model dependence)

All of these are functions of masses and the parent boost only.

## Comparison with the Student's Approach

The student tried to implement biasing **inside** ThreeBodyDecay.
This led to three problems:

1. **The biased SampleFinalState is incomplete.** The `if biased_inj:`
   branch samples lab-frame angles but never computes the rest-frame
   energy or direction. The code falls through to variables defined
   only in the `else:` branch, causing `UnboundLocalError`.

2. **The biased DifferentialDecayWidth returns the wrong thing.** It
   returns `total_width * cone_pdf * jacobian`, which is a generation
   probability estimate, not a differential width. The physical
   `DifferentialDecayWidth` returns `dGamma/(dE dOmega)` in the lab
   frame. These are in different units and different measures. The
   Weighter divides them expecting the same measure.

3. **Frame confusion in the Jacobians.** The unbiased path computes
   `dGamma/(dE_rest dOmega_rest)` then multiplies by lab-frame
   Jacobians. The biased path tries to express the cone density
   directly in lab-frame solid angle. These are inconsistent
   parameterizations within the same function.

The BiasedDecay adapter avoids all three problems:
- Sampling is complete because it follows a clear algorithm
  (sample lab direction + energy, inverse boost, delegate subsystem)
- FinalStateProbability returns the generation density in the same
  measure as the physical decay (rest-frame E_phi, Omega, Omega_star)
- Frame conversions are handled by one Jacobian computation that
  connects lab sampling to rest-frame measure

## Generalization

### 2-Body Decay (P -> A + B)

Degenerate case: E_phi_rest is fixed (no energy DOF). The biased
sampling is just:
1. Sample Omega_lab from g (toward detector)
2. Solve for cos_theta_rest from the boost equations (0, 1, or 2
   solutions)
3. Jacobian is the 2x2 |dOmega_lab / dOmega_rest| (MeVPrtl Eq. 6)

No energy sampling, no inverse CDF, no subsystem.

### Cross-Section Biasing (chi + N -> chi + gamma + N)

The scattering kinematics are already in the lab frame. Direction
biasing means biasing the outgoing particle angle (cos_theta_cm).
The student's `MediatorScattering` class already does this correctly
with `biased_inj=True`. A `BiasedCrossSection` adapter would follow
the same pattern as `BiasedDecay`.

### Arbitrary N-Body Decay

For N-body decays with N > 3, the same factorization applies:
1. Sample the lab direction of the target daughter (2 DOF)
2. Sample the lab energy of the target daughter (1 DOF)
3. Inverse boost to rest frame
4. Sample the remaining subsystem kinematics

The subsystem grows (from 2-body for 3-body decays, to 3-body for
4-body decays, etc.) but the biasing logic is the same.

## Implementation Plan

### C++ Classes

```
DirectionBias (new)
  - SampleDirection(random) -> (theta_lab, phi_lab)
  - Density(theta_lab, phi_lab) -> double
  - Implementations: ConeDirectionBias, BoxDirectionBias

BiasedDecay : public Decay (new)
  - Wraps shared_ptr<Decay> physical_decay_
  - Wraps shared_ptr<DirectionBias> bias_
  - int target_daughter_index_
  - SampleFinalState: biased sampling algorithm above
  - FinalStateProbability: g * h * |J|^-1 * subsystem_density
  - TotalDecayWidth: delegates to physical_decay_
  - TotalDecayLength: delegates to physical_decay_

BiasedCrossSection : public CrossSection (new, same pattern)
```

### Python Interface

```python
import siren

# Physical decay (the student's existing class)
physical_decay = ThreeBodyDecay(mediator_mass=0.01, ...)

# Direction bias toward CCM detector
bias = siren.interactions.ConeDirectionBias(
    axis=[0.9997, 0.0, -0.026],     # detector direction
    half_angle=0.075,                # cone half-angle (rad)
)

# Biased adapter
biased_decay = siren.interactions.BiasedDecay(
    physical_decay,
    bias,
    target_daughter=siren.particles.PhiPrime,  # which daughter to bias
)

# Use in injection process (physical process uses physical_decay)
sim = siren.Simulation(
    ...
    interactions=bundle.primary[siren.particles.NuMu],
    secondary_interactions={
        siren.particles.PhiPrime: [biased_decay],  # injection
    },
    ...
)
# Physical process automatically uses the unbiased decay for weighting
```

### Validation

The biased and unbiased simulations must produce statistically
identical weighted distributions. For each observable O:

```
<O>_biased = sum(w_i * O_i) / sum(w_i)  [biased events]
<O>_unbiased = sum(O_i) / N             [unbiased events, unit weight]
```

These must agree within statistical uncertainty. The biased version
should have lower variance (that is the point of biasing).
