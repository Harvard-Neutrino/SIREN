# Python interaction models

`DecayModel` and `CrossSectionModel` provide signatures, topology, normalized
physical densities, and identity-based equality. Declare the particles and
measure, then implement the total rate, differential rate, and `sample`.
Override `equal(other)` when two distinct instances should compare by value.
The same rules apply to bases created with `decay_model_base(base=...)` and
`cross_section_model_base(base=...)`.

Selected-signature weighting uses the declared measure, including its daughter
indices, for every signature the model advertises through
`GetPossibleSignatures()`; the topology follows each signature's final-state
count. Signatures outside that set resolve to `Unspecified`.
`DensityVariables()` labels do not override the declaration.

The differential hook receives an `InteractionRecord`. The sampler receives a
`CrossSectionDistributionRecord` and writes its secondary particle records.
Sampling must follow the normalized differential rate in the declared measure;
the measure alone does not determine the distribution.

For an isotropic two-body decay, the existing engine sampler is available
through `sample_isotropic`:

```python
import math
import siren


class IsotropicDecay(siren.DecayModel):
    parent = "N4"
    daughters = ("NuLight", "Gamma")
    measure = siren.Measure.SolidAngleRest()

    def total_width(self):
        return 1e-12  # GeV

    def differential_width(self, record):
        return self.total_width() / (4 * math.pi)

    def sample(self, record, random):
        self.sample_isotropic(record, random)
```

`sample_isotropic` preserves configured secondary masses and fills missing
storage from `SecondaryMasses`, including for the legacy Injector interface.
It is suitable only when the physical density is uniform in rest-frame solid
angle. An anisotropic model supplies its own sampler. Direct overrides of
`SampleFinalState` remain supported.

Isotropic two-body decay models that previously relied on implicit sampling
need the explicit `sample` method shown above. Their rate and density functions
need no changes. Cross-section models, including those created with
`cross_section_model_base`, must implement their own scattering sampler;
`sample_isotropic` supplies decay kinematics and is available only on decay bases.
The override audit catches a missing sampler and a missing `equal` in direct
subclasses of the abstract native interaction bases.

## Checking closure

```python
report = siren.check_closure(model, record=initial_record, samples=8000, seed=7)
print(report)
report.raise_if_failed()
```

Supply `record` to test a particular signature, masses, and initial momentum.
Secondary output arrays may be empty; missing masses are resolved through
`SecondaryMasses`, while supplied masses are preserved. The record is left
unchanged. The alternative `primary_energy` and `target` shortcuts
construct a synthetic template; unknown masses receive synthetic defaults.
Each call uses independent validation RNG streams and recomputes the result,
so changing model parameters cannot reuse a stale result.

The built-in reference covers two-body `SolidAngleRest` decays. It checks the
absolute density integral and the joint angular distribution using independent
reference samples. The report names measured coordinates and includes reference
Monte Carlo uncertainty. It tests both individual bins and their combined
deviation, using the full covariance of the model and weighted reference
histograms. This retains sensitivity to correlations spread across many bins.
`report.joint_shape` contains the joint chi-square statistic, degrees of freedom,
and p-value. The joint test uses the two-sided Gaussian tail probability
corresponding to `tol_sigma`. Passing applies to the tested kinematics and histogram
resolution, not to arbitrarily fine structure or untested model parameters.

The density must be evaluable from the reference kinematics and configured
template parameters. If it depends on a coordinate cached in
`interaction_parameters` by the sampler, closure reports incomplete coverage:
the independent reference cannot reconstruct that cache. Compute such density
coordinates from the momenta so both paths evaluate the same physical function.
Sampler-written bookkeeping that does not affect the density is permitted.

`report.checks` distinguishes `passed`, `failed`, and `incomplete` checks.
`report.ok` requires every check to pass. `raise_if_failed()` also rejects
incomplete coverage. Unsupported measures, unresolved angles, sparse reference
coverage, and mixture density probes that have not run cannot certify closure.
Mixture configuration validation is reported separately. Three-body and
scattering models still need an independent reference before full certification.
