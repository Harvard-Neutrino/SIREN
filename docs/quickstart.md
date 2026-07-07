# Quickstart

SIREN is a two-phase Monte Carlo framework: injection draws biased samples of
rare event topologies, and weighting computes the importance weight that maps
each sample back to the physical distribution it represents. This page is the
smallest end-to-end taste of both phases: define a model, verify its sampler
and its density agree, then see the shape of a full simulation run.

## Define a model

A physics model is a small Python class: declare the particles and the
phase-space measure, then implement the width (or cross section) as two
hooks. Here is a two-body decay with an isotropic angular distribution in the
parent rest frame:

```python
import math, siren

class IsoDecay(siren.DecayModel):
    parent = "N4"
    daughters = ("NuLight", "Gamma")
    measure = siren.Measure.SolidAngleRest()

    def total_width(self):
        return 1.0

    def differential_width(self, record):
        return 1.0 / (4.0 * math.pi)

    def density_variables(self):
        return "cost"
```

`siren.DecayModel` derives the interaction signature, the topology, and
`FinalStateProbability` (differential over total) from these hooks. Because
the declared measure (`SolidAngleRest`) has a self-contained default sampler,
you get a working `SampleFinalState` for free -- you only had to describe the
physics.

## Check it closes

Every weight in SIREN is `PhysicalProbability / GenerationProbability`. That
ratio is only meaningful if the sampler's distribution and the model's own
density (`FinalStateProbability`) describe the same thing. `check_closure`
quantifies whether they do:

```python
report = siren.check_closure(IsoDecay())
print(report)
print(report.ok)
```

Actual output from this exact model:

```
Closure report: PASS
  normalization (absolute) E[f/g_ref] = 1.0000 +/- 0.0000 (expect 1.0)
  flatness (shape only) = 1.0000 +/- 0.0224 (self-normalized bin ratios; blind to a uniform scale error)
  moment z-scores:
    cost: z=0.00
  worst region: costheta_secondary0_rest in [-0.40, -0.30): ratio 1.19 (z=1.7)
True
```

Closure means the biased sampler and the physical density agree, so weights
computed against this model are unbiased -- no hidden shape or normalization
error that statistics alone would never reveal.

## Run a simulation

The `Simulation` facade wires a detector, a primary particle, an interaction
model, and injection/physical distributions into one object that runs both
phases:

```python
import siren

sim = siren.Simulation(
    events=100_000,
    detector="IceCube",
    primary="NuMu",
    interactions="CSMSDISSplines",
    target="Nucleon",
    process="CC",
    energy=siren.dist.PowerLaw(2, 1e3, 1e6),
    direction=siren.dist.IsotropicDirection(),
    position=siren.dist.ColumnDepth(
        600, 600.0, siren.distributions.LeptonDepthFunction()),
)
results = sim.run()
len(results)                       # number of generated events
for event, weight in results:
    pass                          # event: InteractionTree, weight: float
sim.describe()                     # human-readable summary of the pipeline
```

This snippet needs a detector model and cross-section splines to run (both
loaded from `resources/`), so it is shown here but not executed on this page.
See `docs/quickstart.md`'s companion tests (`tests/python/test_simulation_api.py`)
for a runnable version against the `IceCube` detector and `CSMSDISSplines`.

`results` is a `siren.Results`: iterate it for `(event, weight)` pairs, index
or slice it, and call `results.summary()` for a weight-quality report.

## Which layer do I use?

Use the `Simulation` facade for standard runs. Drop to the Layer-2
`Injector`/`Weighter` objects (`sim.injector`, `sim.weighter`) only when you
need to inspect or drive generation directly -- each `run()` call builds its
own `Injector`/`Weighter` from the `Simulation`'s stored configuration, so
`sim.injector`/`sim.weighter` only reflect the objects a call to `run()` used
once that call has returned; accessing them beforehand does not let you
customize the objects `run()` will build.

## Next steps

- `docs/units_conventions.md` -- units, geometry argument conventions, frames,
  and the RNG reproducibility contract.
- `docs/errors.md` -- what each SIREN error message means and how to fix it.
