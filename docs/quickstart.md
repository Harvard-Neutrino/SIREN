# A complete data-free authoring workflow

Run the following blocks in order. The model supplies its own sampler;
DecayModel does not provide an implicit physical sampler. The source is a
1 GeV parent with momentum 3 GeV along z. Its decay length is 10 m, and the
sphere accepts decays between 10 and 12 m. The physical weighting includes
survival from the declared creation point.

```python
import math
import tempfile
from pathlib import Path

import siren


class IsotropicDecay(siren.DecayModel):
    parent = "N4"
    daughters = ("NuLight", "Gamma")
    measure = siren.Measure.SolidAngleRest()

    def total_width(self):
        return 3 * 1.973e-16 / 10

    def differential_width(self, record):
        return self.total_width() / (4 * math.pi)

    def sample(self, record, random):
        self.sample_isotropic(record, random)


def run(events=100):
    model = IsotropicDecay()
    volume = siren.geometry.Sphere(
        siren.geometry.Placement(siren.math.Vector3D(0, 0, 11)), 1.0, 0.0
    )
    with tempfile.TemporaryDirectory() as directory:
        source = Path(directory) / "source.csv"
        source.write_text(f"E,px,py,pz,x0,y0,z0,m\n{math.sqrt(10)},0,0,3,0,0,0,1\n")
        primary = siren.Vertex(
            "N4",
            model,
            distributions=[
                siren.dist.External(str(source)),
                siren.dist.BoundedPrimaryVertex(volume, 100),
            ],
            kinematics=siren.channels.physical(),
            weighting=siren.PropagatedFromCreation(),
        )
        simulation = siren.Simulation(
            detector=siren.detector.DetectorModel(), primary=primary,
            events=events, seed=31,
        )
        results = simulation.run(on_failure="raise", on_shortfall="raise")
        results.save(Path(directory) / "events", save_parquet=False)
        weights = results.weights
    probability = float(weights.sum())
    expected = math.exp(-1) * -math.expm1(-0.2)
    assert math.isclose(probability, expected, rel_tol=2e-12)
    return probability


```

Check the angular density and sampler at the specified kinematics with an
independent random stream, then generate and compare with the analytic survival
integral. The closure report is a finite diagnostic, not proof for other models,
energies, or arbitrarily fine angular structure.

```python
from siren._validation import build_template_record
model = IsotropicDecay()
record = build_template_record(model.GetPossibleSignatures()[0],
                               primary_mass=1.0, energy=math.sqrt(10))
report = siren.check_closure(model, record=record, samples=2000, seed=7)
assert report.ok, report
print(run())
```

The script `resources/examples/authoring_vertex.py` runs the same analytic
acceptance example. Its live Python model supports stored event export; it
cannot be archived as a native Results configuration. See [Results](results.md)
for the native snapshot/reweighting path, pooling and attempt normalization.
