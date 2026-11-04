"""Data-free two-body authoring example with an analytic survival check."""

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
        generator = siren.injection.Injector(
            detector=siren.detector.DetectorModel(),
            primary=primary,
            events=events,
            seed=31,
        )
        trees = generator.generate(on_failure="raise", on_shortfall="raise")
        weights = siren.injection.Weighter(generator).weight_all(trees)
    probability = float(weights.sum())
    expected = math.exp(-1) * -math.expm1(-0.2)
    assert math.isclose(probability, expected, rel_tol=2e-12)
    return probability


if __name__ == "__main__":
    print(f"Decay probability between 10 and 12 m: {run():.12g}")
