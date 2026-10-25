"""Regenerate with the unmodified U05 53c64638 runtime, before U06 is installed."""
import argparse
from pathlib import Path
import siren


def generate(directory):
    directory.mkdir(parents=True, exist_ok=True)
    for label, weights in (("positive", [2., 2.]), ("negative", [-1., 2.])):
        csv = directory / f"{label}.csv"
        csv.write_text("E,px,py,pz,m,x0,y0,z0,weight\n" + "".join(
            f"2,0,0,1.7320508075688772,1,0,0,0,{w:g}\n" for w in weights))
        model = siren.interactions.HNLDipoleDecay(
            1., 1e-6, siren.interactions.HNLDipoleDecay.ChiralNature.Dirac)
        process = siren.injection.PrimaryInjectionProcess(
            siren.particles.N4, siren.interactions.InteractionCollection(siren.particles.N4, [model]))
        process.weighting_mode = siren.Fixed()
        process.distributions = [siren.distributions.PrimaryExternalDistribution(str(csv)),
            siren.distributions.SphereVolumePositionDistribution(siren.geometry.Sphere(1., 0.))]
        injector = siren.injection._Injector(
            1, siren.detector.DetectorModel(), process, siren.utilities.SIREN_random(1))
        injector.SetRandom(None)
        injector.SaveInjector(str(directory / f"{label}.injector"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Fresh output directory")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Use a fresh directory to preserve existing fixtures")
    generate(args.output)
