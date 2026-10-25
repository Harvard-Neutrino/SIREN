"""Synthetic source-table demonstration; not a beam-rate prediction."""
import argparse
import math
from pathlib import Path

import numpy as np
import siren
from siren import dk2nu


def run(events=256, output=None):
    # Production energy is retained separately from decay-point momentum.
    data = {"ptype": np.array([111, 111]), "E": np.array([9., 9.]),
            "px": np.zeros(2), "py": np.zeros(2), "pz": np.array([1., 2.]),
            "vx": np.zeros(2), "vy": np.zeros(2), "vz": np.array([0., 100.]),
            "t0": np.array([3., 7.]), "nimpwt": np.array([0.2, 0.8]), "pot": 10.,
            "source_file_id": np.array([123, 123]), "source_entry": np.array([0, 1]),
            "source_sampling_probability": np.ones(2)}
    parents = dk2nu.dk2nu_to_primary_distribution(
        data, None, frame="detector", sampling_bias=data["nimpwt"])
    signature = siren.dataclasses.InteractionSignature()
    signature.primary_type = siren.particles.Pi0
    signature.target_type = siren.particles.Decay
    signature.secondary_types = [siren.particles.Gamma] * 2
    decay = siren.injection.PhaseSpaceDecay(
        signature, [0., 0.], 1e-16, 1e-16, siren.injection.Isotropic2BodyChannel())
    primary = siren.Vertex("Pi0", decay, distributions=[parents], physical=[parents],
                           kinematics=siren.channels.physical(), weighting=siren.Fixed())
    simulation = siren.Simulation(detector=siren.detector.DetectorModel(),
                                  primary=primary, events=events, seed=317)
    results = simulation.run(on_failure="raise", on_shortfall="raise")
    assert math.isclose(float(results.weights.sum()), 0.1, rel_tol=2e-12)
    if output is not None:
        results.save(Path(output), save_parquet=False)
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events", type=int, default=256)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run(args.events, args.output)
    print(f"Synthetic primaries per exposure: {sum(result.weights):.12g}")
