"""Dipole-portal heavy neutral lepton (N4) at MINERvA in the NuMI beam.

NuMI muon neutrinos upscatter to N4 through the DarkNews dipole portal
(primary vertex) and the N4 decays to a photon in the fiducial volume
(secondary vertex). Injection samples the unit-normalized flux shape and
places the primary vertex along the N4 decay range in front of the detector;
weighting uses the physically normalized flux. The primary's ``expand`` rule
recurses only into the N4, whose vertex is terminal.
"""
import argparse

import siren

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--events", type=int, default=100000)
parser.add_argument("--seed", type=int, default=None)
parser.add_argument("--output", default=None,
                    help="output prefix (default: output/MINERvA_Dipole_<params>_example)")
args = parser.parse_args()

model_kwargs = {
    "m4": 0.47,
    "mu_tr_mu4": 2.50e-6,
    "UD4": 0,
    "Umu4": 0,
    "epsilon": 0.0,
    "gD": 0.0,
    "decay_product": "photon",
    "noHC": True,
    "HNLtype": "dirac",
}

NuMu, N4 = siren.particles.NuMu, siren.particles.N4

detector_model = siren.load_detector("MINERvA")

dn_version = siren.utilities.darknews_version()
table_name = f"DarkNewsTables-v{dn_version}/"
table_name += "Dipole_M%2.2e_mu%2.2e" % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"])

bundle = siren.load_processes(
    "DarkNewsTables",
    primary_type=NuMu,
    detector_model=detector_model,
    table_name=table_name,
    **model_kwargs,
)
hnl_decays = bundle.secondary[N4]


def _total_decay_width(decay):
    # TotalDecayWidth(record) only returns a nonzero width when the record's
    # signature matches the decay's own signature; a default-constructed
    # InteractionRecord has an "unknown" signature and always yields zero.
    record = siren.dataclasses.InteractionRecord()
    record.signature = decay.GetPossibleSignatures()[0]
    return decay.TotalDecayWidth(record)


direction = siren.dist.FixedDirection([0, 0, 1])

primary = siren.Vertex(
    NuMu, bundle.primary[NuMu],
    distributions=[
        siren.dist.PrimaryMass(0),
        siren.load_flux(
            "NUMI", tag="FHC_ME_numu",
            min_energy=model_kwargs["m4"], max_energy=20,
            physically_normalized=False),
        direction,
        siren.dist.RangePosition(
            1.24, 5.0,
            siren.dist.DecayRange(
                model_kwargs["m4"],
                min(_total_decay_width(d) for d in hnl_decays),
                3, 240,
            ),
            siren.get_detector_model_targets(detector_model),
        ),
    ],
    physical=[
        siren.load_flux("NUMI", tag="FHC_ME_numu", physically_normalized=True),
        direction,
    ],
    expand=(siren.expand.child("N4"),),
)

hnl = siren.Vertex(
    N4, hnl_decays,
    position=siren.dist.BoundedVertex(siren.get_fiducial_volume("MINERvA"), 25),
    expand=(siren.expand.depth_below(0),),  # terminal: photon and neutrino are final
)

sim = siren.Simulation(
    events=args.events, seed=args.seed, detector=detector_model,
    primary=primary, secondaries=[hnl])
results = sim.run()
results.summary()
results.save(
    args.output or "output/MINERvA_Dipole_M%2.2e_mu%2.2e_example"
    % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"]))
