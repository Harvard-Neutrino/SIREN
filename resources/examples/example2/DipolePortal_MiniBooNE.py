"""Dipole-portal heavy neutral lepton (N4) at MiniBooNE in the BNB beam.

BNB muon neutrinos upscatter to N4 through the DarkNews dipole portal
(primary vertex) and the N4 decays to a photon inside the tank's signal
region (secondary vertex). Injection samples the unit-normalized flux shape
and places the primary vertex along the N4 decay range in front of the
detector; weighting uses the physically normalized flux. The primary's
``expand`` rule recurses only into the N4, whose vertex is terminal.
"""
import argparse

import siren

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--events", type=int, default=100000)
parser.add_argument("--seed", type=int, default=None)
parser.add_argument("--output", default=None,
                    help="output prefix (default: output/MiniBooNE_Dipole_<params>_example)")
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

# MiniBooNE is provided by the SBN composite model, which places the tank
# at its surveyed location (0, 1.896, 541.34) m in the BNB beam frame with
# the proton target at the origin. Detector coordinates, in which every
# distribution below is expressed, have their origin at the tank center;
# the decay-range bounds (3 m to 541 m) are distances upstream of it.
detector_model = siren.load_detector("SBN", detector="MiniBooNE")


def _sector_sphere(detector_model, name):
    """A spherical sector of the model re-placed in detector coordinates."""
    geo = next(s.geo for s in detector_model.Sectors if s.name == name)
    center = detector_model.GeoPositionToDetPosition(
        siren.detector.GeometryPosition(geo.placement.Position)).get()
    return siren.geometry.Sphere(
        siren.geometry.Placement(center, geo.placement.Quaternion),
        geo.Radius, geo.InnerRadius)


# The inner signal region of the tank (the oil inside the optical barrier),
# used to bound the secondary decay vertex.
miniboone_signal_region = _sector_sphere(detector_model, "vol_mb_inner_oil")

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
            "BNB", tag="FHC_numu",
            min_energy=model_kwargs["m4"], max_energy=10,
            physically_normalized=False),
        direction,
        siren.dist.RangePosition(
            6.2, 6.2,
            siren.dist.DecayRange(
                model_kwargs["m4"],
                min(_total_decay_width(d) for d in hnl_decays),
                3, 541,
            ),
            siren.get_detector_model_targets(detector_model),
        ),
    ],
    physical=[
        siren.load_flux("BNB", tag="FHC_numu", physically_normalized=True),
        direction,
    ],
    expand=(siren.expand.child("N4"),),
)

hnl = siren.Vertex(
    N4, hnl_decays,
    position=siren.dist.BoundedVertex(miniboone_signal_region, 25),
    expand=(siren.expand.depth_below(0),),  # terminal: photon and neutrino are final
)

sim = siren.Simulation(
    events=args.events, seed=args.seed, detector=detector_model,
    primary=primary, secondaries=[hnl])
results = sim.run()
results.summary()
results.save(
    args.output or "output/MiniBooNE_Dipole_M%2.2e_mu%2.2e_example"
    % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"]))
