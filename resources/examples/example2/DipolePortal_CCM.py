import os
import numpy as np
import siren

model_kwargs = {
    "m4": 0.0235,
    "mu_tr_mu4": 6e-7,
    "UD4": 0,
    "Umu4": 0,
    "epsilon": 0.0,
    "gD": 0.0,
    "decay_product": "photon",
    "noHC": True,
    "HNLtype": "dirac",
}

detector_model = siren.load_detector("CCM")
fiducial = siren.get_fiducial_volume("CCM")

dn_version = siren.utilities.darknews_version()
table_name = f"DarkNewsTables-v{dn_version}/"
table_name += "Dipole_M%2.2e_mu%2.2e" % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"])

bundle = siren.load_processes(
    "DarkNewsTables",
    primary_type=siren.particles.NuMu,
    detector_model=detector_model,
    table_name=table_name,
    **model_kwargs,
)

target_origin = siren.math.Vector3D(0, 0, -0.241)
detector_origin = siren.math.Vector3D(23, 0, -0.65)
beam_dir = detector_origin - target_origin
beam_dir.normalize()

sim = siren.Simulation(
    n_events=1,
    detector=detector_model,
    primary=siren.particles.NuMu,
    interactions=bundle.primary[siren.particles.NuMu],
    energy=siren.dist.Monoenergetic(0.02965),
    injection_direction=siren.dist.Cone(beam_dir, np.arctan(5 / 23.0)),
    physical_direction=siren.dist.IsotropicDirection(),
    position=siren.dist.PointSource(target_origin - detector_origin, 25),
    secondary_interactions=bundle.secondary,
    secondary_position=siren.dist.BoundedVertex(fiducial, 25),
    stopping_condition=lambda datum, i: (
        datum.record.signature.secondary_types[i] != siren.particles.N4
    ),
)

results = sim.run()

os.makedirs("output", exist_ok=True)
results.save("output/CCM_Dipole")
