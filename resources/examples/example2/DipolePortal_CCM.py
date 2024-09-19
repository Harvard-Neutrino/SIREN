import os
import numpy as np
import siren
from siren import utilities

# Define a DarkNews model
model_kwargs = {
    "m4": 0.0235,
    "mu_tr_mu4": 6e-7,  # GeV^-1
    "UD4": 0,
    "Umu4": 0,
    "epsilon": 0.0,
    "gD": 0.0,
    "decay_product": "photon",
    "noHC": True,
    "HNLtype": "dirac",
}

# Number of events to inject
events_to_inject = 1

# Experiment to run
experiment = "CCM"
detector_model = utilities.load_detector(experiment)

# Particle to inject
primary_type = siren.dataclasses.Particle.ParticleType.NuMu

# Load DarkNews processes
primary_processes, secondary_processes = utilities.load_processes(
    f"DarkNewsTables-v{siren.utilities.darknews_version()}",
    primary_type=primary_type,
    detector_model = detector_model,
    **model_kwargs,
)

# Mass distribution
mass_ddist = siren.distributions.PrimaryMass(0)

# Primary distributions
nu_energy = 0.02965  # from pi+ DAR
edist = siren.distributions.Monoenergetic(nu_energy)
flux_units = siren.distributions.NormalizationConstant(3.76e-2)

# Cone direction distribution
opening_angle = np.arctan(5 / 23.0)
lower_target_origin = siren.math.Vector3D(0, 0, -0.241)
detector_origin = siren.math.Vector3D(23, 0, -0.65)
lower_dir = detector_origin - lower_target_origin
lower_dir.normalize()
lower_inj_ddist = siren.distributions.Cone(lower_dir, opening_angle)
phys_ddist = siren.distributions.IsotropicDirection()

# Position distribution
max_dist = 25
lower_pos_dist = siren.distributions.PointSourcePositionDistribution(
    lower_target_origin - detector_origin, max_dist,
)

primary_injection_distributions = [
    mass_ddist,  # Mass distribution
    edist,  # Energy distribution
    lower_inj_ddist,  # Direction distribution
    lower_pos_dist  # Position distribution
]

primary_physical_distributions = [
    edist,  # Energy distribution
    phys_ddist,  # Direction distribution
]

fiducial_volume = siren.utilities.get_fiducial_volume(experiment)
secondary_injection_distributions = {}
for secondary_type in secondary_processes.keys():
    secondary_injection_distributions[secondary_type] = [
        siren.distributions.SecondaryBoundedVertexDistribution(fiducial_volume, max_dist)
    ]

# Define stopping condition for the injector
def stop(datum, i):
    secondary_type = datum.record.signature.secondary_types[i]
    return secondary_type != siren.dataclasses.Particle.ParticleType.N4

injector = siren.injection.Injector()
injector.number_of_events = 1
injector.detector_model = detector_model
injector.primary_type = primary_type
injector.primary_interactions = primary_processes[primary_type]
injector.primary_injection_distributions = primary_injection_distributions
injector.secondary_interactions = secondary_processes
injector.secondary_injection_distributions = secondary_injection_distributions
injector.stopping_condition = stop


# Generate events
events = [injector.generate_event() for _ in range(events_to_inject)]

# Output the events
os.makedirs("output", exist_ok=True)

weighter = siren.injection.Weighter()
weighter.injectors = [injector]
weighter.detector_model = detector_model
weighter.primary_type = primary_type
weighter.primary_interactions = primary_processes[primary_type]
weighter.secondary_interactions = secondary_processes
weighter.primary_physical_distributions = primary_physical_distributions
weighter.secondary_physical_distributions = {}

weights = [weighter(event) for event in events]
