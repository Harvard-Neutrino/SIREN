import os
import siren
from siren import utilities

# Number of events to inject
events_to_inject = int(1e5)

# Experiment to run
experiment = "DUNEFD"
detector_model = utilities.load_detector(experiment)
primary_type = siren.dataclasses.Particle.ParticleType.NuMu
target_type = siren.dataclasses.Particle.ParticleType.Nucleon

# Primary interactions and distributions
primary_processes, _ = utilities.load_processes(
    "CSMSDISSplines",          # model_name
    primary_types=[primary_type],
    target_types=[target_type],
    isoscalar=True,             # for isoscalar splines
    process_types=["CC"]         # specify the process type, e.g., "CC" for charged current
)

# Energy distribution
edist = siren.distributions.PowerLaw(1, 1e3, 1e6)

# Direction distribution
direction_distribution = siren.distributions.IsotropicDirection()

# Position distribution
muon_range_func = siren.distributions.LeptonDepthFunction()
position_distribution = siren.distributions.ColumnDepthPositionDistribution(
    60, 60.0, muon_range_func)


# Define injection distributions
primary_injection_distributions = {
    "energy": edist,
    "direction": direction_distribution,
    "position": position_distribution
}

# Set up the Injector
injector = siren.injection.Injector()
injector.number_of_events = events_to_inject
injector.detector_model = detector_model
injector.primary_type = primary_type
injector.primary_interactions = primary_processes[primary_type]
injector.primary_injection_distributions = [
    siren.distributions.PrimaryMass(0),
    edist,
    direction_distribution,
    position_distribution
]

# Generate events
event = injector.generate_event()

# Set up the Weighter for event weighting
weighter = siren.injection.Weighter()
weighter.injectors = [injector]
weighter.detector_model = detector_model
weighter.primary_type = primary_type
weighter.primary_interactions = primary_processes[primary_type]
weighter.primary_physical_distributions = [
    edist,
    direction_distribution,
]

# Compute weight
weight = weighter(event)

# Output events and weights
os.makedirs("output", exist_ok=True)
print(str(event))
print(f"Event weight: {weight}")

# Save events
# TODO

